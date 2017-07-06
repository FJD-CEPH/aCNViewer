import os
try:
    import pp
except ImportError:
    # if not os.getenv('SILENT'):
    #    print 'pp modude not installed'
    # pass
    pp = None
import types
import threading
import time
import copy

from Utilities import iterRange, Utilities
from ValueParser import ValueParser
from ObjectBase import defaultdict, runFromTerminal
from FileHandler import FileNameGetter


class T:

    def run(self, cmd):
        os.system(cmd)


''' Function taken from parallel python '''


def _detect_ncpus():
    if os.environ.get('NB_CPUS'):
        return int(os.environ['NB_CPUS'])
    """Detects the number of effective CPUs in the system"""
    # for Linux, Unix and MacOS
    if hasattr(os, "sysconf"):
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux and Unix
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else:
            # MacOS X
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # for Windows
    if "NUMBER_OF_PROCESSORS" in os.environ:
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    # return the default value
    return 1


def _getNbAvailableCpus():
    fh = os.popen('uptime')
    resStr = fh.read()
    fh.close()
    if 'average' not in resStr:
        raise NotImplementedError(
            'uptime result does not contain "average" keyword: "%s"' % resStr)
    try:
        busyCpu = int(float(resStr.split(
            'average')[-1].split(':')[-1].strip().split(' ')[0].rstrip(',').
            replace(',', '.')))
    except:
        print 'Unexpected uptime format: "%s"' % resStr
        raise
    return max(_detect_ncpus() - busyCpu, 1)


class StoppableThread (threading.Thread):
    """Thread class with a stop() method. The thread itself has to check
    regularly for the stopped() condition."""

    def __init__(self):
        self._stop = threading.Event()
        super(StoppableThread, self).__init__()

    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()


class ThreadFunc(threading.Thread):

    def __init__(self, func, *args):
        self.__func = func
        self._args = args
        self._error = None
        threading.Thread.__init__(self)

    def run(self):
        try:
            self._res = self.__func(*self._args)
        except:
            import Utilities
            self._error = Utilities.getLastErrorMessage()
            raise


class ThreadManager:
    RUNNING = '|<*RUNNING*>|'
    TASK_ID_VAR = None

    def __init__(self, maxThreadNb=None, keepResults=False, timeout=None,
                 adjustProcNb=False):
        self._resList = []
        self.__threadList = []
        self.__maxThreadNb = maxThreadNb
        self.__keepResults = keepResults
        self.__timeout = timeout
        self.__nbProc = 1
        self.__adjustProcNb = adjustProcNb
        self.__canProcBeAdjusted = False
        self.__maxProcNb = _detect_ncpus()
        self.__jobDict = {}

    def __del__(self):
        print 'Waiting for all jobs to complete before exiting'
        self.wait()

    def getRunningJobNb(self):
        self.__removeFinishedJobs()
        return len(self.__threadList)

    getNbRunningJobs = getRunningJobNb

    def clear(self):
        self._resList = []

    def __iter__(self):
        if not self.__keepResults:
            raise NotImplementedError
        return iter(self._resList)

    def __removeFinishedJobs(self):
        newThreadList = []
        while len(self.__threadList):
            thread = self.__threadList.pop(0)
            if thread._error:
                if os.environ.get('ALLOW_ERROR') == '1':
                    self.__canProcBeAdjusted = True
                    continue
                msg = 'ERROR in thread:' + thread._error
                try:
                    import Log
                    Log.error(msg)
                except:
                    print msg
                raise NotImplementedError
            if thread.isAlive():
                newThreadList.append(thread)
            else:
                self.__canProcBeAdjusted = True
                if self.__keepResults:
                    # print 'KEEEP', thread._res
                    self._resList.append((thread._args, thread._res))
                if hasattr(thread, '_expectedOutputFile') and \
                   thread._expectedOutputFile:
                    self.__jobDict[thread._expectedOutputFile] = None
        self.__threadList = newThreadList

    def wait(self, callbackFunc=None, *args):
        while len(self.__threadList):
            thread = self.__threadList.pop(0)
            thread.join(self.__timeout)
            if not thread._error and self.__keepResults:
                if thread.isAlive():
                    # thread.stop()
                    # thread.join()
                    self._resList.append((thread._args, self.RUNNING))
                    continue
                else:
                    self._resList.append((thread._args, thread._res))
            if callbackFunc is not None:
                print 'call', thread._res
                callbackFunc(thread._res, *args)

    waitUntilClusterIsFree = wait

    def __getNbUsedThreads(self):
        return sum([thread._nbProc for thread in self.__threadList])

    def __canNewThreadBeCreated(self):
        self.__removeFinishedJobs()
        if self.__canProcBeAdjusted and self.__adjustProcNb:
            oldMax = self.__maxThreadNb
            self.__maxThreadNb = max(
                self.__maxThreadNb, _getNbAvailableCpus() +
                len(self.__threadList))
            if oldMax != self.__maxThreadNb:
                print 'Extending max thread nb from %d to %d' % \
                      (oldMax, self.__maxThreadNb)
        return self.__maxThreadNb is None or \
            self.__getNbUsedThreads() <= self.__maxThreadNb - self.__nbProc

    def __waitForThreadToFinish(self):
        while True:
            if self.__canNewThreadBeCreated():
                return
            time.sleep(1)

    def submit(self, func, *args):
        # print 'Submit', func, args
        # print self.__maxThreadNb, _getNbAvailableCpus(),
        # len(self.__threadList), self.__nbProc
        self.__waitForThreadToFinish()
        thread = ThreadFunc(func, *args)
        thread.daemon = True
        thread._nbProc = self.__nbProc
        if hasattr(self, '__expectedOutputFile'):
            thread._expectedOutputFile = self.__expectedOutputFile
        self.__removeFinishedJobs()
        self.__threadList.append(thread)
        # print 'Running thread'
        thread.start()

    def submitJobAndGetId(self, command, memory=None, nbProc=1,
                          dependentJobId=None, jobName=None,
                          machineToUseList=None, errorFile=None,
                          outputFile=None, queue=None,
                          machineToExcludeList=None, email=None,
                          walltime=604800, node=1, scriptName=None,
                          expectedOutputFile=None, optionList=None):
        appendTouch = False
        if ' touch ' in command:
            cmdList = command.split(' && ')
            if 'touch' in cmdList[-1].split():
                command = ' && '.join(cmdList[:-1])
                cmdList = [command, cmdList[-1]]
                appendTouch = True
        if outputFile and '>' not in command:
            command += ' 1> %s' % outputFile
        if errorFile:
            command += ' 2> %s' % errorFile
        if appendTouch:
            cmdList[0] = command
            command = ' && '.join(cmdList)
        self.__nbProc = nbProc
        self.__expectedOutputFile = expectedOutputFile
        self.submit(Utilities.mySystem, command, scriptName)


class RemoteThreadManager:

    def __init__(self, machineList):
        from sshUtils import Ssh
        _machineList = []
        self.__machineDict = defaultdict(int)
        for machineName in machineList:
            ssh = Ssh(machineName)
            ssh._machineName = machineName
            _machineList.append(ssh)
            self.__machineDict[machineName] += 1
        self._machineList = _machineList
        self._cluster = ThreadManager(len(machineList))
        self._jobId = 0
        self._jobIdDict = {}
        self._jobToRunList = []
        self.__machineToRemoveList = []

    def __updateMachineListAndDict(self):
        machineDict = defaultdict(int)
        if 'USE_REMOTE' not in os.environ:
            return
        for machineName in os.environ.get('USE_REMOTE').split(','):
            machineDict[machineName] += 1
        itemList = copy.deepcopy(self.__machineDict.items())
        for machineName, machineNb in itemList:
            currentMachineNb = machineDict[machineName]
            if machineNb > currentMachineNb:
                for i in range(currentMachineNb - machineName):
                    self.__machineToRemoveList.append(machineName)
            elif machineName < currentMachineNb:
                ssh = Ssh(machineName)
                ssh._machineName = machineName
                self._machineList.append(ssh)
            self.__machineDict[machineName] = currentMachineNb

    def getNbRunningJobs(self):
        return self._cluster.getNbRunningJobs()

    def __run(self, ssh, command):
        ssh.execCommand(command)
        self._machineList.append(ssh)
        self._jobIdDict[ssh._jobId] = 1

    def __appendOutAndErrorFilesToCommand(self, command, errorFile,
                                          outputFile):
        partList = command.split('&&')
        for i, part in enumerate(partList):
            currentPartList = partList[i].strip().split()
            if currentPartList[0] == 'touch':
                if outputFile:
                    partList[i - 1] += ' 1> %s' % outputFile
                if errorFile:
                    partList[i - 1] += ' 2> %s' % errorFile
        return '&&'.join(partList)

    def __getMachineForProcessing(self):
        while True:
            while not self._machineList:
                time.sleep(1)
            ssh = self._machineList.pop(0)
            if ssh._machineName in self.__machineToRemoveList:
                print 'USE_REMOTE has changed: removing machine %s from \
cluster' % ssh._machineName
                self.__machineToRemoveList.remove(ssh._machineName)
                continue
            else:
                break
        return ssh

    def __submitJobAndGetId(self, jobId, command, memory=None, nbProc=1,
                            dependentJobId=None, jobName=None,
                            machineToUseList=None, errorFile=None,
                            outputFile=None, queue=None,
                            machineToExcludeList=None, email=None,
                            walltime=604800, node=1):
        self.__updateMachineListAndDict()
        if dependentJobId:
            for currentJobId in dependentJobId:
                if self._jobIdDict.get(currentJobId) == 0:
                    self._jobToRunList.append((jobId, command, memory, nbProc,
                                               dependentJobId, jobName,
                                               machineToUseList, errorFile,
                                               outputFile, queue,
                                               machineToExcludeList, email))
                    # print '*' * 50
                    # print 'WAIT', command
                    # print self._jobId, self._jobToRunList, self._jobIdDict
                    # print '=' * 30
                    return
        command = self.__appendOutAndErrorFilesToCommand(
            command, errorFile, outputFile)
        ssh = self.__getMachineForProcessing()
        if 'python' in command:
            bashrcDict = {HpcScriptBase.CEPH: '/mnt/isilon1/code/myBashrc2'}
            command = 'unset PYTHONPATH && source %s && %s' % (
                bashrcDict[_guessLocation()], command)
        print 'Submitting cmd "%s" on machine %s' % \
              (command, ssh._remoteConnection.hostName)
        ssh._jobId = jobId
        self._cluster.submit(self.__run, ssh, command)

    def submitJobAndGetId(self, command, memory=None, nbProc=1,
                          dependentJobId=None, jobName=None,
                          machineToUseList=None, errorFile=None,
                          outputFile=None, queue=None,
                          machineToExcludeList=None, email=None,
                          walltime=604800, node=1):
        self._jobId += 1
        self._jobIdDict[self._jobId] = 0
        self.__submitJobAndGetId(self._jobId, command, memory, nbProc,
                                 dependentJobId, jobName, machineToUseList,
                                 errorFile, outputFile, queue,
                                 machineToExcludeList, email)
        return self._jobId

    def wait(self):
        print 'WAITING'
        while self._jobToRunList:
            args = self._jobToRunList.pop(0)
            self.__submitJobAndGetId(*args)
        self._cluster.wait()


class JobServer:

    def submit(self, func, args=(), depfuncs=(), modules=(), callback=None,
               callbackargs=(), group='default', globals=None):
        res = func(*args)
        if callback:
            callback(res, *callbackargs)
        return self.get_ncpus

    def waitForParallelJobs(self, callbackFunc=None, *args):
        return

    def get_ncpus(self):
        return 1


class ClusterBase:
    _nodeTuple = ()

    def __init__(self, nbCpus=None, nodeTuple=None, port=None, sleepTime=1):
        self._nbCpus = nbCpus
        if nodeTuple is not None:
            self._nodeTuple = nodeTuple
        if port:
            self._nodeTuple = tuple(
                [node + ':%d' % port for node in self._nodeTuple])
        # case where we do not want to parallelize i.e. linear case
        if nbCpus == 1:
            self._jobServer = JobServer()
        elif nbCpus != -1:
            if nbCpus is not None:
                # Creates jobserver with ncpus workers
                self._jobServer = pp.Server(
                    nbCpus, ppservers=self._nodeTuple, secret='')
            else:
                self._jobServer = pp.Server(
                    ppservers=self._nodeTuple, secret='')
            # print "Starting pp with", self._jobServer.get_ncpus(), "workers",
            # self._nodeTuple
        self.__sleepTime = sleepTime
        self.__jobList = []

    def getNbSubmittedJobs(self):
        return len(self.__jobList)

    def getNbActiveJobs(self):
        return len([job for job in self.__jobList if not
                    isinstance(job, types.StringType) and not job.finished])

    getNbRunningJobs = getNbActiveJobs

    def __waitForJobs(self):
        while self.getNbActiveJobs() > 10 * self._jobServer.get_ncpus():
            # print 'There'
            time.sleep(self.__sleepTime)
            # newJobList = []
            # self.__jobList = [job for job in self.__jobList if type(job) ==
            # types.StringType or not job.finished]

    def submit(self, func, paramTuple, dependentFuncTuple=(), importTuple=(),
               callBack=None, callBackArgs=(), cacheFileName=None,
               maxNbProcess=None):
        if self._nbCpus == -1:
            if cacheFileName:
                if Utilities.doesFileExist(cacheFileName):
                    res = Utilities.loadCache(cacheFileName)
                else:
                    res = Utilities.getFunctionResultWithCache(
                        cacheFileName, func, *paramTuple)
            else:
                res = func(*paramTuple)
            self.__jobList.append(res)
            if len(res) == 0:
                raise NotImplementedError
            return
        job = cacheFileName
        if cacheFileName is None or not Utilities.doesFileExist(cacheFileName):
            if maxNbProcess is None:
                job = self._jobServer.submit(func, paramTuple,
                                             dependentFuncTuple, importTuple,
                                             callback=callBack,
                                             callbackargs=callBackArgs)
            else:
                job = self._jobServer.submit(
                    func, paramTuple, dependentFuncTuple, importTuple)
            if cacheFileName:
                job._cacheFileName = cacheFileName
        self.__jobList.append(job)
        self.__waitForJobs()
        if maxNbProcess is not None and len(self.__jobList) > maxNbProcess:
            self.wait(callBack, *callBackArgs)
        return job

    def submitList(self, func, paramList, dependentFuncTuple, importTuple,
                   step=None, callBack=None, callBackArgs=()):
        if step is None:
            step = len(paramList)
        resList = []
        for currentParamList in iterRange(paramList, step):
            for paramTuple in currentParamList:
                if not isinstance(paramTuple, types.TupleType):
                    paramTuple = (paramTuple, )
                resList.append(self.submit(
                    func, paramTuple, dependentFuncTuple, importTuple,
                    callBack, callBackArgs))
        return resList

    def submitJobAndGetId(self, command, memory=None, nbProc=1,
                          dependentJobId=None, jobName=None,
                          machineToUseList=None, errorFile=None,
                          outputFile=None, queue=None, arrayIdxStr=None,
                          machineToExcludeList=None):
        if errorFile:
            command += ' 2> %s' % errorFile
        self.submit(T().run, (command, ))

    def removeJob(self, job):
        self.__jobList.remove(job)

    def __getResFromJob(self, job):
        if isinstance(job, types.StringType):
            res = Utilities.loadCache(job)
        else:
            if hasattr(job, '_cacheFileName') and job._cacheFileName:
                res = Utilities.getFunctionResultWithCache(
                    job._cacheFileName, job)
            else:
                res = job()
        return res

    def waitForParallelJobs(self, callbackFunc=None, *args):
        if callbackFunc is None:
            # self._jobServer.wait()
            if self._nbCpus == -1:
                return
            while len(self.__jobList):
                job = self.__jobList.pop(0)
                self.__getResFromJob(job)
            return
        # print 'WAITING'
        while len(self.__jobList):
            job = self.__jobList.pop(0)
            res = job
            if self._nbCpus != -1:
                res = self.__getResFromJob(job)
            # print res
            if callbackFunc:
                # print 'CALLING'
                if res is None:
                    raise NotImplementedError
                callbackFunc(res, *args)
            # print 'JOB DONE'

    wait = waitForParallelJobs
    waitUntilClusterIsFree = waitForParallelJobs


class BeijingCluster(ClusterBase):
    _nodeTuple = tuple(['node%d' % i for i in xrange(1, 20) if i != 8])


class CngCluster(ClusterBase):
    _nodeTuple = ('192.168.38.201', '192.168.38.200')


class XserverCluster(ClusterBase):
    _nodeTuple = ('xserver1.seq.cng.fr', 'xserver2.seq.cng.fr',
                  'xserver3.seq.cng.fr', 'xserver4.seq.cng.fr',
                  'xserver5.seq.cng.fr')


def getClusterClass():
    if os.environ.get('HOSTNAME') == 'node8':
        cluster = BeijingCluster
    else:
        cluster = ClusterBase
    return cluster


class HpcBase:
    JOB_ID_VAR = None
    TASK_ID_VAR = None
    MASTER_VAR = None

    _startJobIdStr = None
    _endJobIdStr = None
    _defaultMemory = 4
    _nbTries = 1
    _deleteJobCmd = None

    def __init__(self):
        self._jobList = []
        if self.__createJobList():
            self.__targetDir, self.__nbParts = os.environ.get(
                'CREATE_JOBS').split(',')
            self.__nbParts = int(self.__nbParts)

    def _extractItemFromStr(self, jobStr, item):
        idx = jobStr.find(item)
        if idx == -1:
            raise NotImplementedError('Can not find "%s" file in [%s]' %
                                      (item, jobStr))
        return jobStr[idx:].split('=')[1].split()[0]

    def _addOptionStrToCmd(self, cmd, optionList):
        for optionName, optionValue in optionList:
            cmd += ' %s %s' % (optionName, optionValue)
        return cmd

    def deleteJob(self, jobId):
        if isinstance(jobId, types.IntType):
            jobId = [jobId]
        if not self._deleteJobCmd:
            raise NotImplementedError('Delete cmd not set')
        Utilities.mySystem('%s %s' % (self._deleteJobCmd,
                                      ' '.join([str(jId) for jId in jobId])))

    def __createJobList(self):
        return os.environ.get('CREATE_JOBS')

    def __del__(self):
        if self.__createJobList() and self._jobList:
            import math
            Utilities.mySystem('mkdir -p %s' % self.__targetDir)
            nbJobs = int(math.ceil(1. * len(self._jobList) / self.__nbParts))
            print 'Creating %d job dump files for %d jobs with %d jobs per \
file in %s' % (self.__nbParts, len(self._jobList), nbJobs, self.__targetDir)
            for i, jobList in enumerate(iterRange(self._jobList, nbJobs)):
                dumpFile = os.path.join(
                    self.__targetDir, '%d-%d.pyDump' % (i, self.__nbParts))
                Utilities.saveCache(jobList, dumpFile)

    def getNbRunningJobs(self):
        raise NotImplementedError

    def getJobIdList(self):
        raise NotImplementedError

    def getJobDetails(self, jobId):
        raise NotImplementedError

    def getOutErrFileForJobId(self, jobId):
        jobStr = self.getJobDetails(jobId)
        return self._extractOutErrFileFromJobStr(jobStr, jobId)

    def _extractOutErrFileFromJobStr(self, jobStr, jobId):
        raise NotImplementedError

    def isMasterNode(self):
        # for key in os.environ:
        # if self.MASTER_VAR in key:
        # return True
        return os.environ.get('HOSTNAME', '') in ['Unicluster', 'master1']

    def areAllVarInEnv(self):
        for varName in [self.JOB_ID_VAR, self.TASK_ID_VAR]:
            if varName not in os.environ:
                return
        return True

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue=None, arrayIdxStr=None, email=None,
                   walltime=604800, node=1, optionList=None):
        raise NotImplementedError

    def _extractJobIdFromStr(self, jobIdStr):
        return jobIdStr

    def _submitRawStrAndgetJobId(self, jobStr):
        from WebExtractor import WebExtractor
        fh = os.popen(jobStr)
        content = fh.read()
        for i in range(self._nbTries):
            try:
                jobIdStr = WebExtractor()._getStrIncludedInTag(
                    content, self._startJobIdStr, self._endJobIdStr)
                jobIdStr = self._extractJobIdFromStr(jobIdStr)
                jobId = int(jobIdStr)
                break
            except:
                print 'Unhandled', [content]
                if i == self._nbTries - 1:
                    raise
                print 'Submission failed i = %d, retrying...' % i
        fh.close()
        return jobId

    def _getMachineList(self):
        raise NotImplementedError

    def submitJobAndGetId(self, command, memory=None, nbProc=1,
                          dependentJobId=None, jobName=None,
                          machineToUseList=None, errorFile=None,
                          outputFile=None, queue=None, arrayIdxStr=None,
                          machineToExcludeList=None, email=None,
                          walltime=604800, node=1, nextJobDict=None,
                          scriptName=None, expectedOutputFile=None,
                          optionList=None):
        if queue and ',' in queue:
            queue = tuple(queue.split(','))
        if not email:
            email = os.environ.get('CLUSTER_MAIL')
        if not machineToExcludeList:
            machineToExcludeList = os.environ.get(
                'MACHINE_TO_EXCLUDE', '').split(',')
            if machineToExcludeList == ['']:
                machineToExcludeList = None
        if memory is None:
            memory = self._defaultMemory
        if machineToExcludeList:
            if machineToUseList:
                raise NotImplementedError
            print 'machineToExcludeList = ', machineToExcludeList
            machineToUseList = list(
                set(self._getMachineList()) - set(machineToExcludeList))
            machineToUseList.sort()
        print 'cmd = [%s]' % command
        if nextJobDict:
            dependentJobId = None
        if scriptName:
            from FileHandler import CsvFileWriter
            outFh = CsvFileWriter(scriptName)
            shaBang = '#!/bin/bash'
            if shaBang not in command:
                command = shaBang + '\n' + command
            outFh.write(command)
            outFh.close()
            Utilities.mySystem('chmod +x %s' % scriptName)
            command = scriptName
        if self.__createJobList():
            if not expectedOutputFile:
                raise NotImplementedError(
                    'expectedOutputFile is None for command: %s' % command)
            self._jobList.append((command, expectedOutputFile))
            return
        jobStr = self._getJobStr(command, memory, nbProc, dependentJobId,
                                 jobName, machineToUseList, errorFile,
                                 outputFile, queue, arrayIdxStr, email,
                                 walltime, node, optionList)
        print 'cmd = [%s]' % jobStr
        if isinstance(jobStr, types.StringType):
            jobList = [jobStr]
        else:
            jobList = jobStr
        currentJobIdList = []
        if nextJobDict and len(jobList) > 1:
            raise NotImplementedError(
                "jobList's length should be <= 1 when nextJobDict is not null:\
 %s %s" % (jobList, nextJobDict))
        for jobStr in jobList:
            # print 'job = [%s]' % jobStr
            jobId = self._submitRawStrAndgetJobId(jobStr)
            print 'Submitted with jobId %d' % jobId
            self._jobList.append(jobId)
            currentJobIdList.append(jobId)
        if len(currentJobIdList) == 1:
            currentJobIdList = currentJobIdList[0]
        return currentJobIdList

    def submitJobs(self, commandAndParamList):
        jobId = None
        for commandAndParams in commandAndParamList:
            memory = 4
            nbProc = 1
            jobName = None
            if len(commandAndParams) == 2:
                command, memory = commandAndParams
            elif len(commandAndParams) == 3:
                if isinstance(commandAndParams[-1], types.IntType):
                    command, memory, nbProc = commandAndParams
                else:
                    command, memory, jobName = commandAndParams
            elif len(commandAndParams) == 4:
                command, memory, nbProc, jobName = commandAndParams
            elif len(commandAndParams) == 5:
                command, memory, nbProc, jobName, jobId = commandAndParams
            paramList = [command, memory, nbProc, jobId, jobName]
            jobId = self.submitJobAndGetId(*paramList)

    def getTaskIdx(self):
        try:
            return int(os.environ.get(self.TASK_ID_VAR, 0)) - 1
        except:
            return

    def getJobId(self):
        return int(os.environ.get(self.JOB_ID_VAR, 0))

    def _isJobFinished(self, jobId):
        raise NotImplementedError

    def waitUntilClusterIsFree(self):
        while True:
            nbJobs = self.getNbRunningJobs()
            if not nbJobs:
                break
            print 'Waiting %ds' % nbJobs
            time.sleep(nbJobs)

    def wait(self):
        while self._jobList:
            jobId = self._jobList.pop()
            if not self._isJobFinished(jobId):
                self._jobList.append(jobId)

    def waitForJob(self, jobId, sleepTime=1):
        while True:
            time.sleep(sleepTime)
            if self._isJobFinished(jobId):
                break


class SGEcluster(HpcBase):
    JOB_ID_VAR = 'JOB_ID'
    TASK_ID_VAR = 'SGE_TASK_ID'
    MASTER_VAR = 'SGE'

    _startJobIdStr = 'Your job'
    _endJobIdStr = ' ('
    _deleteJobCmd = 'qdel'

    def __extractValueForKeywordInJobStr(self, keyword, jobStr):
        if keyword not in jobStr:
            raise NotImplementedError('keyword "%s" not found in:\n"%s"' %
                                      (keyword, jobStr))
        return jobStr.split(keyword)[-1].split('\n')[0].split(':')[-1].strip()

    def _extractOutErrFileFromJobStr(self, jobStr, jobId):
        stdErrFileName = self.__extractValueForKeywordInJobStr(
            'stderr_path_list:', jobStr)
        stdOutFileName = self.__extractValueForKeywordInJobStr(
            'stdout_path_list:', jobStr)
        return stdOutFileName, stdOutFileName

    def getJobIdList(self):
        fh = os.popen('qstat')
        jobIdList = []
        if not fh.readline().strip():
            return jobIdList
        line = fh.readline()
        if '---------------------------------' not in line:
            raise NotImplementedError(
                '"---------------------------------" expected in line "%s"' %
                line)
        for line in fh:
            splittedLine = line.split()
            status = splittedLine[4]
            machine = splittedLine[7]
            jobIdList.append((int(splittedLine[0]), status, machine))
        return jobIdList

    def getJobDetails(self, jobId):
        fh = os.popen('qstat -j %d' % jobId)
        return fh.read()

    def _getMachineList(self):
        return ['c1l%d' % i for i in range(1, 17)]

    def _isJobFinished(self, jobId):
        fh = os.popen('qstat')
        content = fh.read()
        fh.close()
        print 'Content = [%s]' % content
        return ' %d ' % jobId not in content

    def _extractJobIdFromStr(self, jobIdStr):
        if not ValueParser.isNb(jobIdStr):
            jobIdStr = jobIdStr.split()[-1].split('.')[0]
        return int(jobIdStr)

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue=None, arrayIdxStr=None, email=None, walltime=None,
                   node=None, optionList=None):
        if nbProc > 1:
            memory = max(memory, nbProc * 3)
        cmd = 'echo "%s" | qsub -l h_vmem=%dG,virtual_free=%dG ' % (
            command, memory, memory)
        if dependentJobId:
            if not isinstance(dependentJobId, types.ListType):
                dependentJobId = [dependentJobId]
            cmd += '-hold_jid %s ' % (','.join([str(jobId)
                                                for jobId in dependentJobId]))
        if jobName:
            switch = 'N'
            if ValueParser.isNb(jobName.split('-')[0]):
                switch = 't'
            cmd += '-%s %s ' % (switch, jobName)
        if errorFile:
            cmd += '-e %s ' % errorFile
        if outputFile:
            cmd += '-o %s ' % outputFile
        if jobName and '[' in jobName:
            arrayIdxStr = jobName.split('[')[-1].split(']')[0]
            cmd += '-t %s ' % arrayIdxStr
        if machineToUseList:
            cmd += '-q %s ' % (','.join(['*@%s' %
                                         machineName for machineName in
                                         machineToUseList]))
        if email:
            cmd += '-m eas -M %s ' % email
        # if nbProc > 1:
            # cmd += '-pe make %d ' % nbProc
        print 'JOB -> %s' % cmd
        return cmd


class LSFcluster(HpcBase):
    JOB_ID_VAR = 'LSB_JOBID'
    TASK_ID_VAR = 'LSB_JOBINDEX'
    MASTER_VAR = 'LSF'

    _startJobIdStr = '<'
    _endJobIdStr = '>'

    def _isJobFinished(self, jobId):
        fh = os.popen('bjobs')
        content = fh.read()
        fh.close()
        # print 'Content = [%s]' % content
        return '%d ' % jobId not in content

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue=None, arrayIdxStr=None, email=None, walltime=None,
                   node=None, optionList=None):
        cmd = 'bsub -n %d -R "rusage[mem=%d]" -M %d ' % (
            nbProc, memory, memory)
        if jobName:
            cmd += '-J "%s" ' % jobName
        if dependentJobId:
            if not isinstance(dependentJobId, types.ListType):
                dependentJobId = [dependentJobId]
            cmd += '-w "%s" ' % (' && '.join([str(jobId)
                                              for jobId in dependentJobId]))
        if machineToUseList:
            cmd += '-m "%s" ' % (' '.join(machineToUseList))
        if errorFile:
            cmd += '-e %s ' % errorFile
        if outputFile:
            cmd += '-o %s ' % outputFile
        return cmd + '"%s"' % command


class MOABcluster(HpcBase):
    JOB_ID_VAR = 'LSB_JOBID'
    TASK_ID_VAR = 'LSB_JOBINDEX'
    MASTER_VAR = 'PBS'

    _startJobIdStr = '\n'
    _endJobIdStr = '\n'
    _defaultMemory = 3
    _nbTries = 10

    def getNbRunningJobs(self):
        fh = os.popen('showq -u $USER')
        nbJobs = fh.readlines()[-2].strip().split()[-1]
        return int(nbJobs)

    def _isJobFinished(self, jobId):
        fh = os.popen('checkjob -v %d' % jobId)
        content = fh.read()
        fh.close()
        return 'State: Completed' in content

    def __getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                    machineToUseList, errorFile=None, outputFile=None,
                    queue=None, arrayIdxStr=None, email=None, walltime=604800,
                    node=1, optionList=None):
        partList = command.split()
        if len(partList) > 1 and partList[0].split('.')[-1] == 'sh':
            fh = open(partList[0])
            command = fh.read()
            fh.close()
            if '$*' not in command:
                raise NotImplementedError(
                    'Expected to find "$*" in %s but found "%s"' %
                    (partList[0], command))
            command = command[:command.rfind('$*')]
            command += ' ' + ' '.join(partList[1:])
        # print 'CMD [%s]' % command
        dependentStr = ''
        if dependentJobId:
            if not isinstance(dependentJobId, types.ListType):
                dependentJobId = [dependentJobId]
            dependentStr = ',depend=%s' % (
                ':'.join([str(jobId) for jobId in dependentJobId]))
        memoryStr = ''
        if memory != self._defaultMemory:
            memoryStr = ',pmem=%dgb' % memory
        cmd = 'echo "%s" | msub -l nodes=%d:ppn=%d%s,walltime=%d%s ' % (
            command, node, nbProc, memoryStr, walltime, dependentStr)
        if not queue:
            queue = 'sw'
        if queue:
            cmd += '-q %s ' % queue
        print '==cmd = [%s]==' % cmd
        if jobName:
            cmd += '-N %s ' % jobName
        if errorFile:
            cmd += '-e %s ' % errorFile
        if outputFile:
            cmd += '-o %s ' % outputFile
        return cmd

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue='sw', arrayIdxStr=None, email=None, walltime=604800,
                   node=1, optionList=None):
        idxList = None
        if jobName and '[' in jobName and ']' in jobName:
            idxStr = jobName.split('[')[-1].split(']')[0]
            idxList = ValueParser().getIdxListFromIntervalStrOrList(idxStr)
        cmd = self.__getJobStr(command, memory, nbProc, dependentJobId,
                               jobName, machineToUseList, errorFile,
                               outputFile, queue, arrayIdxStr, email,
                               walltime, node, optionList)
        if idxList:
            cmdList = []
            keyword = 'echo "'
            partList = cmd.split(keyword)
            cmd = partList[0] + keyword + '%%s;%s' % partList[1]
            if '[' not in jobName or ']' not in jobName:
                raise NotImplementedError(
                    'Expected to find "[" and "]" in jobName "%s"' % jobName)
            toReplace = '[' + jobName.split('[')[1].split(']')[0] + ']'
            # print 'To Repl ]%s[' % toReplace
            for idx in idxList:
                cmd2 = cmd % 'export %s=%d' % (self.TASK_ID_VAR, idx)
                cmd2 = cmd2.replace('-N %s ' % jobName, '-N %s ' %
                                    (jobName.replace(toReplace, '[%d]' % idx)))
                cmdList.append(cmd2)
            # print '>>>>>>>>>>'
            # print cmdList
            # print '<<<<<<<<<<'
            return cmdList
        return cmd


class SLURMcluster(HpcBase):
    JOB_ID_VAR = 'SLURM_ARRAY_JOB_ID'
    TASK_ID_VAR = 'SLURM_ARRAY_TASK_ID'

    _startJobIdStr = 'Submitted batch job '
    _endJobIdStr = '\n'
    _defaultMemory = 4
    _deleteJobCmd = 'scancel'

    def _extractOutErrFileFromJobStr(self, jobStr, jobId):
        if self._extractItemFromStr(jobStr, 'JobState=').split()[0] == \
           'PENDING':
            return
        return self._extractItemFromStr(jobStr, 'StdOut='), \
            self._extractItemFromStr(jobStr, 'StdErr=')

    def getJobDetails(self, jobId):
        fh = os.popen('sjinfo %d' % jobId)
        return fh.read()

    def getJobIdList(self):
        for jobStr in os.popen('squeue -u %s' % os.environ['USER']):
            partList = jobStr.split()
            if partList[0] == 'JOBID':
                continue
            yield int(partList[0]), partList[4], partList[-1]

    def _getMachineList(self):
        return ['node%02d' % i for i in range(19, 23) + range(1, 7)]

    def __getScriptNameForJob(self, command, scriptName=None):
        if not scriptName:
            import random
            random.seed(time.time())
            scriptName = str(random.random()) + '.sh'
            scriptName = '/env/cng/proj/projet_LIVER_356/scratch/src/victor/\
jobs/%s' % scriptName
        fh = open(scriptName, 'w')
        fh.write('#!/bin/bash\n%s' % command)
        fh.close()
        Utilities.mySystem('chmod +x %s' % scriptName)
        return scriptName

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue='sw', arrayIdxStr=None, email=None, walltime=86400,
                   node=1, optionList=None):
        adjustMemory = True
        if _guessLocation() == HpcScriptBase.CEPH_SLURM:
            # print 'SOKO -> [%s]' % command
            fullCmd = command
            if os.path.isfile(fullCmd):
                fullCmd = open(command).read()
            if ' -p shell -f ' in fullCmd:
                fullCmd = open(fullCmd.split(' -p shell -f ')
                               [1].split()[0]).read()
                # print 'SOKO2 -> [%s]' % fullCmd
            if 'GenomeAnalysisTK' in fullCmd or 'java' in fullCmd:
                adjustMemory = False
                currentMemory = os.environ.get('MEMORY')
                if currentMemory:
                    memory = int(currentMemory)
                factor = 4
                currentFactor = os.environ.get('FACTOR')
                if currentFactor:
                    factor = int(currentFactor)
                nbProc *= factor
                # nbProc = min(_detect_ncpus(), nbProc)
                if memory > 10:
                    memory /= 2
            if walltime:
                walltime *= 4
        if adjustMemory and nbProc:
            memory /= nbProc
        optionStr = ''
        if not errorFile and os.path.isfile(command):
            arrayStr = jobName.split('[')[1].split(']')[0]
            start, end = arrayStr.split('-')
            optionStr = '--array=%d-%d' % (int(start), int(end))
            scriptName = command
        else:
            scriptName = self.__getScriptNameForJob(
                command, FileNameGetter(errorFile).get('_cmd.sh'))
        # scriptName = FileNameGetter(errorFile).get('_cmd.sh')
        cmd = 'sbatch -c %d --mem-per-cpu=%d %s ' % (
            nbProc, memory * 1000, optionStr)
        qos = None
        if isinstance(queue, types.TupleType):
            queue, qos = queue
        if jobName:
            cmd += '-J %s ' % jobName
        if errorFile:
            cmd += '-e %s ' % errorFile
        if outputFile:
            cmd += '-o %s ' % outputFile
        if queue:
            qName = queue
            if isinstance(queue, types.ListType):
                qName, nb = queue[0]
                nb -= 1
                if not nb:
                    queue.pop(0)
                else:
                    queue[0] = qName, nb
            cmd += '-p %s ' % qName
        if qos:
            cmd += '--qos %s ' % qos
        if email:
            cmd += '--mail-type=END '
        if walltime:
            cmd += '-t %d ' % (walltime / 60)
        if dependentJobId:
            if not isinstance(dependentJobId, types.ListType):
                dependentJobId = [dependentJobId]
            else:
                dependentJobId = [jobId for jobId in dependentJobId if jobId]
            if dependentJobId:
                cmd += '-d afterok:%s ' % (':'.join([str(jobId)
                                                     for jobId in
                                                     dependentJobId]))
        if machineToUseList:
            cmd += '-w %s ' % (','.join(machineToUseList))
        if optionList:
            cmd = self._addOptionStrToCmd(cmd, optionList)
        cmd += ' %s' % scriptName
        print [cmd]
        return cmd


class SLURMCCRTcluster(HpcBase):
    JOB_ID_VAR = 'SLURM_ARRAY_JOB_ID'
    TASK_ID_VAR = 'SLURM_ARRAY_TASK_ID'

    _startJobIdStr = 'Submitted Batch Session '
    _endJobIdStr = '\n'
    _defaultMemory = 4
    _deleteJobCmd = 'scancel'

    def _extractOutErrFileFromJobStr(self, jobStr, jobId):
        fh = os.popen('ccc_mstat -r %d' % jobId)
        content = fh.read()
        # print [content]
        # print self._extractItemFromStr(content, 'StdOut='),
        # self._extractItemFromStr(content, 'StdErr=')
        return self._extractItemFromStr(content, 'StdOut='), \
            self._extractItemFromStr(content, 'StdErr=')

    def getJobDetails(self, jobId):
        fh = os.popen('ccc_mstat -b %d' % jobId)
        return fh.read()

    def getJobIdList(self):
        fh = os.popen('ccc_mstat -u $USER')
        content = fh.read()
        lastIdx = content.rfind('------')
        lineList = content[lastIdx:].strip().split('\n')
        print lineList.pop(0)
        return [int(line.split()[0]) for line in lineList]

    def submitJobAndGetId(self, command, memory=None, nbProc=1,
                          dependentJobId=None, jobName=None,
                          machineToUseList=None, errorFile=None,
                          outputFile=None, queue=None, arrayIdxStr=None,
                          machineToExcludeList=None, email=None,
                          walltime=604800, node=1, nextJobDict=None,
                          scriptName=None, expectedOutputFile=None,
                          optionList=None):
        jobId = HpcBase.submitJobAndGetId(self, command, memory, nbProc,
                                          dependentJobId, jobName,
                                          machineToUseList, errorFile,
                                          outputFile, queue, arrayIdxStr,
                                          machineToExcludeList, email,
                                          walltime, node, nextJobDict,
                                          optionList=optionList,
                                          scriptName=scriptName)
        if walltime == -1:
            os.system('scontrol hold %d' % jobId)
            os.system('scontrol release %d' % jobId)
        return jobId

    def _getJobStr(self, command, memory, nbProc, dependentJobId, jobName,
                   machineToUseList, errorFile=None, outputFile=None,
                   queue='broadwell', arrayIdxStr=None, email=None,
                   walltime=86400, node=1, optionList=None):
        if not queue:
            queue = 'broadwell'
        qos = None
        maxWallTime = 86400
        projectName = 'fg0016'
        if os.environ.get('PROJECT_NAME'):
            projectName = os.environ.get('PROJECT_NAME')
        fileSystem = None
        if os.environ.get('FILE_SYSTEM'):
            fileSystem = os.environ.get('FILE_SYSTEM')
        if walltime > maxWallTime and not isinstance(queue, types.TupleType):
            queue = queue, 'long'
        if isinstance(queue, types.TupleType):
            queue, qos = queue
        if qos == 'long':
            if walltime <= maxWallTime:
                qos = None
            maxWallTime *= 3
        elif qos == 'test':
            walltime = -1
        if not walltime:
            walltime = maxWallTime
        if walltime > maxWallTime:
            print 'Warning: walltime is too high: [%s] %d' % (command,
                                                              walltime)
            walltime = maxWallTime
        if walltime == -1:
            walltime = None
        # memory = min(memory, 8)
        cmd = 'echo "%s" | ccc_msub -c %d -M %d -A %s -q %s ' % (
            command, nbProc, memory * 1000, projectName, queue)
        if qos:
            cmd += '-Q %s ' % qos
        if walltime:
            cmd += '-T %d ' % walltime
        if fileSystem:
            cmd += '-m %s ' % fileSystem
        if jobName:
            if '[' in jobName:
                arrayIdxStr = jobName.split('[')[-1].split(']')[0]
                step = os.environ.get('STEP')
                if step:
                    arrayIdxStr += ':%d' % int(step)
                jobName = jobName.split('[')[0]
                cmd += '-E "--array=%s" ' % arrayIdxStr
            cmd += '-r %s ' % jobName
        if errorFile:
            cmd += '-e %s ' % errorFile
        if outputFile:
            cmd += '-o %s ' % outputFile
        if email:
            cmd += '-@ %s:end ' % email
        if dependentJobId:
            if not isinstance(dependentJobId, types.ListType):
                dependentJobId = [dependentJobId]
            cmd += '-a %s' % (','.join([str(jobId)
                                        for jobId in dependentJobId]))
        print [cmd]
        return cmd


def __isMasterNode():
    return 'vmaster' in os.popen('hostname').read()
    # stdin, stdout, stderr = os.popen3('qconf -sconf %s' %
    # os.popen('hostname').read().strip())
    # return stdout.read() and 'not defined' not in stderr.read()


def guessHpcOld(allowLocalRun=False, nbCpus=None, *args):
    if not nbCpus:
        nbCpus = _getNbAvailableCpus()
    if _guessLocation() == HpcScriptBase.CNG_MAC:
        return ClusterBase(nbCpus)
    msg = 'Warning cluster not found, using ThreadManager instead with %d \
cpus' % nbCpus
    if os.environ.get('USE_REMOTE'):
        print 'REMOTE: %s' % msg
        machineList = os.environ['USE_REMOTE'].split(',')
        print 'Machines', machineList
        return RemoteThreadManager(machineList)
    if os.environ.get('USE_CLUSTER') == '0':
        print msg
        return ThreadManager(nbCpus, *args)
    if 'TORTUGA_ROOT' in os.environ:
        return LSFcluster()
    for clusterClass in [SGEcluster, LSFcluster]:
        cluster = clusterClass()
        if cluster.areAllVarInEnv():
            if isinstance(cluster, SGEcluster) and not __isMasterNode():
                cluster = ThreadManager(nbCpus, *args)
            return cluster
    if os.environ.get('CLUSTER_NAME') == 'lirac':
        return SLURMcluster()
    for envName in os.environ:
        if envName[:4] == 'SGE_':
            if not __isMasterNode():
                return ThreadManager(nbCpus, *args)
            return SGEcluster()
        elif envName[:4] == 'LSF_':
            return LSFcluster()
        elif envName[:4] == 'PBS_' or envName == 'MODULEPATH':
            return MOABcluster()
    if allowLocalRun:
        print msg
        return ThreadManager(nbCpus, *args)


def __isMasterNode():
    return 'vmaster' in os.popen('hostname').read()
    # stdin, stdout, stderr = os.popen3('qconf -sconf %s' %
    # os.popen('hostname').read().strip())
    # return stdout.read() and 'not defined' not in stderr.read()


def guessHpc(allowLocalRun=False, nbCpus=None, useFromMasterNode=True, *args):
    cmdAndClusterTypeList = [('qsub', SGEcluster),
                             ('ccc_msub', SLURMCCRTcluster),
                             ('sbatch', SLURMcluster), ('bsub', LSFcluster),
                             ('msub', MOABcluster)]
    if not nbCpus:
        nbCpus = _getNbAvailableCpus()
    msg = 'Warning cluster not found, using ThreadManager instead with %d \
cpus' % nbCpus
    if os.environ.get('USE_REMOTE'):
        print 'REMOTE: %s' % msg
        machineList = os.environ['USE_REMOTE'].split(',')
        print 'Machines', machineList
        return RemoteThreadManager(machineList)
    for cmd, clusterClass in cmdAndClusterTypeList:
        # fh = os.popen('%s -h 2> /dev/null' % cmd)
        fh = os.popen('which %s 2> /dev/null' % cmd)
        if fh.read():
            if clusterClass == SGEcluster and useFromMasterNode and not \
               __isMasterNode():
                return ThreadManager(nbCpus, *args)
            return clusterClass()
    if allowLocalRun:
        print msg
        return ThreadManager(nbCpus, *args)
    return ThreadManager(nbCpus, *args)


def guessHpcOld2(allowLocalRun=False, nbCpus=None, useFromMasterNode=True,
                 *args):
    if not nbCpus:
        nbCpus = _getNbAvailableCpus()
    if _guessLocation() == HpcScriptBase.CNG_MAC:
        return ClusterBase(nbCpus)
    if _guessLocation() == HpcScriptBase.CEPH_SLURM:
        return SLURMcluster()
    msg = 'Warning cluster not found, using ThreadManager instead with %d \
cpus' % nbCpus
    if os.environ.get('USE_REMOTE'):
        print 'REMOTE: %s' % msg
        machineList = os.environ['USE_REMOTE'].split(',')
        print 'Machines', machineList
        return RemoteThreadManager(machineList)
    if os.environ.get('USE_CLUSTER') == '0':
        print msg
        return ThreadManager(nbCpus, *args)
    if 'TORTUGA_ROOT' in os.environ:
        return LSFcluster()
    for clusterClass in [SGEcluster, LSFcluster]:
        cluster = clusterClass()
        if cluster.areAllVarInEnv():
            if isinstance(cluster, SGEcluster) and useFromMasterNode and \
               not __isMasterNode():
                cluster = ThreadManager(nbCpus, *args)
            return cluster
    if os.environ.get('CLUSTER_NAME') == 'lirac':
        return SLURMcluster()
    elif os.environ.get('CLUSTER_NAME') == 'cobalt':
        return SLURMCCRTcluster()
    for envName in os.environ:
        if envName[:4] == 'SGE_':
            if useFromMasterNode and not __isMasterNode():
                return ThreadManager(nbCpus, *args)
            return SGEcluster()
        elif envName[:4] == 'LSF_':
            return LSFcluster()
        elif envName[:4] == 'PBS_' or envName == 'MODULEPATH':
            return MOABcluster()
    if allowLocalRun:
        print msg
        return ThreadManager(nbCpus, *args)


def _getTaskIdx():
    for clusterClass in [SGEcluster, SLURMCCRTcluster, LSFcluster,
                         MOABcluster]:
        taskIdx = clusterClass().getTaskIdx()
        if taskIdx is not None and taskIdx != -1:
            return taskIdx


def _guessLocation():
    cluster = None
    location = HpcScriptBase.EBI
    hostName = os.popen('hostname').read().strip()
    if '.cng.fr' in os.environ.get('HOSTNAME', '') or '193.50' in \
       os.popen('ifconfig 2> /dev/null').read():
        if os.path.isfile('/mach_kernel'):
            location = HpcScriptBase.CNG_MAC
        else:
            location = HpcScriptBase.CNG
    elif 'cobalt' in os.environ.get('HOSTNAME', ''):
        location = HpcScriptBase.CCRT_SLURM
    elif 'mcgill' in os.environ.get('HOSTNAME', ''):
        location = HpcScriptBase.MC_GILL
    elif hostName == 'vmaster1':
        location = HpcScriptBase.CEPH_SLURM
    elif cluster is None or isinstance(cluster, SGEcluster):
        location = HpcScriptBase.CEPH
    elif 'MODULEPATH' in os.environ:
        location = HpcScriptBase.GUILLIMIN
    return location


class HpcScriptBase:
    _idx = None
    _jobId = None

    CNG = 0
    EBI = 1
    CEPH = 2
    CNG_MAC = 3
    MC_GILL = 4
    GUILLIMIN = 5
    CNG_SLURM = 6
    CCRT_SLURM = 7
    CEPH_SLURM = 8

    def __init__(self, nbCpus=30):
        self._nbCpus = nbCpus
        self._cluster = guessHpc()
        self._location = _guessLocation()
        if self._cluster and not isinstance(self._cluster, ThreadManager):
            self._idx = self._cluster.getTaskIdx()
            self._jobId = self._cluster.getJobId()
        if self._location in [self.CNG, self.CEPH]:
            self._clusterClass = SGEcluster
        else:
            self._clusterClass = LSFcluster

    def _getParamList(self, *args):
        raise NotImplementedError

    def _processAll(self, *args):
        threadManager = ThreadManager(self._nbCpus)
        for param in self._getParamList(*args[:-1]):
            args = list(args[:-1]) + [param]
            print args
            threadManager.submit(self._action, *args)
        threadManager.wait()

    def process(self, *args):
        if not self._cluster:
            self._processAll(*args)
        else:
            self._process(*args)

    def _process(self, *args):
        paramList = self._getParamList(*args[:-1])
        idx = args[-1]
        print 'IDX', idx, args
        if idx is not None:
            self._idx = idx
        param = paramList[self._idx]
        args = list(args[:-1]) + [param]
        self._action(*args)

    def _action(self, *args):
        raise NotImplementedError


def run(options, args):
    cluster = ThreadManager(4)
    cluster.submitJobAndGetId('sleep 4 && touch job1_done', nbProc=4)
    cluster.submitJobAndGetId('touch job2_done', nbProc=1)

runFromTerminal(__name__, [], run)
