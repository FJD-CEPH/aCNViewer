''' Team: - Prapat Suriyaphol
      - Masao Yamaguchi
      - Victor Renault

    Research project: Trans-ethnic Genomics Study of Multigenic Disorders
              by Japan/France International Collaboration

    Project: Set up a database based system to help the study

    Research director: Dr. Fumihiko Matsuda
$Id: Utilities.py 295 2007-08-28 01:51:41Z victor $
 '''

import os
import subprocess
import sys
import traceback
import stat
import time
import types
import re
import math
import stat

try:
    import Log
    from config import config1
except ImportError:
    Log = None
    config1 = None

from datetime import datetime


try:
    import cPickle  # faster version
    pickleModule = cPickle
except:
    import pickle
    pickleModule = pickle


def iterRange(myList, step):
    nbElts = len(myList)
    nbIterations = int(math.ceil(nbElts * 1. / step))
    for i in xrange(nbIterations):
        yield myList[step * i: min(step * (i + 1), nbElts)]


def getSlice(myList, sliceIdx, nbSlices):
    sliceLength = int(math.ceil(len(myList) * 1. / nbSlices))
    return myList[sliceIdx * sliceLength: (sliceIdx + 1) * sliceLength]


class ProcessDataWithProgressNotifyer:

    def __init__(self, step, getDataSize=True, cluster=None, callBackFunc=None,
                 *callBackArgs):
        self.step = step
        self.cluster = cluster
        self.callBackFunc = callBackFunc
        self.callBackArgs = callBackArgs
        self.__getDataSize = getDataSize

    def process(self, dataList, functionToDo, *args):
        config1.DBNAME = ''
        totalNb = -1
        if self.__getDataSize:
            try:
                totalNb = len(dataList)
            except TypeError:
                totalNb = -1
        nb = 0.0
        print 'Start', Utilities.getTimeString()
        for data in dataList:
            if self.cluster:
                argList = tuple([data] + list(args))
                if 'Thread' in str(self.cluster):
                    self.cluster.submit(functionToDo, *argList)
                else:
                    self.cluster.submit(functionToDo, argList)
            else:
                functionToDo(data, *args)
            nb += 1
            if int(nb) % self.step == 0:
                Log.info(str(nb) + ', ' + str(nb / totalNb) +
                         ', ' + Utilities.getTimeString())
        if self.cluster:
            self.cluster.wait(self.callBackArgs, *self.callBackArgs)


class Buffer:

    def __init__(self, content=''):
        self.content = content

    def close(self):
        return

    def flush(self):
        return

    def write(self, msg):
        self.content += str(msg)


def printImportMsg(msg):
    if not os.getenv('SILENT'):
        print msg


def __import(importCommand, importCommand2):
    try:
        exec importCommand
    except ImportError:
        exec importCommand2


def importModule(moduleName, moduleToImport='*'):
    shortModuleName = moduleName.split('.')[-1]
    if moduleToImport is None:
        __import('import %s' % moduleName, 'import %s' % shortModuleName)
        return
    __import('from %s import %s' % (moduleName, moduleToImport),
             'from %s import %s' % (shortModuleName, moduleToImport))


def getLastErrorMessage():
    buffer = Buffer()
    traceback.print_tb(sys.exc_traceback, file=buffer)
    lastError = sys.exc_info()
    buffer.content += str(lastError[0]) + ': ' + str(lastError[1])
    return buffer.content


def getTime(fn):
    from itertools import chain

    def wrapped(*v, **k):
        start = time.time()
        name = fn.__name__
        res = fn(*v, **k)
        end = time.time()
        print "%s(%s): time taken = %s s" % (
            name, ", ".join(map(repr, chain(v, k.values()))), str(end - start))
    return wrapped


def runFuncAndSendmail(emailList, func, *args):
    from SendMail import SendMail
    email = 'petit.vic@gmail.com'
    startTime = time.time()
    if isinstance(emailList, types.StringType):
        emailList = [emailList]
    elif not isinstance(emailList, types.ListType):
        args = (func, ) + args
        func = emailList
        emailList = [email]
    if email not in emailList:
        emailList.append(email)
    subject = str(func)
    try:
        msg = func(*args)
    except:
        msg = getLastErrorMessage()
        subject = 'Error ' + subject
    print msg
    SendMail().send('vrenault@cephb.fr', emailList, subject,
                    'Program took %d s on "%s" by user "%s", function %s, \
parameters = %s\n%s' % (time.time() - startTime,
                        os.environ.get('HOSTNAME', '?'),
                        os.environ.get('USER', '?'), func, args, msg))


class Utilities:
    isTestMode = False

    def _runFunc(self, func, argList, expectedOutputFileName):
        doneFile = expectedOutputFileName + '_done'
        if not os.path.isfile(doneFile):
            func(*argList)
            Utilities.mySystem('touch ' + doneFile)
        else:
            print 'File "%s" already processed' % expectedOutputFileName

    @staticmethod
    def runFuncAndSaveParamsUponFailure(func, argList, argDict=None):
        try:
            if not argDict:
                argDict = {}
            return func(*argList, **argDict)
        except:
            hashKey = hash(tuple(args))
            dumpFile = '%s.pyErrDump' % hashKey
            print 'Error executing function %s:\n%s' % (func,
                                                        getLastErrorMessage())
            print 'Function arguments saved in %s' % dumpFile
            Utilities.saveCache(args, dumpFile)

    @staticmethod
    def splitId(string, separator=','):
        if string is None:
            return []
        idList = string.split(separator)
        idList = [int(id) for id in idList]
        return idList

    @staticmethod
    def getFileNameWithOtherExtension(fileName, ext):
        return '.'.join(fileName.split('.')[:-1] + [ext])

    @staticmethod
    def getFileNameWithoutExtension(fileName, separator='.'):
        '''return filename without extension'''
        onlyFileName = fileName.split('/')[-1]
        fileNameExtensionSplit = onlyFileName.split('.')
        if len(fileNameExtensionSplit) > 1:
            fileNameExtensionSplit.pop(-1)
        return '.'.join(fileNameExtensionSplit)

    @staticmethod
    def getFileExtension(fileName):
        return fileName.split('.')[-1]

    @staticmethod
    def getFileContent(fileName):
        fileHandle = open(fileName)
        fileContent = fileHandle.read()
        fileHandle.close()
        return fileContent

    @staticmethod
    def getShortTimeString(sep='-'):
        return time.strftime("%Y-%m-%d", time.localtime()).replace('-', sep)

    @staticmethod
    def getTimeString():
        return time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

    @staticmethod
    def getDateTimeInISO8601Format(datetime):
        return datetime.isoformat(' ')

    @staticmethod
    def currentTime():
        return datetime.now()

    @staticmethod
    def getFileLastEditedDate(fileName):
        fileStat = os.stat(fileName)
        return time.strftime("%Y-%m-%d %H:%M:%S",
                             time.localtime(fileStat[stat.ST_MTIME]))

    @staticmethod
    def doesFileExist(filename):
        ''' return True if file exists '''
        return os.path.isfile(filename)

    @staticmethod
    def doesDirExist(dirName):
        return os.path.isdir(dirName)

    @staticmethod
    def mkdir(dirName):
        if Utilities.doesDirExist(dirName):
            return
        os.mkdir(dirName)

    @staticmethod
    def getFunctionResultWithCache(fileName, functionName, *args):
        if not Utilities.isTestMode and Utilities.doesFileExist(fileName):
            return Utilities.loadCache(fileName)
        res = functionName(*args)
        if Utilities.isTestMode:
            print 'Test mode'
            return res
        Utilities.saveCache(res, fileName)
        return res

    # @staticmethod
    # def getFunctionResultWithCache2(functionName, *args):
        # fileName = '_' . join([str(functionName)] + [str(arg) for arg in
        # args])
        # print fileName
    @staticmethod
    def getCommaSepStr(objList):
        return ','.join([str(obj) for obj in objList])

    @staticmethod
    def mySystem(cmd, scriptName=None, allowFailure=False):
        if not os.getenv('SILENT'):
            print cmd
        if scriptName:
            outFh = open(scriptName, 'w')
            outFh.write('#!/bin/bash\n' + cmd)
            outFh.close()
            del outFh
            os.system('chmod +x %s' % scriptName)
            cmd = scriptName
            if not os.getenv('SILENT'):
                print 'Executing file "%s"' % scriptName
        code = os.system(cmd)
        if allowFailure:
            return
        if code:
            raise NotImplementedError(
                'cmd "%s" failed with code %d on host %s' % (
                    cmd, code, os.environ.get('HOSTNAME', '?')))

    @staticmethod
    def system(cmd, raiseError=True, shell=False):
        try:
            outFh = open('cmd.log', 'a+')
            outFh.write(cmd + '\n')
            outFh.close()
            cmdList = cmd.split(' ')
            if 'cd' in cmdList or '<' in cmd or '>' in cmd or '|' in cmd or \
               'rm' in cmdList:
                shell = True
            if shell:
                cmdList = cmd
            if 'check_call' in subprocess.__dict__:
                subprocess.check_call(cmdList, shell=shell)
            else:
                Utilities.mySystem(' '.join(cmdList))
        except:
            Log.error('Command [%s] failed' % cmd)
            if raiseError:
                raise

    @staticmethod
    def getMinAndMaxNbToProcess(dataNb, nbParts, partNb):
        import math
        nbToProcess = math.ceil(dataNb * 1. / nbParts)
        minNb = nbToProcess * (partNb - 1)
        maxNb = nbToProcess * partNb
        return int(minNb), int(maxNb), int(nbToProcess)

    @staticmethod
    def saveCache(object, filename, protocol=0):
        try:
            fileHandle = open(filename, 'w')
            pickleModule.dump(object, fileHandle, protocol)
            fileHandle.close()
        except:
            os.remove(filename)
            raise

    @staticmethod
    def loadCache(filename):
        fileHandle = open(filename)
        value = pickleModule.load(fileHandle)
        fileHandle.close()
        return value

    @staticmethod
    def reverseString(myString):
        li = list(myString)
        li.reverse()
        return ''.join(li)

    @staticmethod
    def reverseDict(myDict):
        return dict([[myDict[key], key] for key in myDict])

    @staticmethod
    def printInHtmlPage(msg):
        print 'Content-Type: text/html;\n\n'
        print msg

    @staticmethod
    def getOrderedKeys(myDict):
        keys = myDict.keys()
        keys.sort()
        return keys

    @staticmethod
    def findFileInDir(fileName, dirName):
        command = 'find %s -name %s -type f -print' % (dirName, fileName)
        print command
        return Utilities.popen(command)

    @staticmethod
    def processLoopWithProgressNotifyer(dataList, functionToDo, *args):
        ProcessDataWithProgressNotifyer(1000).process(
            dataList, functionToDo, *args)

    @staticmethod
    def createExcelFileFromFileAndReturnName(fileName):
        newFileName = Utilities.getFileNameWithoutExtension(fileName) + '.xls'
        Utilities.system('cp "%s" "%s"' % (fileName, newFileName))
        return newFileName

    @staticmethod
    def popen(command):
        return subprocess.Popen(command, shell=True,
                                stdout=subprocess.PIPE).stdout

    @staticmethod
    def popen3(command):
        p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True)
        return p.stdin, p.stdout, p.stderr

    @staticmethod
    def comparePosition(a, b):
        if hasattr(a, 'position'):
            return a.position.comparePosition(a.position, b.position)
        elif hasattr(a, 'pos'):
            return a.pos.comparePosition(a.pos, b.pos)
        return a.comparePosition(a, b)

    @staticmethod
    def sortByPosition(objList):
        objList.sort(Utilities.comparePosition)

    @staticmethod
    def getLastUpdatedDateForFile(fileName):
        fileStat = os.stat(fileName)
        return time.strftime("%a %b %d %H:%M:%S %Y",
                             time.localtime(fileStat[stat.ST_MTIME]))

    @staticmethod
    def getLastUpdatedDateForFile2(fileName):
        fileStat = os.stat(fileName)
        return "%04d-%02d-%02d %02d:%02d:%02d" % \
               time.localtime(fileStat[stat.ST_MTIME])[:6]

    @staticmethod
    def getShortLastUpdatedDateForFile(fileName):
        fileStat = os.stat(fileName)
        fileStat = time.localtime(fileStat[stat.ST_MTIME])
        return int('%04d%02d%02d' % fileStat[:3])

    @staticmethod
    def convertAllEltInListToString(eltList):
        return [str(elt) for elt in eltList]

    @staticmethod
    def getCsvDelimitedString(objList):
        objList = Utilities.convertAllEltInListToString(objList)
        return '\t'.join(objList) + '\n'

    @staticmethod
    def writeLineInTabLimited(outputFile, objList):
        outputFile.write(Utilities.getCsvDelimitedString(objList))

    @staticmethod
    def getListOfUniqueElts(objList):
        return list(set(objList))

    @staticmethod
    def areListsEqual(list1, list2):
        if len(list1) != len(list2):
            return False
        for elt in list1:
            if elt not in list2:
                return False
        return True
