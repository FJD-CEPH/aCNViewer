from ObjectBase import *
from FolderIterator import FolderIterator
from Cluster import guessHpc, HpcScriptBase, HpcBase, ClusterBase, MOABcluster, ThreadManager
from FileHandler import FileNameGetter

IlluminaRun = None

class Cmd:
	def __init__(self, cmd, nbCpus = 1, memory = 4, jobName = None, targetDir = None, jobIdList = None, machineToUseList = None, machineToExcludeList = None,
	             onSubmitAction = False, walltime = None, queue = None, email = None, scriptName = None, optionList = None):
		self.cmd = cmd
		nodeNb = 1
		if type(nbCpus) == types.TupleType:
			nodeNb, nbCpus = nbCpus
		self.nbCpus = nbCpus
		self.nodeNb = nodeNb
		self.memory = memory
		self.jobName = jobName
		self.targetDir = targetDir
		if jobIdList and type(jobIdList) != types.ListType:
			jobIdList = [jobIdList]
		self.jobIdList = jobIdList
		self.machineToUseList = machineToUseList
		self.machineToExcludeList = machineToExcludeList
		self.onSubmitAction = onSubmitAction
		self.walltime = walltime
		self.queue = queue
		self.email = email
		self.scriptName = scriptName
		self.optionList = optionList


class ProcessFileFromCluster:
	_fileExt = None
	_outputFilePattern = None
	_nbCpus = None
	_machineList = None
	_test = False
	_recursive = False
	_maxNbJobs = None
	_laneList = None
	_clusterMachineList = None
	_cpuFactor = 1
	
	def __init__(self, nbCpus = None, binPath = None, fileExt = None, machineList = None, test = False, machineToExcludeList = None, ramPerCpu = None, partNbAndTotalParts = None, dirName = None):
		self._binPath = binPath
		self._nbCpus = nbCpus
		if binPath:
			self._bwaPath = os.path.join(binPath, 'bwa')
		if fileExt:
			self._fileExt = fileExt
		self._machineList = machineList
		self._test = test
		self._ramPerCpu = ramPerCpu
		self._machineToExcludeList = machineToExcludeList
		self._partNbAndTotalParts = partNbAndTotalParts
		self._initFileToProcessList(partNbAndTotalParts, dirName)
			
	def _initFileToProcessList(self, partNbAndTotalParts, dirName):
		if not partNbAndTotalParts:
			self._fileToProcessList = None
			return
		partNb, totalNbParts = partNbAndTotalParts
		self._fileToProcessList = list(FolderIterator().getDataIterator(dirName, self._fileExt))
		totalNbFiles = len(self._fileToProcessList)
		self._fileToProcessList = getSlice(self._fileToProcessList, partNb-1, totalNbParts)
		print 'partNb = %d, totalNbParts = %d, totalNbFiles = %d, nbFilesToProcess = %d' % (partNb, totalNbParts, totalNbFiles, len(self._fileToProcessList))
		
	def _getTargetFileNameForFileAndTargetDir(self, fileName, targetDir):
		sampleName, fcName, uplex = IlluminaRun().getSampleNameRunNameAndUplexFromPath(fileName)
		currentTargetDir = os.path.join(targetDir, sampleName, fcName)
		if uplex:
			currentTargetDir = os.path.join(currentTargetDir, uplex)
		Utilities.mySystem('mkdir -p %s' % currentTargetDir)
		return os.path.join(currentTargetDir, os.path.basename(fileName))
	
	def _createFileLink(self, fileName, targetDir):
		if targetDir:
			targetFileName = self._getTargetFileNameForFileAndTargetDir(fileName, targetDir)
			cmd = 'ln -s %s %s' % (fileName, targetFileName)
			fileName = targetFileName
			os.system(cmd)
		return fileName
	
	def _isFileToProcess(self, fileName):
		if self._fileToProcessList:
			print 'KOKOJO data: %d' % len(self._fileToProcessList)
			print fileName in self._fileToProcessList
			return fileName in self._fileToProcessList
		return True
	
	def _getCmdToRunFromFile(self, fileName, *args):
		raise NotImplementedError
	
	def __getFileWithExtensionInDir(self, dirName, fileExt):
		if os.path.sep in fileExt:
			return fileExt
		elif '.' in fileExt:
			return os.path.join(dirName, fileExt)
		fileList = glob.glob(os.path.join(dirName, '*%s' % fileExt))
		if len(fileList) != 1:
			raise NotImplementedError('Expecting one file in %s with ext %s but found %s' % (dirName, fileExt, fileList))
		return fileList[0]
	
	def _onSubmitAction(self, jobId, inputFile, outputFile, cmd):
		return
	
	def _replaceInputAndOutputInStr(self, cmd, inputFile, outputFile):
		if type(inputFile) != types.ListType:
			cmd = cmd.replace('[input]', inputFile)
		return cmd.replace('[output]', outputFile)
	
	def __createCmdListDumpFileAndGetName(self, fileName, cmdList, idxToProcess):
		if type(fileName) == types.ListType:
			dirName = os.path.dirname(fileName[0])
			fileName = os.path.join(dirName, '_'.join([os.path.basename(currentFile) for currentFile in fileName]))
		dumpFile = fileName + '_cmdList_%s.pyDump' % Utilities.getTimeString()
		cmdList = [[cmd, False] for cmd in cmdList]
		for i in range(idxToProcess):
			cmdList[i][-1] = True
		Utilities.saveCache(cmdList, dumpFile)
		return dumpFile
	
	def __getCmdFromParamList(self, paramList):
		currentCmd = paramList[-1]
		if type(currentCmd) != types.StringType:
			currentCmd = currentCmd.cmd
		return currentCmd
	
	def _resumeRecursiveCmd(self, dumpFileName, lastIdx = None, test = False):
		self._recursive = True
		#time.sleep(15)
		cmdList = Utilities.loadCache(dumpFileName)
		if lastIdx is None:
			for i, cmd in enumerate(cmdList):
				if not cmd[-1]:
					nextIdx = i
					break
		else:
			cmdList[lastIdx][1] = True
			nextIdx = lastIdx + 1
		Utilities.saveCache(cmdList, dumpFileName)
		if len(cmdList) == nextIdx:
			print 'Finished recursive cmd'
			return
		while True:
			nextCmdTuple = cmdList[nextIdx]
			if not nextCmdTuple[1]:
				paramList = list(nextCmdTuple[0])
				break
			nextIdx += 1
			print 'Passing done job: [%s]' % self.__getCmdFromParamList(list(nextCmdTuple[0]))
		if test:
			print nextIdx, paramList
			currentCmd = self.__getCmdFromParamList(paramList)
			print nextIdx, currentCmd
		cmd, fileName, outputFile, onSubmitAction, walltime, queue, machineToUseList, machineToExcludeList, nbCpus, nodeNb, memory, jobName, expectedFormat, otherExpectedFormatList, errorFile, outFile, jobId, dumpFileName2, email, scriptName, optionList = \
		   self._getSubmitParamListFromCmd(*[paramList[0]] + paramList + [dumpFileName, None, None, nextIdx])
		if test:
			print cmd
			return
		cluster = guessHpc()
		jobId = cluster.submitJobAndGetId(cmd, errorFile = errorFile, outputFile = outFile, dependentJobId = jobId, nbProc = nbCpus, jobName = jobName, memory = memory, machineToUseList = machineToUseList, node = nodeNb,
		                                  machineToExcludeList = machineToExcludeList, walltime = walltime, queue = queue, scriptName = scriptName, optionList = optionList)
	
	def _getSubmitParamListFromCmd(self, fileName, expectedFormat, outputFilePattern, cmd, dumpFileName = None, fileToJobIdDict = None, cmdList = None, i = None):
		if fileToJobIdDict is None:
			fileToJobIdDict = {}
		baseFile = onSubmitAction = walltime = queue = email = None
		machineToUseList = self._machineList
		machineToExcludeList = self._machineToExcludeList
		nbCpus = nodeNb = 1
		memory = 4
		jobName = otherExpectedFormatList = None
		if type(expectedFormat) == types.TupleType:
			expectedFormat, baseExt = expectedFormat
			baseFile = self.__getFileWithExtensionInDir(os.path.dirname(fileName), baseExt)
		elif type(expectedFormat) == types.ListType:
			otherExpectedFormatList = expectedFormat[1:]
			expectedFormat = expectedFormat[0]
		print 'F', fileName, expectedFormat#, myStr(cmd)
		#print cmdList
		if type(fileName) != types.ListType and expectedFormat != Utilities.getFileExtension(fileName):
			fileName = self.__getFileWithExtensionInDir(os.path.dirname(fileName), expectedFormat)
		outputFile = outputFilePattern
		if '%s' in outputFilePattern:
			outputFile = outputFilePattern % fileName
		doneFile = outputFile + '_done'
		if os.path.isfile(doneFile):
			print 'Already processed file: %s' % fileName
			if isinstance(cmd, Cmd) and cmd.onSubmitAction:
				self._onSubmitAction(None, fileName, outputFile, cmd)
			return outputFile
		cluster = guessHpc()
		logSuffix = ''
		if not jobName and type(cmd) != types.StringType:
			jobName = cmd.jobName
		if cluster.TASK_ID_VAR and ((jobName and '[' in jobName) or (type(cmd) != types.StringType and cmd.scriptName and '[' in cmd.scriptName)):
			logSuffix = '_${%s}' % cluster.TASK_ID_VAR
		errorFile = outputFile + '%s.err' % logSuffix
		outFile = outputFile + '%s.out' % logSuffix
		cmdJobIdList = []
		scriptName = optionList = None
		if type(cmd) != types.StringType:
			nbCpus = cmd.nbCpus
			jobName = cmd.jobName
			memory = cmd.memory
			if cmd.targetDir:
				errorFile = os.path.join(cmd.targetDir, os.path.basename(errorFile))
				outFile = os.path.join(cmd.targetDir, os.path.basename(outFile))
			if cmd.jobIdList:
				cmdJobIdList = cmd.jobIdList
			if cmd.machineToUseList:
				machineToUseList = cmd.machineToUseList
			if cmd.machineToExcludeList:
				machineToExcludeList = cmd.machineToExcludeList
			if cmd.walltime:
				walltime = cmd.walltime
			if cmd.queue:
				queue = cmd.queue
			if cmd.email:
				email = cmd.email
			if cmd.scriptName:
				scriptName = cmd.scriptName
			if cmd.optionList:
				optionList = cmd.optionList
			nodeNb = cmd.nodeNb
			onSubmitAction = cmd.onSubmitAction
			cmd = cmd.cmd
		if jobName and '[' in jobName and not scriptName:
			#scriptName = FileNameGetter(outputFile).get('_cmd.sh')
			scriptName = outputFile + '_cmd.sh'
			cmd += ' 1> %s 2> %s' % (outFile, errorFile)
			errorFile = outFile = None
			print 'New cmd = [%s]' % cmd
		cmd = self._replaceInputAndOutputInStr(cmd, fileName, outputFile)
		if baseFile:
			cmd = cmd.replace('[base]', baseFile)
		if cmd[-1] == ';':
			cmd = cmd[:-1]
		if logSuffix:
			doneFile = doneFile.replace('_done', '%s_done' % logSuffix)
		cmd = cmd + ' && touch %s' % doneFile
		#print cmd
		#if isinstance(cluster, HpcBase):
		if type(fileName) == types.StringType:
			jobId = fileToJobIdDict.get(fileName)
		else:
			jobId = [fileToJobIdDict.get(currentFileName) for currentFileName in fileName if fileToJobIdDict.has_key(currentFileName)]
		if not jobId:
			jobId = []
		if type(jobId) != types.ListType:
			jobId = [jobId]
		if otherExpectedFormatList:
			jobId += [fileToJobIdDict.get(otherExpectedFormat) for otherExpectedFormat in otherExpectedFormatList if fileToJobIdDict.has_key(otherExpectedFormat)]
		jobId += cmdJobIdList
		if jobId:
			otherExpectedFormatStr = ''
			if otherExpectedFormatList:
				otherExpectedFormatStr += ' && ' + ' && '.join(['[ -f %s_done ]' % otherExpectedFormat for otherExpectedFormat in otherExpectedFormatList])
			cmd = 'if [ -f %s_done ]%s; then %s; fi' % (expectedFormat, otherExpectedFormatStr, cmd)
		if self._recursive:
			if not dumpFileName:
				dumpFileName = self.__createCmdListDumpFileAndGetName(fileName, cmdList, i)
			cmd += ' && python %s -p recursiveCmd -f %s -i %d' % (os.path.abspath(__file__), dumpFileName, i)
		return cmd, fileName, outputFile, onSubmitAction, walltime, queue, machineToUseList, machineToExcludeList, nbCpus, nodeNb, memory, jobName, expectedFormat, otherExpectedFormatList, \
		       errorFile, outFile, jobId, dumpFileName, email, scriptName, optionList
	
	def _runCmdList(self, cmdList, fileName, cluster = None, jobId = None, fileToJobIdDict = None, targetDumpFile = None):
		if not fileToJobIdDict:
			fileToJobIdDict = {}
		if type(fileName) == types.TupleType:
			if fileName[1]:
				fileName = fileName[1]
			else:
				fileName = fileName[-1]
		if jobId:
			fileToJobIdDict[cmdList[0][0]] = jobId
		jobNb = 0
		newCmdList = []
		for i, (expectedFormat, outputFilePattern, cmd) in enumerate(cmdList):
			paramList = self._getSubmitParamListFromCmd(fileName, expectedFormat, outputFilePattern, cmd, fileToJobIdDict = fileToJobIdDict, cmdList = cmdList, i = i)
			if type(paramList) == types.StringType:
				fileName = paramList
				continue
			cmd, fileName, outputFile, onSubmitAction, walltime, queue, machineToUseList, machineToExcludeList, nbCpus, nodeNb, memory, jobName, expectedFormat, otherExpectedFormatList, errorFile, outFile, jobId, dumpFileName, email, scriptName, optionList = paramList
			#if not jobId:
				#jobId = fileToJobIdDict.get(os.path.basename(fileName))
			print 'FILENAME %s, jobId %s' % (fileName, jobId), fileToJobIdDict
			if self._test:
				jobNb += 1
				jobId = jobNb
				print 'cmd = [%s]' % cmd
				continue
			if self._ramPerCpu:
				nbCpus = max(nbCpus, int(math.ceil(1. * memory / self._ramPerCpu)))
			if self._clusterMachineList:
				self._simulateCluster(cmd, nbCpus, memory)
				continue
			if targetDumpFile:
				cmd += ' 1> %s 2> %s' % (outFile, errorFile)
				newCmdList.append((cmd, outFile))
			else:
				#print 'NB cpus = %d, memory = %d' % (nbCpus, memory)
				jobId = cluster.submitJobAndGetId(cmd, errorFile = errorFile, outputFile = outFile, dependentJobId = jobId, nbProc = nbCpus, jobName = jobName, memory = memory, machineToUseList = machineToUseList, node = nodeNb,
					                                  machineToExcludeList = machineToExcludeList, walltime = walltime, queue = queue, email = email, scriptName = scriptName, expectedOutputFile = outputFile, optionList = optionList)
			fileToJobIdDict[outputFile] = jobId
			if onSubmitAction:
				self._onSubmitAction(jobId, fileName, outputFile, cmd)
			#else:
				#Utilities.mySystem(cmd)
				#self._onSubmitAction(None, fileName, outputFile)
			if self._recursive:
				break
			fileName = outputFile
		if targetDumpFile:
			print 'runCmdInParallel with %d threads' % self._nbCpus
			Utilities.saveCache(newCmdList, targetDumpFile)
			cmd = 'python %s -c %s -n %d' % (os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runCmdInParallel.py'), targetDumpFile, int(self._nbCpus * self._cpuFactor))
			jobId = cluster.submitJobAndGetId(cmd, errorFile = FileNameGetter(targetDumpFile).get('err'), outputFile = FileNameGetter(targetDumpFile).get('out'), nbProc = self._nbCpus, jobName = jobName, memory = memory, machineToUseList = machineToUseList, node = nodeNb,
			                                  machineToExcludeList = machineToExcludeList, walltime = walltime, queue = queue, email = email, scriptName = scriptName, expectedOutputFile = outputFile, optionList = optionList)
			for cmd, outputFile in newCmdList:
				fileToJobIdDict[outputFile] = jobId
		return fileToJobIdDict
	
	def _finishProcess(self, cluster, *args):
		return
	
	def __getFileSizeAndDumpFileFromCmd(self, cmd):
		keyword = 'recursiveCmd -f '
		if keyword not in cmd:
			raise NotImplementedError('Keyword "%s" not in cmd [%s]' % (keyword, cmd))
		dumpFile = cmd.split(keyword)[-1].split()[0]
		cmdList = Utilities.loadCache(dumpFile)
		for cmd in cmdList:
			cmd = cmd[0][2]
			if type(cmd) != types.StringType:
				cmd = cmd.cmd
			if 'bwa ' in cmd:
				partList = cmd.split()
				fileSize = 0
				for fileName in partList[5:7]:
					fileSize += SpaceManager().getFileSize(fileName, True)
				return fileSize, dumpFile
	
	def __getJobListSortedByFileSize(self, jobList):
		fileSizeAndJobList = []
		for cmd, outputFile in jobList:
			fileSize, dumpFile = self.__getFileSizeAndDumpFileFromCmd(cmd)
			fileSizeAndJobList.append((fileSize, dumpFile, cmd, outputFile))
		fileSizeAndJobList.sort()
		fileSizeAndJobList.reverse()
		return fileSizeAndJobList
	
	def __runMergedJobs(self, jobList, cluster, dirName):
		nbJobsToMerge = int(math.ceil(1. * len(jobList) / self._maxNbJobs))
		jobList = self.__getJobListSortedByFileSize(jobList)
		step = int(math.ceil(len(jobList)) / nbJobsToMerge)
		jobToMergeDict = defaultdict(list)
		for i in range(step):
			for j in range(nbJobsToMerge):
				idx = i * nbJobsToMerge + j
				jobToMergeDict[i].append(jobList[idx])
		print '%d jobs to merge with an average of %d lanes per job' % (len(jobToMergeDict), nbJobsToMerge)
		for i in range(step):
			targetDumpFile = os.path.join(dirName, 'mergedCmd%d.pyDump' % i)
			currentJobList = jobToMergeDict[i]
			cmdList = []
			for fileSize, dumpFile, cmd, outputFile in currentJobList:
				cmdList += Utilities.loadCache(dumpFile)
			Utilities.saveCache(cmdList, targetDumpFile)
			cmd = 'python %s -p recursiveCmd -f %s -i %d' % (os.path.abspath(__file__), targetDumpFile, idx)
			cluster.submitJobAndGetId(cmd, errorFile = targetDumpFile + '.err', outputFile = targetDumpFile + '.out', machineToUseList = self._machineList)
			
	def process(self, dirName, *args):
		#if self._machineList:
			#cluster = RemoteThreadManager(self._machineList)
		#else:
		cluster = guessHpc(nbCpus = self._nbCpus)
		print cluster
		roundNb = 1
		hasJobBeenSubmitted = False
		jobList = []
		while True:
			print 'Round nb %d' % roundNb
			for fileName in FolderIterator().getDataIterator(dirName, self._fileExt):
				print 'C = [%s]' % str(fileName)
				if not self._isFileToProcess(fileName):
					print 'Passing file "%s"' % str(fileName)
					continue
				if self._laneList:
					laneNb, seqFile1, seqFile2 = fileName
					sampleName, fcName, uplex = IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile1)
					laneNb, readNb = IlluminaRun().getLaneAndReadNbfromFile(seqFile1)
					print fileName, sampleName, fcName, laneNb
					if (sampleName, fcName, laneNb) not in self._laneList and sampleName not in self._laneList:
						print 'Passing file %s' % str(fileName)
						continue
				#self._runCmdList(self, cmdList, fileName, cluster = None)
				print 'Processing', fileName
				cmd = self._getCmdToRunFromFile(fileName, *args)
				print 'CMD', cmd
				if type(cmd) == types.StringType:
					outputFile = self._outputFilePattern % fileName
					doneFile = outputFile + '_done'
					if os.path.isfile(doneFile):
						print 'Already processed file: %s' % fileName
						continue
					cmd += ' && touch %s' % doneFile
					if self._recursive and self._maxNbJobs:
						jobList.append((cmd, outputFile))
						continue
					cluster.submitJobAndGetId(cmd, errorFile = outputFile + '.err', outputFile = outputFile + '.out', machineToUseList = self._machineList)
				else:
					if self._recursive and self._maxNbJobs:
						jobList.append((cmd, outputFile))
						continue
					if isinstance(cluster, ClusterBase):
						cluster.submit(self._runCmdList, (cmd, fileName))
					else:
						self._runCmdList(cmd, fileName, cluster)
						#break
				hasJobBeenSubmitted = True
			if isinstance(cluster, HpcBase):
				break
			if not hasJobBeenSubmitted or not cluster.getNbRunningJobs():
				print 'Done'
				break
			roundNb += 1
			if not isinstance(cluster, MOABcluster):
				break
			cluster.waitUntilClusterIsFree()
		if self._recursive and self._maxNbJobs:
			self.__runMergedJobs(jobList, cluster, dirName)
		if isinstance(cluster, ThreadManager):
			cluster.wait()
		self._finishProcess(cluster, dirName, *args)

			
