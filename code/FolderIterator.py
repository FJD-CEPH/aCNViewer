import glob
import os
import copy
import types

from Utilities import Utilities
from ValueParser import ValueParser

try:
	from IlluminaRun import IlluminaRun
except ImportError:
	IlluminaRun = None

class FolderIterator:
	def __init__(self, followLinks = True, isDir = False):
		self.followLinks = followLinks
		self._isDir = isDir
		
	def _updateFileListFromDirNameAndFileExt(self, dirName, fileList, isFileToProcessFunc, *args):
		#print 'd', dirName
		for fileName in glob.glob(os.path.join(dirName, '*')):
			#print 'fi', fileName
			if not self._isDir and os.path.isdir(fileName):
				if not self.followLinks and os.path.islink(fileName):
					print 'Passing %s' % fileName
					continue
				self._updateFileListFromDirNameAndFileExt(fileName, fileList, isFileToProcessFunc, *args)
			elif self._isDir and not isFileToProcessFunc(fileName, *args):
				if not self.followLinks and os.path.islink(fileName):
					print 'Passing %s' % fileName
					continue
				self._updateFileListFromDirNameAndFileExt(fileName, fileList, isFileToProcessFunc, *args)
			elif isFileToProcessFunc(fileName, *args):
				fileList.append(fileName)
				
	def __doesFileHasExtensionInList(self, fileName, fileExtList):
		for fileExt in fileExtList:
			if os.path.basename(fileName)[-len(fileExt):] == fileExt:
				return True
				
	def __doesFileHasExtension(self, fileName, fileExt):
		return (not self._isDir and os.path.isfile(fileName) and ((type(fileExt) == types.ListType and self.__doesFileHasExtensionInList(fileName, fileExt)) \
		                                     or ('.' not in fileExt and Utilities.getFileExtension(fileName) == fileExt or fileExt == '*') or \
		                                     (os.path.basename(fileName)[-len(fileExt):] == fileExt))) or \
		       (self._isDir and os.path.basename(fileName)[-len(fileExt):] == fileExt)
				
	def getFileListInDirWithExt(self, dirName, fileExt, cacheFileName = None):
		if cacheFileName:
			return Utilities.getFunctionResultWithCache(cacheFileName, self._getFileList, dirName, self.__doesFileHasExtension, fileExt)
		return self._getFileList(dirName, self.__doesFileHasExtension, fileExt)
	
	def _getFileList(self, dirName, isFileToProcessFunc, *args):
		#print 'DIR', dirName
		fileList = []
		if not self._isDir and os.path.isfile(dirName):
			if isFileToProcessFunc(dirName, *args):
				fileList.append(dirName)
		elif self._isDir and os.path.isdir(dirName):
			if isFileToProcessFunc(dirName, *args):
				fileList.append(dirName)
			else:
				self._updateFileListFromDirNameAndFileExt(dirName, fileList, isFileToProcessFunc, *args)
		else:
			if not self.followLinks and os.path.islink(dirName):
				print 'Passing %s' % dirName
			else:
				self._updateFileListFromDirNameAndFileExt(dirName, fileList, isFileToProcessFunc, *args)
		return fileList
	
	def __isSeqFile(self, fileName, allowFailed):
		seqEnd = '_sequence.txt.gz'
		seqEnd2 = '_sequence.txt.bz2'
		seqEnd3 = '.fastq.gz'
		isFileToProcess = os.path.basename(fileName)[-len(seqEnd):] == seqEnd or os.path.basename(fileName)[-len(seqEnd2):] == seqEnd2 or \
		                os.path.basename(fileName)[-len(seqEnd3):] == seqEnd3
		if not allowFailed and 'BAD' in fileName:
			isFileToProcess = False
		return isFileToProcess
	
	def __doesSeqFileNameContainReadNb(self, seqFile):
		#return ValueParser.isNb(os.path.basename(seqFile).split('_')[2])
		laneNb, readNb = IlluminaRun().getLaneAndReadNbfromFile(seqFile, False)
		return readNb is not None
	
	def __getBaseNameForFile(self, fileName):
		fileName = os.path.basename(fileName).split('.')[0]
		partList = fileName.split('_')
		readIdx = None
		if '_L00' in fileName and '_R' in fileName:
			for i, part in enumerate(partList):
				if len(part) == 2 and part[0] == 'R':
					readIdx = i
		else:
			laneIdx = None
			for i, part in enumerate(partList):
				if len(part) == 1 and ValueParser().isNb(part):
					if laneIdx:
						readIdx = i
						break
					else:
						laneIdx = i
		if not readIdx:
			raise NotImplementedError('Could not find readIdx for file %s' % fileName)
		return '_'.join(partList[:readIdx] + partList[readIdx+1:])
	
	def __getSortedSeqFileList(self, seqFileList):
		newSeqFileList = []
		for seqFile in seqFileList:
			sampleName, fcName, uplex = IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile)
			laneNb, readNb = IlluminaRun().getLaneAndReadNbfromFile(seqFile)
			newSeqFileList.append((sampleName, fcName, laneNb, readNb, seqFile))
		newSeqFileList.sort()
		return [seqInfo[-1] for seqInfo in newSeqFileList]
	
	def _iterSeqFiles(self, seqFileList, allowReadNbToBeTheSame = False, sort = True):
		if sort:
			seqFileList.sort()
		seqFileList = self.__getSortedSeqFileList(seqFileList)
		#print 'SS', seqFileList
		seqFileList = copy.deepcopy(seqFileList)
		prevSeqFile = seqFileList.pop(0)
		laneNb, readNb = IlluminaRun().getLaneAndReadNbfromFile(prevSeqFile)
		prevLaneAndReadNb = IlluminaRun().getLaneAndReadNbfromFile(prevSeqFile)
		prevSampleInfo = IlluminaRun().getSampleNameRunNameAndUplexFromPath(prevSeqFile)
		lastReturnedSeqFile = None
		if not self.__doesSeqFileNameContainReadNb(prevSeqFile) and allowReadNbToBeTheSame:
			lastReturnedSeqFile = prevSeqFile
			yield prevLaneAndReadNb[0], prevSeqFile
		isOldFormattedSeqName = False
		if 'sequence' in os.path.basename(prevSeqFile):
			isOldFormattedSeqName = True
		for seqFile in seqFileList:
			laneNb, readNb = IlluminaRun().getLaneAndReadNbfromFile(seqFile)
			isSingleEnd = False
			currentSampleInfo = IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile)
			#print '))', currentSampleInfo, laneNb, readNb
			if currentSampleInfo[0] == prevSampleInfo[0] and laneNb == prevLaneAndReadNb[0]:
				if not allowReadNbToBeTheSame and readNb == prevLaneAndReadNb[1]:
					raise NotImplementedError('Read nb should be different for seqFile %s' % seqFile)
				#print laneNb, allowReadNbToBeTheSame, prevSeqFile, seqFile,lastReturnedSeqFile
				if allowReadNbToBeTheSame and os.path.dirname(prevSeqFile) != os.path.dirname(seqFile):
					if lastReturnedSeqFile != prevSeqFile:
						lastReturnedSeqFile = prevSeqFile
						yield laneNb, prevSeqFile
				else:
					if isOldFormattedSeqName or self.__getBaseNameForFile(prevSeqFile) == self.__getBaseNameForFile(seqFile):
						lastReturnedSeqFile = seqFile
						prevSampleInfo = IlluminaRun().getSampleNameRunNameAndUplexFromPath(prevSeqFile)
						sampleInfo = IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile)
						if prevSampleInfo != sampleInfo:
							raise NotImplementedError('Sample information do not match for files %s and %s: %s != %s' % (prevSeqFile, seqFile, prevSampleInfo, sampleInfo))
						yield laneNb, prevSeqFile, seqFile
			else:
				if not self.__doesSeqFileNameContainReadNb(seqFile):
					isSingleEnd = True
					lastReturnedSeqFile = seqFile
					yield laneNb, seqFile
			prevSeqFile = seqFile
			prevLaneAndReadNb = laneNb, readNb
			prevSampleInfo = currentSampleInfo
		if allowReadNbToBeTheSame and prevSeqFile != lastReturnedSeqFile:
			yield laneNb, prevSeqFile
		#if isSingleEnd:
			#yield laneNb, prevSeqFile
	
	def __trimmedSeqIterator(self, seqFileList):
		while len(seqFileList):
			seqFile = seqFileList.pop(0)
			seqFile2 = seqFileList.pop(0)
			if os.path.dirname(seqFile) != os.path.dirname(seqFile2):
				if 'single.fq.gz' in seqFile:
					seqFileList.insert(0, seqFile2)
					singleSeqFile = seqFile
					seqFile = seqFile2 = None
				else:
					raise NotImplementedError('%s and %s should be paired' % (seqFile, seqFile2))
			else:
				singleSeqFile = None
				if seqFileList:
					singleSeqFile = seqFileList.pop(0)
					if 'single.fq.gz' not in singleSeqFile:
						seqFileList.insert(0, singleSeqFile)
						singleSeqFile = None
			if singleSeqFile:
				yield seqFile, seqFile2, singleSeqFile
			else:
				yield seqFile, seqFile2
	
	def getDataIterator(self, seqDir, fileExt, isSeq = False):
		print 'seqDir, fileExt', seqDir, fileExt
		if fileExt:
			iterator = self.getFileListInDirWithExt(seqDir, fileExt)
			if isSeq:
				iterator = self._iterSeqFiles(iterator, True, False)
		else:
			fh = os.popen('find %s -name "*.fq.gz" -follow | head -1' % seqDir)
			if fh.read().strip():
				seqFileList = list(self.getFileListInDirWithExt(seqDir, '.fq.gz'))
				seqFileList.sort()
				iterator = self.__trimmedSeqIterator(seqFileList)
			else:
				seqFileList = self.getSeqFileListInDir(seqDir)
				seqFileList = [(IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile), IlluminaRun().getLaneAndReadNbfromFile(seqFile), seqFile) for seqFile in seqFileList]
				seqFileList.sort()
				seqFileList = [seqFile for info1, info2, seqFile in seqFileList]
				if seqFileList:
					dirName = os.path.basename(os.path.dirname(seqFileList[0]))
					if ValueParser().isNb(dirName):
						seqFileList = [(IlluminaRun().getSampleNameRunNameAndUplexFromPath(seqFile), IlluminaRun().getSampleNbFromFile(seqFile), IlluminaRun().getLaneAndReadNbfromFile(seqFile), seqFile) for seqFile in seqFileList]
						seqFileList.sort()
						seqFileList = [seqFile for info1, info2, info3, seqFile in seqFileList]
				print '%d files' % len(seqFileList)
				print '\n'.join(seqFileList)
				iterator = self._iterSeqFiles(seqFileList, True, False)
		return iterator
	
	def getSeqFileListInDir(self, dirName, allowFailed = False):
		return self._getFileList(dirName, self.__isSeqFile, allowFailed)
	
	def getSymFileListInDir(self, dirName):
		return self._getFileList(dirName, self.__isSymLink)
	
	def __isSymLink(self, fileName):
		return os.path.islink(fileName)
	
	def __getFileListFromRemoteServer(self, srcDir, ssh, withinMinutes):
		hostname, srcDir = srcDir.split(':')
		baseName = os.path.basename(srcDir)
		originalSrcDir = srcDir
		fileExtList = []
		if '*' in baseName:
			baseName = baseName.split('*')[-1]
			print 'B=[%s]' % baseName
			srcDir = os.path.dirname(srcDir)
			if '{' in baseName:
				from WebExtractor import WebExtractor
				baseName = WebExtractor()._getStrIncludedInTag(baseName, '{', '}')
				print 'SOKOK', baseName
				fileExtList = baseName.split(',')
			else:
				fileExtList.append(baseName)
		if fileExtList in [[''], []] and '*' in srcDir:
			print [hostname, srcDir]
			stdin, stdout, stderr = ssh.execCommand('ls %s' % originalSrcDir)
			fileList = ['%s:%s' % (hostname, fileName) for fileName in stdout.read().strip().split('\n')]
			print 'FILES', fileList
		else:
			print 'Getting file list from %s in %s with ext %s' % (hostname, srcDir, fileExtList)
			cmd = 'find %s -type f -follow' % (srcDir + os.path.sep)
			if fileExtList:
				cmd += ' %s' % ' -o '.join(['-name "*.%s"' % fileExt for fileExt in fileExtList])
			if withinMinutes:
				cmd += ' -mmin -%d' % withinMinutes
			print cmd
			stdin, stdout, stderr = ssh.execCommand(cmd)
			print 'Exec'
			fileList = stdout.read().strip().split('\n')
			print 'Read from stdout'
			fileList = ['%s:%s' % (hostname, fileName) for fileName in fileList if not fileExtList or (fileExtList and self.__doesFileHasExtension(fileName, fileExtList))]
		print '%d files' % len(fileList)
		return fileList
		
	def getFileListInDir(self, srcDir, ssh = None, withinMinutes = None):
		hostname = None
		if ':' in srcDir:
			fileList = self.__getFileListFromRemoteServer(srcDir, ssh, withinMinutes)
		else:
			if os.path.isdir(srcDir):
				fileList = glob.glob(os.path.join(srcDir, '*'))
			else:
				fileList = glob.glob(srcDir)
		for fileName in fileList:
			if hostname:
				yield fileName
			else:
				if os.path.isfile(fileName) or ':' in fileName:
					yield fileName
				else:
					for fileName2 in self._getFileListToCopy(fileName):
						yield fileName2
	

from ObjectBase import *

def run(options, args):
	for f in FolderIterator().getDataIterator(options.dirName, options.pattern):
		print f
		if len(f) != 3:
			raise NotImplementedError

runFromTerminal(__name__, [CommandParameter('d', 'dirName', 'string'),
                           CommandParameter('p', 'pattern', 'string')], run)
						