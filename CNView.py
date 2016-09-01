import math
import types
import copy
import os
import glob

from MultiDict import *
from FileHandler import *
from commandLineUtils import *
from structure import Position
from Cluster import ThreadManager, _getNbAvailableCpus, guessHpc
from Utilities import Utilities
#from stats import getMedianFromList
from annovar import ProcessFileFromCluster, Cmd
from R import R
from color import Color
from coverage import Coverage, OrientedPosition

class DefaultPloidyDict(dict):
	def has_key(self, key):
		return True
	
	def __getitem__(self, key):
		return 2


class RunAscat:
	def __init__(self, binDir = None):
		self.__binDir = binDir
	
	def __createGroupFilesAndGetHeaderAndPloidyIdx(self, ploidyFile, targetDir):
		fh = ReadFileAtOnceParser(ploidyFile)
		header = [columnName.lower() for columnName in fh.getSplittedLine()]
		ploidyIdx = header.index('ploidy')
		outFhDict = {}
		groupSetDict = defaultdict(set)
		for columnName in header[ploidyIdx+1:]:
			outFileName = os.path.join(targetDir, os.path.basename(FileNameGetter(ploidyFile).get('_%s.txt' % columnName.replace(' ', ''))))
			outFh = CsvFileWriter(outFileName)
			outFh.write(['sample', 'ploidy', 'group'])
			outFhDict[columnName] = outFh
		for splittedLine in fh:
			for i, columnName in enumerate(header[ploidyIdx+1:]):
				outFh = outFhDict[columnName]
				group = splittedLine[ploidyIdx+1+i]
				groupSetDict[columnName].add(group)
				outFh.write([splittedLine[0], splittedLine[ploidyIdx], group])
		for columnName in header[ploidyIdx+1:]:
			print 'Column "%s" has %d values' % (columnName, len(groupSetDict[columnName]))
		return header, ploidyIdx, groupSetDict, dict([[columnName, outFh._fileName] for columnName, outFh in outFhDict.iteritems()])
		
	def _createDendrogramForEachFeature(self, ploidyFile, targetDir, windowSize, percentage, segmentFile, chrFile):
		Utilities.mySystem('mkdir -p %s' % targetDir)
		header, ploidyIdx, groupSetDict, outFhDict = Utilities.getFunctionResultWithCache(FileNameGetter(os.path.join(targetDir, os.path.basename(ploidyFile))).get('pyDump'), self.__createGroupFilesAndGetHeaderAndPloidyIdx, ploidyFile, targetDir)
		# filter segment file to keep only samples with BCLC Staging information
		# add option in createDendrogram to replace null value in group column with a default value
		# add option in createDendrogram to replace sample name with group value
		for columnName in header[ploidyIdx+1:]:
			print 'Plotting dendrogram for column "%s"' % columnName
			groupSet = groupSetDict[columnName]
			useShape = False
			if len(groupSet) <= 5:
				useShape = True
			#obj = CNView(windowSize, percentage, self.__binDir, useShape)
			#groupDict, group2Dict, colorDict, shapeDict, shapeDict2 = obj._getGroupAndColorCodeDictFromFile(ploidyFile)
			#obj._createDendrogram(matrixFile, groupDict, colorDict, shapeDict, ploidyFile2, shapeDict2, coeff, keyword = columnName.replace(' ', ''))
			CNView(windowSize, percentage, self.__binDir, useShape).process(segmentFile, chrFile, targetDir, outFhDict[columnName], histogram = 0, merge = 1, dendrogram = 1, groupColumnName = 'test', keyword = columnName.replace(' ', ''), defaultGroupValue = 'UNKNOWN')
			#break
	
	def __getSampleDictFromFh(self, fh, sampleIdx):
		sampleDict = {}
		for splittedLine in fh:
			sampleDict[splittedLine[sampleIdx]] = splittedLine
		return sampleDict
	
	def _mergePloidyFileWithSampleInfoFile(self, ploidyFile, sampleFile, sampleAliasFile):
		sampleDict = self._getSampleDictFromFile(sampleAliasFile)
		outFh = CsvFileWriter(FileNameGetter(ploidyFile).get('_merged.txt'))
		fh = ReadFileAtOnceParser(ploidyFile)
		sampleHeader, sampleFh, sampleIdx = self._getHeaderFhAndSampleIdxForSampleFile(sampleFile)
		sampleInfoDict = self.__getSampleDictFromFh(sampleFh, sampleIdx)
		header = fh.getSplittedLine() + sampleHeader
		outFh.write(header)
		for splittedLine in fh:
			sampleName = splittedLine[sampleIdx]
			sampleAlias = sampleDict[sampleName]
			outFh.write(splittedLine + sampleInfoDict[sampleAlias])
		
	def __getBuildFromPfbFile(self, fileName):
		for part in os.path.basename(fileName).split('.'):
			if part[:2] == 'hg':
				return part
		
	def __getProbeToPosDictFromFile(self, fileName):
		probeDict = {}
		for splittedLine in ReadFileAtOnceParser(fileName):
			chrName, start, end, probeName = splittedLine
			probeDict[probeName] = chrName.replace('chr', ''), int(end)
		return probeDict
	
	def _createFileWithUpdatedPositions(self, fileName, rawProfileFile, targetBuild = None):
		fileExt = Utilities.getFileExtension(fileName)
		if not targetBuild:
			targetBuild = self.__getBuildFromPfbFile(fileName)
		liftedOverFile = Utilities.getFunctionResultWithCache(FileNameGetter(rawProfileFile).get('_%s.pyDump' % targetBuild), self._liftOverRawProbeFile, rawProfileFile, targetBuild)
		probeDict = Utilities.getFunctionResultWithCache(FileNameGetter(liftedOverFile).get('pyDump2'), self.__getProbeToPosDictFromFile, liftedOverFile)
		outFileName = FileNameGetter(fileName).get('_liftedOver.%s' % fileExt)
		outFh = CsvFileWriter(outFileName)
		fh = ReadFileAtOnceParser(fileName)
		outFh.write(fh.getSplittedLine())
		for splittedLine in fh:
			probeName, chrName, pos = splittedLine[:3]
			probePos = probeDict.get(probeName)
			if not probePos:
				print 'Passing probe %s' % probeName
				continue
			chrName, pos = probePos
			splittedLine[1] = chrName
			splittedLine[2] = pos
			outFh.write(splittedLine)
			
	def __getChrNameFromStr(self, chrName):
		if chrName[:3] != 'chr':
			chrName = 'chr' + chrName
		return chrName
	
	def _liftOverRawProbeFile(self, fileName, targetBuild):
		fh = ReadFileAtOnceParser(fileName)
		for splittedLine in fh:
			if splittedLine[0][0] != '#':
				header = splittedLine
				break
		outFileName = FileNameGetter(fileName).get('bed')
		outFh = CsvFileWriter(outFileName)
		probeIdx = header.index('ID')
		buildIdx = header.index('Genome Version')
		chrIdx = header.index('Chromosome')
		posIdx = header.index('Physical Position')
		buildName = None
		for splittedLine in fh:
			if not splittedLine[chrIdx]:
				break
			currentBuild = splittedLine[buildIdx]
			if not buildName:
				buildName = currentBuild
			elif buildName != currentBuild:
				raise NotImplementedError('Inconsistent build as "%s" and "%s" were found. Current line = [%s]' % (buildName, currentBuild, str(splittedLine)))
			probeName = splittedLine[0]
			chrName = splittedLine[chrIdx]
			pos = int(splittedLine[posIdx])
			outFh.write([self.__getChrNameFromStr(chrName), pos-1, pos, probeName])
		outFh.close()
		from liftOver import LiftOver
		return LiftOver(self.__binDir).process(outFileName, (buildName, targetBuild))
		
	def __getApproximatedPloidy(self, ploidy):
		if ploidy % 2 >= 1:
			ploidy = int(ploidy) + 1
		else:
			ploidy = int(ploidy)
		return ploidy
	
	def _getSampleToGroupDictFromFile(self, sampleFile, columnName, sampleAliasDict = None):
		header, fh, sampleIdx = self._getHeaderFhAndSampleIdxForSampleFile(sampleFile)
		sampleToGroupDict = {}
		groupIdx = header.index(columnName)
		sampleAliasColumnName = 'sampleAlias'
		sampleAliasIdx = None
		if sampleAliasColumnName in header:
			sampleAliasIdx = header.index(sampleAliasColumnName)
		for splittedLine in fh:
			sampleName = splittedLine[0]
			group = splittedLine[groupIdx].strip()
			if not group:
				continue
			if sampleAliasDict:
				if not sampleAliasDict.has_key(sampleName):
					print 'Passing sample %s because no alias' % sampleName
					continue
				sampleName = sampleAliasDict[sampleName]
			sampleToGroupDict[sampleName] = group
			if sampleAliasIdx is not None:
				sampleToGroupDict[splittedLine[sampleAliasIdx]] = group
		return sampleToGroupDict
	
	def _getHeaderFhAndSampleIdxForSampleFile(self, sampleFile):
		fh = ReadFileAtOnceParser(sampleFile)
		sampleIdx = header = None
		while fh.hasLinesLeft():
			splittedLine = fh.getSplittedLine()
			splittedLine2 = [columnName.lower() for columnName in splittedLine]
			if 'sample' in splittedLine2:
				header = splittedLine
				sampleIdx = splittedLine2.index('sample')
				break
		return header, fh, sampleIdx
	
	def _createPloidyFile(self, fileName, sampleAliasFile, sampleFile, columnName):
		sampleDict = None
		if sampleAliasFile:
			sampleDict = self._getSampleDictFromFile(sampleAliasFile)
		sampleToGroupDict = self._getSampleToGroupDictFromFile(sampleFile, columnName, sampleDict)
		#print sampleToGroupDict
		outFh = CsvFileWriter(FileNameGetter(fileName).get('_approx.txt'))
		fh = ReadFileAtOnceParser(fileName)
		header = fh.getSplittedLine()
		if header[1] == 'x':
			header[1] = 'ploidy'
		outFh.write(header + ['group'])
		for splittedLine in fh:
			sampleName = splittedLine[0]
			#if not sampleToGroupDict.has_key(sampleName):
				#print 'Passing sample %s' % sampleName
				#continue
			ploidy = float(splittedLine[1])
			outFh.write([splittedLine[0], self.__getApproximatedPloidy(ploidy), sampleToGroupDict.get(sampleName, '')])
		
	def _getSampleDictFromFile(self, sampleAliasFile):
		sampleDict = {}
		fh = ReadFileAtOnceParser(sampleAliasFile)
		for splittedLine in fh:
			sampleName1, sampleName2 = splittedLine
			sampleDict[sampleName1] = sampleName2
			sampleDict[sampleName2] = sampleName1
		return sampleDict
	
	def __getSampleDictFromClinicalFile(self, sampleFile):
		header, fh, sampleIdx = self._getHeaderFhAndSampleIdxForSampleFile(sampleFile)
		sampleDict = {}
		sampleAliasIdx = None
		if 'sampleAlias' in header:
			sampleAliasIdx = header.index('sampleAlias')
		for splittedLine in fh:
			sampleName = splittedLine[sampleIdx]
			if sampleAliasIdx is not None:
				sampleAlias = splittedLine[sampleAliasIdx]
				sampleDict[sampleName] = sampleAlias
				sampleDict[sampleAlias] = sampleName
			else:
				sampleDict[sampleName] = sampleName
		return sampleDict
	
	def __getTumorSampleListFromFile(self, platform, tumorSampleFile, sampleAliasFile = None):
		if sampleAliasFile:
			sampleDict = self._getSampleDictFromFile(sampleAliasFile)
		else:
			sampleDict = self.__getSampleDictFromClinicalFile(tumorSampleFile)
			#print 'sampleDict = ', sampleDict
		header, fh, sampleIdx = self._getHeaderFhAndSampleIdxForSampleFile(tumorSampleFile)
		tumorSampleList = []
		for splittedLine in fh:
			sampleName = splittedLine[sampleIdx]
			if not sampleDict.has_key(sampleName):
				print 'Passing sample %s' % sampleName
				continue
			sampleName = sampleDict[sampleName]
			tumorSampleList.append(sampleName)
		return tumorSampleList
	
	def __getTumorAndNormalLogRandBafIdxDictFromFile(self, lrrBafFile, tumorSampleList):
		fh = ReadFileAtOnceParser(lrrBafFile, 1)
		header = fh.getSplittedLine()
		idxDict = {'Normal': defaultdict(list), 'Tumor': defaultdict(list)}
		currentTumorSampleList = []
		normalSampleList = []
		for i, colName in enumerate(header):
			if i <= 2:
				continue
			sampleName = colName.split('.')[0]
			#print i, sampleName
			if sampleName in tumorSampleList:
				sampleType = 'Tumor'
				if sampleName not in currentTumorSampleList:
					currentTumorSampleList.append(sampleName)
			else:
				sampleType = 'Normal'
				if sampleName not in normalSampleList:
					normalSampleList.append(sampleName)
			dataType = None
			if '.Log R Ratio' in colName:
				dataType = 'LogR'
			elif '.B Allele Freq' in colName:
				dataType = 'BAF'
			if not dataType:
				continue
			idxDict[sampleType][dataType].append(i)
		return idxDict, currentTumorSampleList, normalSampleList
	
	def __createRscript(self, sampleFile, lrrBafFile, sampleAliasFile, gcFile, platform):
		tumorSampleList = self.__getTumorSampleListFromFile(platform, sampleFile, sampleAliasFile)
		print '%d tumor samples' % len(tumorSampleList)
		#print tumorSampleList
		idxDict, tumorSampleList, normalSampleList = self.__getTumorAndNormalLogRandBafIdxDictFromFile(lrrBafFile, tumorSampleList)
		#print idxDict
		snpPosFile = os.path.join(os.path.dirname(lrrBafFile), 'SNPpos_%s.txt' % os.path.basename(lrrBafFile).split('.')[0])
		cmd = 'cut -f 1-3 %s > %s' % (lrrBafFile, snpPosFile)
		Utilities()._runFunc(Utilities.mySystem, [cmd], snpPosFile)
		print len(tumorSampleList), len(normalSampleList), len(idxDict['Tumor']['LogR']), len(idxDict['Tumor']['BAF']), len(idxDict['Normal']['LogR']), len(idxDict['Normal']['BAF'])
		baseName = '.'.join(lrrBafFile.split('.')[:-1])
		rStr = '''X11.options(colortype="pseudo.cube")

baseName <- "%s"
lrrbaf = read.table("%s", header = T, sep = "\\t", row.names=1)
SNPpos = read.table("%s", header=T, sep="\\t", row.names=1)

normalSampleList <- %s
tumorSampleList <- %s

isMatchedNormalTumor <- length(normalSampleList) == length(tumorSampleList)
print(paste(c("isMatchedNormalTumor=", isMatchedNormalTumor)))

tumorLogRidxList <- %s
tumorBAFidxList <- %s
normalLogRidxList <- %s
normalBAFidxList <- %s

Tumor_LogR = lrrbaf[rownames(SNPpos),tumorLogRidxList,drop=F]
colnames(Tumor_LogR) = tumorSampleList

Tumor_BAF = lrrbaf[rownames(SNPpos),tumorBAFidxList,drop=F]
colnames(Tumor_BAF) = tumorSampleList

if (isMatchedNormalTumor){
	Normal_LogR = lrrbaf[rownames(SNPpos),normalLogRidxList,drop=F]
	colnames(Normal_LogR) = normalSampleList

	Normal_BAF = lrrbaf[rownames(SNPpos),normalBAFidxList,drop=F]
	colnames(Normal_BAF) = normalSampleList
}

#replace 2's by NA
Tumor_BAF[Tumor_BAF==2]=NA

# Tumor_LogR: correct difference between copy number only probes and other probes
CNprobes = substring(rownames(SNPpos),1,2)=="CN"

Tumor_LogR[CNprobes,1] = Tumor_LogR[CNprobes,1]-mean(Tumor_LogR[CNprobes,1],na.rm=T)
Tumor_LogR[!CNprobes,1] = Tumor_LogR[!CNprobes,1]-mean(Tumor_LogR[!CNprobes,1],na.rm=T)

if (isMatchedNormalTumor) {
	Normal_BAF[Normal_BAF==2]=NA
	Normal_LogR[CNprobes,1] = Normal_LogR[CNprobes,1]-mean(Normal_LogR[CNprobes,1],na.rm=T)
	Normal_LogR[!CNprobes,1] = Normal_LogR[!CNprobes,1]-mean(Normal_LogR[!CNprobes,1],na.rm=T)
}

# limit the number of digits:
Tumor_LogR = round(Tumor_LogR,4)

file.tumor.LogR <- paste(baseName, ".tumor.LogR.txt", sep="")
file.tumor.BAF <- paste(baseName, ".tumor.BAF.txt", sep="")
write.table(cbind(SNPpos,Tumor_BAF),file.tumor.BAF,sep="\t",row.names=T,col.names=NA,quote=F)
write.table(cbind(SNPpos,Tumor_LogR),file.tumor.LogR,sep="\t",row.names=T,col.names=NA,quote=F)

if (isMatchedNormalTumor){
	Normal_LogR = round(Normal_LogR,4)
	file.normal.LogR <- paste(baseName, ".normal.LogR.txt", sep="")
	file.normal.BAF <- paste(baseName, ".normal.BAF.txt", sep="")
	write.table(cbind(SNPpos,Normal_BAF),file.normal.BAF,sep="\t",row.names=T,col.names=NA,quote=F)
	write.table(cbind(SNPpos,Normal_LogR),file.normal.LogR,sep="\t",row.names=T,col.names=NA,quote=F)
}

#run ASCAT functions

library(ASCAT)
if (isMatchedNormalTumor){
	ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, file.normal.LogR, file.normal.BAF, chrs=1:22)
} else {
	ascat.bc <- ascat.loadData(file.tumor.LogR, file.tumor.BAF, chrs=1:22)
}

#GC correction
ascat.bc <- ascat.GCcorrect(ascat.bc, "%s")
ascat.plotRawData(ascat.bc)

if (isMatchedNormalTumor){
	ascat.bc <- ascat.aspcf(ascat.bc)
} else {
	ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, "%s")
	ascat.bc = ascat.aspcf(ascat.bc,ascat.gg=ascat.gg) 
}

ascat.plotSegmentedData(ascat.bc)
save(ascat.bc, baseName, file = "%s")

ascat.output <- ascat.runAscat(ascat.bc)
save(ascat.output, file = "%s")

#save ASCAT results

write.table(ascat.output$segments, file=paste(baseName,".segments.txt",sep=""), sep="\\t", quote=F, row.names=F)
write.table(ascat.output$aberrantcellfraction, file=paste(baseName,".acf.txt",sep=""), sep="\\t", quote=F, row.names=F)
write.table(ascat.output$ploidy, file=paste(baseName,".ploidy.txt",sep=""), sep="\\t", quote=F, row.names=F)

save.image(paste(baseName,".RData",sep=""))
''' % (baseName, lrrBafFile, snpPosFile, R()._getStrFromList(normalSampleList), R()._getStrFromList(tumorSampleList), R()._getStrFromList(idxDict['Tumor']['LogR']),
	   R()._getStrFromList(idxDict['Tumor']['BAF']), R()._getStrFromList(idxDict['Normal']['LogR']), R()._getStrFromList(idxDict['Normal']['BAF']), gcFile, platform, FileNameGetter(lrrBafFile).get('Ro'), FileNameGetter(lrrBafFile).get('_ASCAT.Ro'))
		rFileName = FileNameGetter(lrrBafFile).get('R')
		outFh = CsvFileWriter(rFileName)
		outFh.write(rStr)
		return rFileName, baseName + '.segments.txt'
	
	def __createCelFileListFile(self, celDirName):
		celFilePattern = os.path.join(celDirName, '*.cel')
		celFileList = glob.glob(os.path.join(celDirName, '*.cel'))
		if not celFileList:
			celFilePattern2 = os.path.join(celDirName, '*.cel.gz')
			celFileList = glob.glob(celFilePattern2)
			if celFileList:
				cmd = 'gunzip %s' % celFilePattern2
				Utilities.mySystem(cmd)
				celFileList = glob.glob(celFilePattern)
		if not celFileList:
			raise NotImplementedError('Could not find CEL files in dir %s' % celDirName)
		targetFileName = os.path.join(celDirName, 'list.nsp')
		cmd = 'echo "cel_files" > %s' % targetFileName
		Utilities.mySystem(cmd)
		cmd = 'ls %s >> %s' % (celFilePattern, targetFileName)
		Utilities.mySystem(cmd)
		return targetFileName
	
	def __getAptBinFolder(self):
		dirList = glob.glob(os.path.join(self.__binDir, 'apt-*', 'bin'))
		if len(dirList) != 1:
			raise NotImplementedError('Could not find APT binary folder in %s' % self.__binDir)
		return dirList[0]
	
	def __getCdfAndChrXfileFromDir(self, libDir):
		fileList = glob.glob(os.path.join(libDir, '*.cdf'))
		if len(fileList) != 1:
			raise NotImplementedError('Expecting one cdf file in dir %s but found %d:\n%s' % (libDir, len(fileList), '\n'.join(fileList)))
		cdfFile = fileList[0]
		fileList = glob.glob(os.path.join(libDir, '*.chrx'))
		if len(fileList) != 1:
			raise NotImplementedError('Expecting one chrx file in dir %s but found %d:\n%s' % (libDir, len(fileList), '\n'.join(fileList)))
		chrXfile = fileList[0]
		return cdfFile, chrXfile
	
	def __getGw6LibDirFromPlatform(self, gw6Dir, platform):
		if platform in ["Affy250k_sty", "Affy250k_nsp", "Affy500k"]:
			libDir = 'lib500k'
		elif platform == 'AffySNP6':
			libDir = 'lib'
		else:
			raise NotImplementedError('Unsupported plaform "%s"' % platform)
		libDir = os.path.join(gw6Dir, libDir)
		if not os.path.isdir(libDir):
			raise NotImplementedError('libDir %s does not exist...' % libDir)
		return libDir
	
	def __getFilteredFileList(self, fileList, platform):
		if len(fileList) > 1:
			if platform == 'Affy250k_sty':
				fileList = [fileName for fileName in fileList if '.sty.' in fileName]
			elif platform == 'Affy250k_nsp':
				fileList = [fileName for fileName in fileList if '.nsp.' in fileName]
			else:
				raise NotImplementedError('Unhandled platform "%s"' % platform)
		if len(fileList) != 1:
			raise NotImplementedError('Could not find targetSketchFile in %s' % gw6LibDir)
		return fileList
	
	def __getTargetSketchFileFromLibDirAndPlatform(self, gw6LibDir, platform):
		fileList = glob.glob(os.path.join(gw6LibDir, '*.normalization-target.txt'))
		fileList = self.__getFilteredFileList(fileList, platform)
		return fileList[0]
	
	def __getGenoClusterPfbFileTargetSketchFileAndBinDirFromDir(self, gw6Dir, platform):
		if platform not in ["Affy250k_sty", "Affy250k_nsp", "AffySNP6"]:
			raise NotImplementedError('Unsupported platform "%s"' % platform)
		gw6LibDir = self.__getGw6LibDirFromPlatform(gw6Dir, platform)
		fileList = glob.glob(os.path.join(gw6LibDir, '*.genocluster'))
		fileList = self.__getFilteredFileList(fileList, platform)
		genoClusterFile = fileList[0]
		fileList = glob.glob(os.path.join(gw6LibDir, '*.pfb'))
		if len(fileList) > 1:
			fileList = [fileName for fileName in fileList if 'liftedOver' in fileName]
		if len(fileList) != 1:
			raise NotImplementedError('Expecting one cdf file in dir %s but found %d:\n%s' % (gw6LibDir, len(fileList), '\n'.join(fileList)))
		pfbFile = fileList[0]
		binDir = os.path.join(gw6Dir, 'bin')
		if not os.path.isdir(binDir):
			raise NotImplementedError('gw6 bin dir %s not found' % binDir)
		targetSketchFileName = self.__getTargetSketchFileFromLibDirAndPlatform(gw6LibDir, platform)
		return genoClusterFile, pfbFile, targetSketchFileName, binDir
	
	def _runPennCnvAndGetLrrBafFile(self, celDirName, libDir, gw6Dir, platform):
		if not os.path.isdir(self.__binDir):
			raise NotImplementedError('Please specify bin dir where Affymetrix Power Tools (APT) is installed. The bin folder should contain a folder named "apt-*" which represents the APT folder. This APT folder should contain a bin folder containing all the necessary binaries.')
		cdfFile, chrXfile = self.__getCdfAndChrXfileFromDir(libDir)
		aptBinDir = self.__getAptBinFolder()
		celFileListFile = self.__createCelFileListFile(celDirName)
		targetDir = os.path.join(celDirName, 'apt_out')
		genoClusterFile, pfbFile, targetSketchFileName, gw6BinDir = self.__getGenoClusterPfbFileTargetSketchFileAndBinDirFromDir(gw6Dir, platform)
		
		print 'Step 1.1: Extracting genotypes from CEL files'
		cmd = '%s/apt-probeset-genotype -c %s --chrX-snps %s --out-dir %s --cel-files %s' % (aptBinDir, cdfFile, chrXfile, targetDir, celFileListFile)
		Utilities()._runFunc(Utilities.mySystem, [cmd], os.path.join(celDirName, 'step1.1'))
		
		print 'Step 1.2: Allele-specific signal extraction from CEL files'
		cmd = '%s/apt-probeset-summarize --cdf-file %s --out-dir %s --cel-files %s -a quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch %s' % (aptBinDir, cdfFile, targetDir, celFileListFile, targetSketchFileName)
		Utilities()._runFunc(Utilities.mySystem, [cmd], os.path.join(celDirName, 'step1.2'))
		
		print 'Step 1.4: LRR and BAF calculation'
		quantNormFile = os.path.join(targetDir, 'quant-norm.pm-only.med-polish.expr.summary.txt')
		if not os.path.isfile(quantNormFile):
			raise NotImplementedError('quantNormFile %s not found' % quantNormFile)
		lrrBafFile = os.path.join(targetDir, 'lrr_baf.txt')
		cmd = '%s/normalize_affy_geno_cluster.pl %s %s -locfile %s -out %s' % (gw6BinDir, genoClusterFile, quantNormFile, pfbFile, lrrBafFile)
		Utilities()._runFunc(Utilities.mySystem, [cmd], os.path.join(celDirName, 'step1.4'))
		return lrrBafFile
	
	def __checkAndInstallRpackagesIfNecessary(self, rDir):
		#RColorBrewer
		if R(rDir).isPackageInstalled('ASCAT'):
			return
		fileName = 'ASCAT.tar.gz'
		for cmd in ['wget https://github.com/Crick-CancerGenomics/ascat/archive/master.zip',
		            'unzip master.zip && rm master.zip',
		            'cd ascat-master && tar czf ../%s ASCAT && cd .. && rm -rf ascat-master' % fileName]:
			Utilities.mySystem(cmd)
		R(rDir).installPackage(fileName)
		Utilities.mySystem('rm %s' % fileName)
		
	def process(self, lrrBafFile, sampleFile, sampleAliasFile, gcFile, platform, libDir = None, gw6Dir = None):
		if os.path.isdir(lrrBafFile):
			lrrBafFile = Utilities.getFunctionResultWithCache(os.path.join(lrrBafFile, 'pennCNV.pyDump'), self._runPennCnvAndGetLrrBafFile, lrrBafFile, libDir, gw6Dir, platform)
			print 'lrrBafFile = %s' % lrrBafFile
		rFileName, ascatFile = self.__createRscript(sampleFile, lrrBafFile, sampleAliasFile, gcFile, platform)
		rDir = self.__binDir
		if self.__binDir and not os.path.isfile(os.path.join(self.__binDir, 'R')):
			rDir = None
		self.__checkAndInstallRpackagesIfNecessary(rDir)
		R(rDir).runScript(rFileName)
		return ascatFile


class CNView:
	def __init__(self, windowSize, percent, binDir = None, useShape = False, sampleFile = None, sampleAliasFile = None, groupColumnName = None):
		self.__windowSize = windowSize
		self.__percent = percent
		self.__binDir = binDir
		self.__binStr = ''
		if binDir:
			self.__binStr = binDir + '/'
		self.__useShape = useShape
		self.__sampleFile = sampleFile
		self.__sampleAliasFile = sampleAliasFile
		self.__groupColumnName = groupColumnName
		self.__sampleDict = None
		if sampleAliasFile:
			self.__sampleDict = RunAscat()._getSampleDictFromFile(sampleAliasFile)
		if sampleFile:
			self.__sampleToGroupDict = RunAscat()._getSampleToGroupDictFromFile(sampleFile, groupColumnName, self.__sampleDict)
	
	def processAll(self, ascatFile, chrFile, targetDir, ploidyFile, percentList, baseList, histogram = True, merge = False, dendrogram = False, plotAll = False, centromereFile = None, groupColumnName = None):
		Utilities.mySystem('mkdir -p %s' % targetDir)
		cluster = guessHpc()
		for baseNb in baseList:
			if not baseNb:
				continue
			cmdList = []
			cmd = 'python %s -f %s --ploidyFile %s -w %d -c %s --dendrogram %d --histogram %d -t %s -u %d -m %d --plotAll %d'
			if centromereFile:
				cmd += ' --centromereFile %s' % centromereFile
			if groupColumnName:
				cmd += ' -G "%s"' % groupColumnName
			if self.__binDir:
				cmd += ' -b %s' % self.__binDir
			cmdList.append((ascatFile, os.path.join(targetDir, os.path.basename(ascatFile) + '_%db' % baseNb), Cmd(cmd % (os.path.abspath(__file__), ascatFile, ploidyFile, baseNb, chrFile, dendrogram, histogram, targetDir, self.__useShape, merge, plotAll), nbCpus = 1, memory = 40)))
			ProcessFileFromCluster()._runCmdList(cmdList, ascatFile, cluster)
		for percent in percentList:
			cmdList = []
			cmd = 'python %s -f %s --ploidyFile %s -p %f -c %s --dendrogram %d --histogram %d -t %s -u %d -m %d --plotAll %d'
			if centromereFile:
				cmd += ' --centromereFile %s' % centromereFile
			if groupColumnName:
				cmd += ' -G "%s"' % groupColumnName
			if self.__binDir:
				cmd += ' -b %s' % self.__binDir
			cmdList.append((ascatFile, os.path.join(targetDir, os.path.basename(ascatFile) + '_%fp' % percent), Cmd(cmd % (os.path.abspath(__file__), ascatFile, ploidyFile, percent, chrFile, dendrogram, histogram, targetDir, self.__useShape, merge, plotAll), nbCpus = 1, memory = 40)))
			ProcessFileFromCluster()._runCmdList(cmdList, ascatFile, cluster)
		
	def __getChrSizeDictFromFile(self, fileName):
		chrSizeDict = {}
		for splittedLine in ReadFileAtOnceParser(fileName):
			chrName = splittedLine[0]
			if chrName[:3] == 'chr':
				chrName = chrName[3:]
			chrSizeDict[chrName] = int(splittedLine[1])
		return chrSizeDict
	
	def __getSampleListFromFile(self, fileName):
		sampleList = []
		fh = ReadFileAtOnceParser(fileName)
		fh.getSplittedLine()
		for splittedLine in fh:
			sampleList.append(splittedLine[0])
		return sampleList
	
	def __getNextSampleLineDictAndSampleName(self, fh):
		lineDict = SimpleMultiDict()
		sampleName = None
		while fh.hasLinesLeft():
			splittedLine = fh.getSplittedLine()
			currentSampleName = splittedLine[0]
			if not sampleName:
				sampleName = currentSampleName
			if currentSampleName != sampleName:
				fh.restore(splittedLine)
				break
			lineDict[splittedLine[1]] = splittedLine[2:]
		return lineDict, sampleName
	
	def _getLineListInBetween(self, lineList, start, endSegment):
		fragPos = Position(1, start, endSegment)
		selectedLineList = []
		while lineList:
			splittedLine = lineList.pop(0)
			currentStart = int(splittedLine[0])
			currentEnd = int(splittedLine[1])
			currentPos = Position(1, currentStart, currentEnd)
			overlapPos = fragPos.getOverlapPosition(currentPos)
			if overlapPos:
				selectedLineList.append([overlapPos.start, overlapPos.end] + splittedLine[2:])
				if currentEnd > fragPos.end:
					lineList.insert(0, [fragPos.end+1, currentEnd] + splittedLine[2:])
					break
			elif currentStart > fragPos.end:
				lineList.insert(0, splittedLine)
				break
		return selectedLineList
	
	def __getCurrentPloidyFromLineList(self, lineList):
		totalSize = totalPloidy = 0.
		for splittedLine in lineList:
			start = int(splittedLine[2])
			end = int(splittedLine[3])
			fragSize = end - start + 1
			totalSize += fragSize
			ploidy = int(splittedLine[4]) + int(splittedLine[5])
			totalPloidy += ploidy * fragSize
		return int(round(totalPloidy / totalSize))
	
	def __getCurrentPloidyFromLineList2(self, lineList):
		fragDict = defaultdict(int)
		for splittedLine in lineList:
			start = int(splittedLine[2])
			end = int(splittedLine[3])
			fragSize = end - start + 1
			ploidy = int(splittedLine[4]) + int(splittedLine[5])
			fragDict[ploidy] += fragSize
			#valueList.append(ploidy)
		#valueList.sort()
		#return getMedianFromList(valueList)
		fragSizeToPloidyDict = SimpleMultiDict([[fragSize, ploidy] for ploidy, fragSize in fragDict.iteritems()])
		maxFragSize = max(fragSizeToPloidyDict)
		ploidyList = fragSizeToPloidyDict.getall(maxFragSize)
		if len(ploidyList) > 1:
			print 'WARNING: several fragments with same size: %d -> %s' % (maxFragSize, str(ploidyList))
		return ploidyList[0]
		
	def __getNewSegment(self, prevSegment, currentSegment, defaultPloidy):
		prevEnd = int(prevSegment[1])
		currentStart = int(currentSegment[0])
		newStart = prevEnd + 1
		newEnd = currentStart - 1
		if newEnd >= newStart:
			return [newStart, newEnd, defaultPloidy, 0]
	
	def _getGappedFilledSegmentListFromList(self, sampleName, chrName, segmentList, startSegment, endSegment, defaultPloidy = 2):
		#print '+' * 100
		start = int(segmentList[0][0])
		end = int(segmentList[-1][1])
		newSegmentList = []
		if startSegment < start:
			newSegmentList.append([sampleName, chrName, startSegment, start-1, defaultPloidy, 0])
		segment = segmentList.pop(0)
		newSegmentList.append([sampleName, chrName] + segment)
		for currentSegment in segmentList:
			newSegment = self.__getNewSegment(segment, currentSegment, defaultPloidy)
			#print '@', currentSegment, newSegment
			if newSegment:
				newSegmentList.append([sampleName, chrName] + newSegment)
			newSegmentList.append([sampleName, chrName] + currentSegment)
			segment = currentSegment
		if end < endSegment:
			newSegmentList.append([sampleName, chrName, end+1, endSegment, defaultPloidy, 0])
		return newSegmentList
		
	def _getSegmentListFromList(self, lineList, sampleName, chrName, chrSizeDict, defaultPloidy = 2):
		middleSegmentList = []
		start = int(lineList[0][0])
		chrEnd = chrSizeDict[str(chrName)]
		end = min(int(lineList[-1][1]), chrEnd)
		startSegmentList = self._getSegmentStartList(start, sampleName, chrName, defaultPloidy)
		
		endSegmentList = self._getEndSegmentList(end, sampleName, chrName, chrEnd, defaultPloidy)
		startSegmentNb = self.__getNbSegmentsFromPos(start)
		startSegment = startSegmentNb * self.__windowSize + 1
		endSegment = startSegment + self.__windowSize - 1
		while lineList:
			totalSize = 0
			currentLineList = self._getLineListInBetween(lineList, startSegment, endSegment)
			#print 'C = %s, S = %d, E = %d, L = %d' % (chrName, startSegment, endSegment, len(currentLineList))
			#if currentLineList:
				#print currentLineList
			if currentLineList:
				currentLineList = self._getGappedFilledSegmentListFromList(sampleName, chrName, currentLineList, startSegment, endSegment, defaultPloidy)
				middleSegmentList.append([sampleName, chrName, startSegment, endSegment, self.__getCurrentPloidyFromLineList2(currentLineList)])
			else:
				middleSegmentList.append([sampleName, chrName, startSegment, endSegment, defaultPloidy])
			startSegment += self.__windowSize
			endSegment = min(endSegment + self.__windowSize, chrEnd)
			if startSegment > endSegment:
				break
		return startSegmentList + middleSegmentList + endSegmentList
	
	def __getNbSegmentsFromPos(self, pos):
		residue = pos % self.__windowSize
		nbSegments = (pos - residue) / self.__windowSize
		if not residue:
			nbSegments -= 1
		return nbSegments
	
	def _getSegmentStartList(self, start, sampleName, chrName, defaultPloidy = 2):
		nbSegments = self.__getNbSegmentsFromPos(start)
		return [[sampleName, chrName, 1 + i * self.__windowSize, self.__windowSize * (i+1), defaultPloidy] for i in range(nbSegments)]
	
	def _getEndSegmentList(self, end, sampleName, chrName, chrEnd, defaultPloidy = 2):
		nbSegments = self.__getNbSegmentsFromPos(end)
		endSegmentNb = self.__getNbSegmentsFromPos(chrEnd)
		segmentList = [[sampleName, chrName, 1 + i * self.__windowSize, self.__windowSize * (i+1), defaultPloidy] for i in range(nbSegments+1, endSegmentNb)]
		lastSegmentStart = endSegmentNb * self.__windowSize + 1
		if lastSegmentStart > end:
			segmentList.append([sampleName, chrName, lastSegmentStart, chrEnd, defaultPloidy])
		return segmentList
	
	def __mergeLastTwoSegmentsIfNecessary(self, segmentList, i = -1, reverse = False, mergeAnyway = False):
		sampleName1, chrName1, lastSegmentStart1, chrEnd1, defaultPloidy1 = segmentList[i]
		sampleName0, chrName0, lastSegmentStart0, chrEnd0, defaultPloidy0 = segmentList[i-1]
		if sampleName1 != sampleName0:
			raise NotImplementedError('Expecting same samples but found %s and %s:%s\n%s' % (sampleName0, sampleName1, segmentList[i-1], segmentList[i]))
		if chrName0 != chrName1:
			raise NotImplementedError('Expecting same chr but found %s and %s: %s\n%s' % (chrName0, chrName1, segmentList[i-1], segmentList[i]))
		segmentSize1 = chrEnd1 - lastSegmentStart1 + 1
		segmentSize0 = chrEnd0 - lastSegmentStart0 + 1
		if not mergeAnyway and ((not reverse and segmentSize1 >= self.__windowSize * 1. / 2) or (reverse and segmentSize0 >= self.__windowSize * 1. / 2)):
			#print 'No merge', sampleName0, chrName0, lastSegmentStart0, chrEnd0, lastSegmentStart1, chrEnd1
			return
		segmentList.pop(i-1)
		newPloidy = int(round((1. * defaultPloidy0 * segmentSize0 + defaultPloidy1 * segmentSize1) / (segmentSize0 + segmentSize1)))
		lastIndex = i-1
		if i == -1:
			lastIndex = i
		segmentList[lastIndex] = [sampleName0, chrName0, lastSegmentStart0, chrEnd1, newPloidy]
		#print 'There', sampleName0, chrName0, lastSegmentStart0, chrEnd0, lastSegmentStart1, chrEnd1, defaultPloidy0, defaultPloidy1, newPloidy
	
	def __getNewLineStartAndEnd(self, linePos, overlapPos):
		if linePos.start < overlapPos.start:
			start = linePos.start
			end = overlapPos.start - 1
		else:
			start = overlapPos.end + 1
			end = linePos.end
		if start > end:
			print linePos, overlapPos, start, end
			raise NotImplementedError
		return start, end
	
	def _removeEventOverCentromereFromList(self, lineList, centroStart, centroEnd):
		#if (self.__sampleName, self.__chrName) == ('GSM248805', 9):
			#print '~' * 55
			#print '_removeEventOverCentromereFromList'
		chrPos = Position(None, centroStart, centroEnd)
		toRemoveList = []
		shift = 0
		if lineList and not ValueParser().isNb(lineList[0][0]):
			shift = 2
		#print '@' * 100
		toReplaceList = []
		for line in lineList:
			lineStart, lineEnd = line[shift:shift+2]
			linePos = Position(None, lineStart, lineEnd)
			overlapPos = linePos.getOverlapPosition(chrPos)
			#print 'linePos', linePos
			#print 'overrl', overlapPos
			if overlapPos:
				if overlapPos == linePos:
					#print '=='
					toRemoveList.append(line)
				else:
					start, end = self.__getNewLineStartAndEnd(linePos, overlapPos)
					#print 'not eql', start, end
					endLine = None
					if int(line[shift+1]) > centroEnd and centroStart > int(line[shift+0]):
						endLine = [centroEnd+1, line[shift+1]] + line[shift+2:]
						if shift:
							endLine = line[:2] + endLine
						#print 'ENDline', endLine, line[1], end
					line[shift+0] = start
					line[shift+1] = end
					if endLine:
						toReplaceList.append((line, endLine))
		#if (self.__sampleName, self.__chrName) == ('GSM248805', 9):
			#print '%d toRemove, %d toReplace' % (len(toRemoveList), len(toReplaceList))
			#print 'toReplace:'
		for line, endLine in toReplaceList:
			#if (self.__sampleName, self.__chrName) == ('GSM248805', 9):
				#print line, endLine
			idx = lineList.index(line)
			lineList.insert(idx+1, endLine)
		#if (self.__sampleName, self.__chrName) == ('GSM248805', 9):
			#print '@' * 20
			#print 'toRemove:'
		for line in toRemoveList:
			#if (self.__sampleName, self.__chrName) == ('GSM248805', 9):
				#print line
			lineList.remove(line)
	
	def _mergeCentromericSegmentsIfNecessary(self, segmentList, mergeCentromereSegments = None, nbExpectedSegments = None):
		prevSegment = None
		for i, segment in enumerate(segmentList):
			if prevSegment and not mergeCentromereSegments and prevSegment[3] != segment[2] - 1:
				#print 'Merge', i
				nbSegments = len(segmentList)
				self.__mergeLastTwoSegmentsIfNecessary(segmentList, i-1)
				#print 1, segmentList
				shift = nbSegments - len(segmentList)
				if shift:
					idx = i
				else:
					idx = i+1
				self.__mergeLastTwoSegmentsIfNecessary(segmentList, idx, True)
				#print 2, segmentList
				break
			elif prevSegment and mergeCentromereSegments and prevSegment[3] != segment[2] - 1:
				nbSegments = len(segmentList)
				if nbSegments == nbExpectedSegments:
					continue
				#print '~' * 50
				#print prevSegment
				#print segment
				#print ']' * 30
				self.__mergeLastTwoSegmentsIfNecessary(segmentList, i)
				#print 1, segmentList
				#shift = nbSegments - len(segmentList)
				#if shift:
					#idx = i+1
				#else:
					#idx = i+2
				#self.__mergeLastTwoSegmentsIfNecessary(segmentList, idx, True)
				break
			prevSegment = segment
	
	def _fillCentromericGapsIfNecessary(self, segmentList, defaultPloidy):
		prevSegment = None
		for i, segment in enumerate(segmentList):
			if prevSegment and prevSegment[3] != segment[2] - 1:
				#print 'Merge', i
				segmentList.insert(i, segment[:2]+[prevSegment[3]+1, segment[2] - 1, defaultPloidy])
				break
			prevSegment = segment
			
	def __getCentromericSegmentFromCentromerePosition(self, centroPos):
		start = 1
		while True:
			pos = Position(centroPos.ctgId, start, start + self.__windowSize - 1)
			if pos.getOverlapPosition(centroPos):
				return pos
			start = pos.end + 1
			
	def __mergeCentromericSegmentsIntoOne(self, segmentList, chrName, centroStart, centroEnd):
		#print '__mergeCentromericSegmentsIntoOne ' + '@' * 100
		centroSegmentDict = SimpleMultiDict()
		centroPos = Position(str(chrName), centroStart, centroEnd)
		centroSegment = self.__getCentromericSegmentFromCentromerePosition(centroPos)
		#print centroPos, centroSegment, len(segmentList)
		for i, (sampleName, chrName, startSegment, endSegment, ploidy) in enumerate(segmentList):
			pos = Position(str(chrName), startSegment, endSegment)
			#print 'pos', pos
			if pos.getOverlapPosition(centroSegment):
				#print 'POS', pos
				centroSegmentDict[sampleName] = i
		for sampleName in centroSegmentDict:
			currentSegmentIdxList = centroSegmentDict.getall(sampleName)
			if len(currentSegmentIdxList) > 2:
				raise NotImplementedError('#segments > 2 for sample = %s\nSegments = %s' % (sampleName, currentSegmentIdxList))
			idx = currentSegmentIdxList[0]
			for i in currentSegmentIdxList[1:]:
				if i - idx != 1:
					raise NotImplementedError('Indexes are not consecutive for sample %s, idx = %s' % (sampleName, currentSegmentIdxList))
				idx = i
			if len(currentSegmentIdxList) == 2:
				#print '@' * 100
				#print 'Sample %s: merging segmentList %s' % (sampleName, [segmentList[i] for i in currentSegmentIdxList])
				maxIdx = max(currentSegmentIdxList)
				self.__mergeLastTwoSegmentsIfNecessary(segmentList, maxIdx)
				#print segmentList[maxIdx-1:maxIdx+1]
				#print '=' * 50
		
	def __getChrSegmentListFromPosition(self, start, end, chrName):
		endSegmentNb = int(math.ceil(1. * end / self.__windowSize))
		startSegmentNb = int(1. * start / self.__windowSize)
		posList = []
		if endSegmentNb == startSegmentNb:
			endSegmentNb += 1
		for i in range(startSegmentNb, endSegmentNb):
			posList.append(Position(chrName, i * self.__windowSize + 1, (i+1) * self.__windowSize))
		return posList
	
	def __cutSegmentsAccordingToStep(self, segmentList, chrName):
		newSegmentList = []
		chrName = str(chrName)
		segmentDict = SimpleMultiDict()
		for segment in segmentList:
			segChrName = str(segment[1])
			if segChrName != chrName:
				raise NotImplementedError('Inconsistent chrs: %s != %s' % (segChrName, chrName))
			chrSegmentList = self.__getChrSegmentListFromPosition(segment[2], segment[3], chrName)
			pos = Position(segChrName, *segment[2:4])
			#print 'Current seg', segment, self.__windowSize
			#print 'chr pos:'
			#for chrPos in chrSegmentList:
				#print chrPos
			#print '@' * 10
			for chrPos in chrSegmentList:
				overlapPos = pos.getOverlapPosition(chrPos)
				currentSegment = copy.copy(segment)
				#print chrName, chrPos, segment, overlapPos
				currentSegment[2] = overlapPos.start
				currentSegment[3] = overlapPos.end
				#newSegmentList.append(currentSegment)
				segmentDict[(chrPos.start, chrPos.end)] = currentSegment
		for chrStart, chrEnd in Utilities.getOrderedKeys(segmentDict):
			#print '@' * 50
			#print chrStart, chrEnd
			currentSegmentList = segmentDict.getall((chrStart, chrEnd))
			#print currentSegmentList
			if len(currentSegmentList) > 1:
				self.__mergeLastTwoSegmentsIfNecessary(currentSegmentList, mergeAnyway = True)
			if len(currentSegmentList) != 1:
				print currentSegmentList
				raise NotImplementedError
			currentSegment = currentSegmentList[0]
			currentSegment[2] = chrStart
			currentSegment[3] = chrEnd
			newSegmentList.append(currentSegment)
		segmentList[:] = newSegmentList
	
	def __createMatrixFileAndGetDicts(self, ascatFile, chrFile, targetDir, ploidyFile, centromereFile = None, mergeCentromereSegments = None):
		centromereDict = None
		if centromereFile:
			centromereDict = self.__getCentromereDictFromFile(centromereFile)
			print 'CENTRO:', centromereDict
			ascatFile = Utilities.getFunctionResultWithCache(FileNameGetter(ascatFile).get('_%s.pyDump' % os.path.basename(centromereFile).split('.')[0]), self.__createFileWithTruncatedCentromereData, ascatFile, centromereFile, centromereDict)
		if self.__windowSize:
			self.__suffix = str(self.__windowSize) + 'nt'
		elif self.__percent:
			self.__suffix = str(self.__percent) + 'pc'
		ploidyDict = DefaultPloidyDict()
		if ploidyFile:
			ploidyDict = self.__getPloidyDictFromFile(ploidyFile)
		#print ascatFile, len(ploidyDict), ploidyDict
		Utilities.mySystem('mkdir -p %s' % targetDir)
		chrSizeDict = self.__getChrSizeDictFromFile(chrFile)
		#sampleList = self.__getSampleListFromFile(ascatFile)
		#print str(len(sampleList))+" samples found : "+", ".join(sampleList)
		outFileName, matrixDict = Utilities.getFunctionResultWithCache(ascatFile + '_%s_%d.pyDump5' % (self.__suffix, mergeCentromereSegments), self.__createSegmentPloidyFile, ascatFile, chrSizeDict, targetDir, ploidyDict, centromereDict, mergeCentromereSegments)
		
		return outFileName, matrixDict, ploidyDict, chrSizeDict, centromereDict
	
	def __calculatePloidyFromLine(self, splittedLine):
		ploidyDict = defaultdict(int)
		for ploidy in splittedLine[1:]:
			ploidyDict[int(ploidy)] += 1
		maxNb = max(ploidyDict.values())
		#print ploidyDict, maxNb
		maxPloidyList = []
		for ploidy, nb in ploidyDict.iteritems():
			if nb == maxNb:
				maxPloidyList.append(ploidy)
		if len(maxPloidyList) > 1:
			raise NotImplementedError('Could not deduce ploidy for sample %s as there are as many segments with ploidies %s' % (splittedLine[0], maxPloidyList))
		return maxPloidyList[0]
	
	def _createPloidyFile2(self, ascatFile, ploidyFile = None):
		fh = ReadFileAtOnceParser(ascatFile)
		header = fh.getSplittedLine()
		outFh = CsvFileWriter(FileNameGetter(ascatFile).get('_ploidy2.txt'))
		outFh.write(['Sample', 'Ploidy'])
		ploidyDict = defaultdict(int)
		baseDict = defaultdict(int)
		for splittedLine in fh:
			sampleName, chrName, start, end, nMajor, nMinor = splittedLine
			nbBases = int(end) - int(start) + 1
			baseDict[sampleName] += nbBases
			ploidyDict[sampleName] += (int(nMajor) + int(nMajor)) * nbBases
		if ploidyFile:
			otherPloidyDict = self.__getPloidyDictFromFile(ploidyFile)
		nbErrors = 0
		for sampleName in Utilities.getOrderedKeys(ploidyDict):
			ploidy = ploidyDict[sampleName]
			nbBases = baseDict[sampleName]
			ploidy = int(1. * ploidy / nbBases)
			if ploidyFile and otherPloidyDict.has_key(sampleName):
				otherPloidy = otherPloidyDict[sampleName]
				if ploidy != otherPloidy:
					print 'Different ploidies for sample %s: %d != %d' % (sampleName, ploidy, otherPloidy)
					print nbBases, ploidyDict[sampleName]
					nbErrors += 1
			outFh.write([sampleName, ploidy])
		if ploidyFile:
			print 'nbErrors = %d' % nbErrors
		
	def _createPloidyFile(self, ascatFile, chrFile, targetDir, ploidyFile, centromereFile = None, mergeCentromereSegments = None, otherPloidyFile = None):
		outFileName, matrixDict, ploidyDict, chrSizeDict, centromereDict = self.__createMatrixFileAndGetDicts(ascatFile, chrFile, targetDir, ploidyFile, centromereFile, mergeCentromereSegments)
		matrixFile = self.__createMatrixFile(matrixDict, targetDir, None, ploidyFile)
		print 'Matrix file: ', matrixFile
		fh = ReadFileAtOnceParser(matrixFile)
		targetFileName = FileNameGetter(outFileName).get('_ploidy.txt')
		outFh = CsvFileWriter(targetFileName)
		header = fh.getSplittedLine()
		nbErrors = 0
		header += ['ploidy']
		if self.__sampleFile:
			header += ['group']
		outFh.write(header)
		if otherPloidyFile:
			otherPloidyDict = self.__getPloidyDictFromFile(otherPloidyFile)
		for splittedLine in fh:
			ploidy = self.__calculatePloidyFromLine(splittedLine)
			sampleName = splittedLine[0]
			if otherPloidyFile and otherPloidyDict.has_key(sampleName):
				otherPloidy = otherPloidyDict[sampleName]
				if ploidy != otherPloidy:
					print 'Different ploidies for sample %s: %d != %d' % (sampleName, ploidy, otherPloidy)
					nbErrors += 1
				#else:
					#print 'Ok for sample %s' % splittedLine[0]
			splittedLine += [ploidy]
			if self.__sampleFile:
				splittedLine += [self.__sampleToGroupDict.get(sampleName, '')]
			outFh.write(splittedLine)
		if otherPloidyFile:
			print 'nbErrors = %d' % nbErrors
		return targetFileName
	
	def __updateSegmentLengths(self, segmentList, chrName):
		segmentList[:] = self.__cutSegmentsAccordingToStep(segmentList, chrName)
		#segmentDict = 
				
	def __getSegmentLineListFromLineDict(self, lineDict, sampleName, chrSizeDict, defaultPloidy = 2, centromereDict = None, isLoh = False, mergeCentromereSegments = False):
		segmentList = []
		nbExpectedSegments = None
		if self.__percent:
			nbExpectedSegments = int(100 / self.__percent)
		for chrName in range(1, 23):
			chrLineList = lineDict.getall(str(chrName))
			if self.__percent:
				self.__windowSize = int(math.ceil(chrSizeDict[str(chrName)] * self.__percent / 100.))
				print 'Using %d window size for chr %s' % (self.__windowSize, chrName)
			#print 'Chr %s' % chrName
			#print len(chrLineList)
			if not chrLineList:
				continue
			#if isLoh and centromereDict:
				#if 'P51T' in sampleName:
					#print '@' * 100
					#print chrLineList
				#self._removeEventOverCentromereFromList(chrLineList, *centromereDict[str(chrName)])
				#if 'P51T' in sampleName:
					#print ':' * 50
					#print chrLineList
			
			currentSegmentList = self._getSegmentListFromList(chrLineList, sampleName, chrName, chrSizeDict, defaultPloidy)
			#if sampleName == 'GSM248805' and chrName == 9:
				#print isLoh, ':' * 50
				#print centromereDict
				#print chrLineList
				#print currentSegmentList
			if not isLoh and centromereDict:
				#if 'P51T' in sampleName:
					#print '>' * 50
					#print currentSegmentList
				self.__sampleName = sampleName
				self.__chrName = chrName
				self._removeEventOverCentromereFromList(currentSegmentList, *centromereDict[str(chrName)])
				#if sampleName == 'GSM248805' and chrName == 9:
					#print ':' * 50
					#print currentSegmentList
			#if 'P51T' in sampleName:
				#print '=' * 50
				#print currentSegmentList
				if mergeCentromereSegments:
					#print 'MERGING'
					#print currentSegmentList
					self.__cutSegmentsAccordingToStep(currentSegmentList, chrName)
					#print currentSegmentList
					#print '~' * 20
				else:
					self._mergeCentromericSegmentsIfNecessary(currentSegmentList)
			##if not isLoh and centromereDict:
				##self._fillCentromericGapsIfNecessary(currentSegmentList, defaultPloidy)
			#if 'P51T' in sampleName:
				#print '>' * 50
				#print currentSegmentList
			if isLoh and centromereDict:
				self._removeEventOverCentromereFromList(currentSegmentList, *centromereDict[str(chrName)])
				self._mergeCentromericSegmentsIfNecessary(currentSegmentList)
				##self._fillCentromericGapsIfNecessary(currentSegmentList, 0)
				
			#if 'P51T' in sampleName:
				#print '<' * 50
				#print currentSegmentList
			self.__mergeLastTwoSegmentsIfNecessary(currentSegmentList)
			#if sampleName == 'GSM248805' and chrName == 3:
				#print '<' * 50
				#print currentSegmentList
			segmentList += currentSegmentList
			#if mergeCentromereSegments:
				#self.__mergeCentromericSegmentsIntoOne(segmentList, chrName, *centromereDict[str(chrName)])
		return segmentList
	
	def __createSegmentPloidyFile(self, ascatFile, chrSizeDict, targetDir, ploidyDict, centromereDict = None, mergeCentromereSegments = False):
		fh = ReadFileAtOnceParser(ascatFile)
		header = fh.getSplittedLine()
		outFileName = os.path.join(targetDir, os.path.basename(ascatFile).split('.')[0] + '_segments_%s.txt' % self.__suffix)
		outFh = CsvFileWriter(outFileName)
		matrixDict = {}
		isLoh = sum(ploidyDict.values()) == 0
		while fh.hasLinesLeft():
			lineDict, sampleName = self.__getNextSampleLineDictAndSampleName(fh)
			defaultPloidy = ploidyDict.get(sampleName.split('.')[0], 2)
			if not ploidyDict.has_key(sampleName):
				print 'Passing sample %s' % sampleName
				continue
			print 'Processing sample %s with ploidy %d' % (sampleName, defaultPloidy)
			segmentList = self.__getSegmentLineListFromLineDict(lineDict, sampleName, chrSizeDict, defaultPloidy, centromereDict, isLoh, mergeCentromereSegments)
			outFh.writeAllLinesAtOnce(segmentList)
			matrixDict[sampleName.split('.')[0]] = [['%s:%d-%d' % (chrName, start, end), alleleNb] for currentSampleName, chrName, start, end, alleleNb in segmentList]
		return outFileName, matrixDict
	
	def __getSegmentListFromFile(self, fileName):
		fh = ReadFileAtOnceParser(fileName)
		segmentList = []
		for splittedLine in fh:
			segmentList.append((splittedLine[0], splittedLine[-2:]))
		return segmentList
	
	def _compareSegmentFiles(self, fileName, fileName2):
		fh1 = ReadFileAtOnceParser(fileName)
		fh2 = ReadFileAtOnceParser(fileName2)
		for splittedLine1 in fh1:
			splittedLine2 = fh2.getSplittedLine()
			if splittedLine1[-2:] != splittedLine2[-2:]:
				if abs(int(splittedLine1[-1]) - int(splittedLine2[-1])) > 1:
					print 'difference', splittedLine1[:2], splittedLine2[:2], splittedLine1[-2:], splittedLine2[-2:]
	
	def __getPloidyDictFromFile(self, fileName):
		ploidyDict = {}
		fh = ReadFileAtOnceParser(fileName)
		header = fh.getSplittedLine()
		ploidyIdx = len(header) - 1
		if 'ploidy' in header:
			ploidyIdx = header.index('ploidy')
		header = [columnName.lower() for columnName in header]
		sampleIdx = header.index('sample')
		if not header[ploidyIdx].strip():
			ploidyIdx -= 1
		#print len(header), header
		for splittedLine in fh:
			#print splittedLine
			ploidyDict[splittedLine[sampleIdx].split('_')[0]] = int(splittedLine[ploidyIdx])
		return ploidyDict
	
	def __getPloidyKeyFromPloidy(self, currentPloidy):
		return max(-4, currentPloidy) if currentPloidy < 0 else min(6, currentPloidy)
	
	def __getLohDataDictFromList(self, dataList):
		dataDict = {}
		for sampleName, in dataList:
			return
		return dataDict
	
	def __getPosFromKey(self, key):
		chrName, key = key.split(':')
		start, end = key.split('-')
		return Position(chrName, start, end)
	
	def __getKeyListOverlappingPos(self, currentLohDataDict, chrPos):
		keyList = []
		for key in currentLohDataDict:
			pos = self.__getPosFromKey(key)
			if chrPos.getOverlapPosition(pos):
				keyList.append(key)
		return keyList
	
	def __getLohNbDictAndKeyListFromChrPos(self, lohDataDict, chrPos):
		lohNbDict = defaultdict(int)
		chrPos = self.__getPosFromKey(chrPos)
		keySet = set()
		for sampleName, currentLohDataDict in lohDataDict.iteritems():
			currentPloidy = currentLohDataDict.get(chrPos)
			for key in self.__getKeyListOverlappingPos(currentLohDataDict, chrPos):
				keySet.add(key)
				currentPloidy = currentLohDataDict.get(chrPos)
				if currentPloidy:
					lohNbDict[key] += 1
		keyList = []
		for key in keySet:
			pos = self.__getPosFromKey(key)
			keyList.append(((pos.ctgId, pos.start, pos.end), key))
		keyList.sort()
		return lohNbDict, [key for pos, key in keyList]
	
	def __getLohNbFromDictAndKey(self, lohDataDict, chrPos):
		nbLoh = 0
		for sampleName, currentLohDataDict in lohDataDict.iteritems():
			currentPloidy = currentLohDataDict.get(chrPos)
			if currentPloidy:
				nbLoh += 1
		return nbLoh
	
	def __writeLOH(self, lohDataDict, outFh2, nbVal, chrPos, chrName, start, end):
		nbLoh = self.__getLohNbFromDictAndKey(lohDataDict, chrPos)
		if not nbLoh:
			lohNbList = []
			lohNbDict, keyList = self.__getLohNbDictAndKeyListFromChrPos(lohDataDict, chrPos)
			for key in keyList:
				nbLoh = self.__getLohNbFromDictAndKey(lohDataDict, key)
				lohNbList.append(nbLoh)
		else:
			lohNbList = [nbLoh]
		#lohPointList.append((chrName, start, end, ))
		for nbLoh in lohNbList:
			percent = nbLoh * 100. / nbVal
			#print 'CHR', chrPos, chrName, start, end, nbLoh, percent
			outFh2.write([start, chrName, -percent])
			outFh2.write([end, chrName, -percent])
	
	def __createHistDataFileAndGetMaxValue(self, dataList, targetDir, ploidyDict, keyword, lohDataDict = None):
		fileName = os.path.join(targetDir, '%s_hist_%s.txt' % (keyword, self.__suffix))
		outFh = CsvFileWriter(fileName)
		maxValue = lohFileName = None
		if lohDataDict:
			lohFileName = os.path.join(targetDir, '%s_hist_%s_loh.txt' % (keyword, self.__suffix))
			outFh2 = CsvFileWriter(lohFileName)
		dataList0 = dataList[0][-1]
		lohPointList = []
		for i, (chrPos, currentPloidy) in enumerate(dataList0):
			ploidyList = []
			countDict = defaultdict(int)
			nbVal = nbLoh = 0
			for sampleName, currentDataList in dataList:
				currentPloidy = currentDataList[i][1]
				defaultPloidy = ploidyDict[sampleName.split('_')[0]]
				currentPloidy -= defaultPloidy
				ploidyKey = self.__getPloidyKeyFromPloidy(currentPloidy)
				nbVal += 1
				countDict[ploidyKey] += 1
			chrName, pos = chrPos.split(':')
			start, end = pos.split('-')
			start = int(start)
			end = int(end)
			if lohDataDict:
				self.__writeLOH(lohDataDict, outFh2, nbVal, chrPos, chrName, start, end)
			for key, val in countDict.iteritems():
				val = val * 100. / nbVal
				countDict[key] = val
			#countDict = {valueDict[currentPloidy]: 100}
			
			histStart = start+(end-start)/2
			length = end-start+1
			# Order of outputs bellow is VERY important for color coding in the histogram !
			currentValue1 = currentValue2 = 0
			for i in range(1, 7):
				value = countDict.get(i, 0)
				outFh.write([i, chrName, histStart, length, value])
				currentValue1 += value
			for i in range(-1, -5, -1):
				value = countDict.get(i, 0)
				outFh.write([i, chrName, histStart, length, -value])
				currentValue2 += value
			maxValue = max(maxValue, currentValue1, currentValue2)
		return fileName, maxValue, lohFileName
	
	def __createReorderedPhenotypeFile(self, fileName, outFileName, colorDict):
		rColorDict = {'red': 2, 'green': 3, 'blue': 4, 'grey': 8}
		colorList = rColorDict.values()
		colorList.sort()
		colorList.sort()
		fh = ReadFileAtOnceParser(fileName)
		outFh = CsvFileWriter(outFileName)
		header = fh.getSplittedLine()
		outFh.write(header)
		outFh2 = CsvFileWriter(outFileName + '2')
		outFh2.write(header)
		groupIdx, sampleIdx = self.__getGroupAndSampleIdxFromHeader(header)
		lineDict = dict([[splittedLine[sampleIdx].split('_')[0], splittedLine] for splittedLine in fh])
		for sampleName in Utilities.getOrderedKeys(lineDict):
			splittedLine = lineDict[sampleName]
			outFh.write(splittedLine)
			if not splittedLine[groupIdx]:
				continue
			group = splittedLine[groupIdx].split()[-1]
			#color = rColorDict[colorDict[group]]
			#for i in range(len(splittedLine)):
				#if i != 2:		
					#splittedLine[i] = colorList.index(color) + 1
			#outFh2.write([sampleName] + [splittedLine[i] for i in range(len(splittedLine)) if i != 2])
			outFh2.write([sampleName])
		
	def __createHistRscriptFile(self, targetDir, keyword, histDataFile, maxValue, lohHistFileName, centromereDict):
		print 'MAX = [%f]' % maxValue
		#print lohPointList
		maxValue = min(100, maxValue + 5)
		lohStr = lohGraphStr = ''
		if lohHistFileName:
			lohGraphStr = '+ geom_line(data=dataLoh, aes(x=V1, y=V3, color = "black"), stat = "identity", size = 0.3) + scale_colour_manual(values = c("black"), labels = c(""), name = "LOH")'
			lohStr = 'dataLoh = read.table("%s")' % lohHistFileName
		#if lohPointList:
			#lohStr = 'lohDataFrame = data.frame(pos=c(%s), )' % (','.join([str(start) for chrName, start, end, percent in lohPointList]))
		centroPosList = [125000000,93300000,91000000,50400000,48400000,61000000,59900000,45600000,49000000,40200000,53700000,35800000,17900000,17600000,19000000,36600000,24000000,17200000,26500000,27500000,13200000,14700000]
		if centromereDict:
			centroPosList = []
			for chrName in range(1, 23):
				start, end = centromereDict[str(chrName)]
				centroPosList.append((start + end) / 2)
		rFileName = os.path.join(targetDir, '%s_hist_%s.R' % (keyword, self.__suffix))
		rStr = """
		library(ggplot2)
		data = read.table("%(histDataFile)s")
		%(lohStr)s
		subdat1 = subset(data, V1>0)
		subdat2 = subset(data, V1<0)
		subdat2$V1 <- factor(subdat2$V1, levels = levels(factor(subdat2$V1)))
		subdat1$V1 <- factor(subdat1$V1, levels = rev(levels(factor(subdat1$V1))))
		centro = data.frame(pos=%(centroPosStr)s, V2=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22))
		png("%(filePattern)s.png", width=4000, height=1800, res=300)
		cols <- c("6" = "magenta", "5" = "purple", "4" = "royalblue4", "3" = "steelblue4", "2" = "darkolivegreen4", "1" = "mediumseagreen", "-1" = "darkorange1", "-2" = "red", "-3" = "darkred", "-4" = "black")
		if ("expand" %%in%% names(formals(coord_cartesian))) {
			yLim <- coord_cartesian(ylim=c(-%(yMax)d,%(yMax)d), expand = FALSE)
		} else{
			yLim <- coord_cartesian(ylim=c(-%(yMax)d,%(yMax)d))
		}
		gr = ggplot(data=data)+ geom_hline(yintercept=-25, colour="white", size=0.5) + facet_grid(~V2, space="free_x", scales="free_x", labeller=label_value)+ theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),text = element_text(size=15), axis.text.y = element_text(size=15), legend.text = element_text(size=10), legend.position = "bottom", legend.box = "horizontal") + geom_bar(data=subdat1,aes(x=V3, y=V5, fill=factor(V1), width=V4), stat="identity")+geom_bar(data=subdat2,aes(x=V3, y=V5, fill=V1, width=V4), stat="identity") + guides(fill=guide_legend(title="CNV",keywidth = 1.5, keyheight = 1.5, label.position="bottom", label.hjust=0.4, title.vjust=0.7, reverse=TRUE, nrow = 1))+ylab("%%copy number gain / loss, copy neutral LOH")+scale_fill_manual(values = cols, breaks = c(6, 5, 4, 3, 2, 1, -1, -2, -3, -4), labels=c(expression("" >= 6),"5","4","3","2","1","-1","-2","-3",expression("" <= -4))) + yLim +scale_x_continuous(breaks = NULL)+geom_point(aes(x=pos, y=0),centro, size=1.5) + geom_hline(yintercept=0, colour="black", size=0.5) %(lohGraphStr)s + scale_y_continuous(breaks = c(50, 100, -50, -100, 0), labels = c(50, 100, 50, 100, 0))
		plot(gr)
		dev.off()
		""" % {'lohStr': lohStr, 'centroPosStr': R()._getStrFromList(centroPosList), 'yMax': maxValue, 'lohGraphStr': lohGraphStr, 'histDataFile': histDataFile, 'filePattern': rFileName.replace('.R', '')}
		outFh = CsvFileWriter(rFileName)
		outFh.write(rStr)
		return rFileName
	
	def __createHistogramForSample(self, dataList, ploidyDict, sampleName, targetDir, lohDataDict, centromereDict):
		histFileName, maxValue, lohHistFileName = self.__createHistDataFileAndGetMaxValue([(sampleName, dataList)], targetDir, ploidyDict, sampleName, lohDataDict)
		rFileName = self.__createHistRscriptFile(targetDir, sampleName, histFileName, maxValue, lohHistFileName, centromereDict)
		Utilities.mySystem('xvfb-run -w 60 --auto-servernum  %sRscript --vanilla %s' % (self.__binStr, rFileName))
	
	def __createMergedHistogram(self, matrixDict, ploidyDict, targetDir, lohMatrixDict, keyword = '', getMaxValueOnly = False, maxValueToUse = None, centromereDict = None):
		print 'Merged hist for "%s": %d samples: %s' % (keyword, len(matrixDict), matrixDict.keys())
		dataList = []
		#NFPET=['P1T','P22T','P26T','P27T','P47T','P4T','P58T','P60T','P6T','P71T','P72T','P73T','P7T','P8T']
		lohDataDict = {}
		dataSize = None
		for sampleName, chrDict in matrixDict.iteritems():
			sampleName = sampleName.split('.')[0]
			if sampleName[:2] == 'CS':# or sampleName in NFPET:
				continue
			currentDataList = matrixDict[sampleName]
			dataList.append((sampleName, currentDataList))
			if not dataSize:
				dataSize = len(currentDataList)
			if dataSize != len(currentDataList):
				raise NotImplementedError
			if not lohMatrixDict:
				continue
			lohList = lohMatrixDict.get(sampleName)
			if lohList:
				lohDataDict[sampleName] = {}
				for chrPos, ploidy in lohList:
					lohDataDict[sampleName][chrPos] = ploidy
				#lohDataList.append((sampleName, lohList))
		histFileName, maxValue, lohHistFileName = self.__createHistDataFileAndGetMaxValue(dataList, targetDir, ploidyDict, 'merged%s' % keyword, lohDataDict)
		if getMaxValueOnly:
			return maxValue
		if maxValueToUse:
			maxValue = maxValueToUse
		rFileName = self.__createHistRscriptFile(targetDir, 'merged%s' % keyword, histFileName, maxValue, lohHistFileName, centromereDict)
		Utilities.mySystem('xvfb-run --auto-servernum  %sRscript --vanilla %s' % (self.__binStr, rFileName))
		
	def __createHistogramBySample(self, matrixDict, ploidyDict, targetDir, lohMatrixDict, centromereDict):
		cluster = ThreadManager(_getNbAvailableCpus())
		for sampleName, chrDict in matrixDict.iteritems():
			if sampleName[:2] == 'CS' or 'P51T' not in sampleName:
				continue
			sampleName = sampleName.split('.')[0]
			print 'Processing sample %s' % sampleName
			lohDataDict = {}
			lohList = lohMatrixDict.get(sampleName)
			if lohList:
				lohDataDict[sampleName] = {}
				for chrPos, ploidy in lohList:
					lohDataDict[sampleName][chrPos] = ploidy
				print 'LOH'
				print lohDataDict
			cluster.submit(Utilities()._runFunc, self.__createHistogramForSample, [matrixDict[sampleName], ploidyDict, sampleName, targetDir, lohDataDict], os.path.join(targetDir, '%s_hist2_%s' % (sampleName, self.__suffix)))
			#break
		cluster.wait()
		
	def __getGroupAndSampleIdxFromHeader(self, header):
		groupIdx = None
		header = [columnName.lower() for columnName in header]
		for i, columnName in enumerate(header):
			if 'classification' in columnName or columnName == 'group':
				groupIdx = i
				break
		return groupIdx, header.index('sample')
	
	def __getShapeList(self, getRcode = False):
		# check http://www.endmemo.com/program/R/pchsymbols.php for full list of R shapes
		shapeDict = {'diamond': 18, 'circle': 16, 'triangle': 17, 'square': 15, 'triangle2': 25}#, 'plus': 3, 'cross': 4, 'star': 8}
		if getRcode:
			return shapeDict.values()
		return shapeDict.keys()
	
	def _getGroupAndColorCodeDictFromFile(self, fileName, defaultGroupValue = None):
		fh = ReadFileAtOnceParser(fileName)
		header = fh.getSplittedLine()
		isCustom = False
		if 'group' not in header:
			isCustom = True
		groupIdx, sampleIdx = self.__getGroupAndSampleIdxFromHeader(header)
		if groupIdx is None:
			raise NotImplementedError('Could not find group column from header = %s in file %s' % (str(header), fileName))
		groupDict = SimpleMultiDict()
		group2Dict = SimpleMultiDict()
		for splittedLine in fh:
			group = splittedLine[groupIdx]
			if not group:
				if defaultGroupValue:
					group = defaultGroupValue
				else:
					continue
			group = group.split()[-1]
			#groupDict[group] = splittedLine[1].split('_')[0]
			groupDict[group] = splittedLine[sampleIdx].strip()
			if isCustom:
				group2Dict[(4, splittedLine[groupIdx].strip())] = splittedLine[sampleIdx].strip()
				group2Dict[(3, splittedLine[3].strip())] = splittedLine[sampleIdx].strip()
		for group in groupDict:
			sampleList = groupDict.getall(group)
			print group, len(sampleList), sampleList
		if isCustom:
			colorDict = {'Gastrinoma': 'green4', 'NF-PET': 'blue', 'SIET': 'grey', 'Insulinoma': 'red'}
			shapeDict = {'Gastrinoma': 'diamond', 'NF-PET': 'circle', 'SIET': 'triangle', 'Insulinoma': 'square'}
			shapeDict2 = {'Gastrinoma': 18, 'NF-PET': 16, 'SIET': 17, 'Insulinoma': 15}
		else:
			groupList = groupDict.keys()
			colorDict = dict([[groupList[i], color] for i, color in enumerate(Color().getColorList(len(groupList)))])
			shapeDict = dict([[groupList[i], marker] for i, marker in enumerate(self.__getShapeList()[:len(groupList)])])
			shapeList = self.__getShapeList(True)[:len(groupList)]
			shapeDict2 = dict([[groupList[i], shape] for i, shape in enumerate(shapeList)])
		return groupDict, group2Dict, colorDict, shapeDict, shapeDict2
	
	def _createLOHneutralSegments(self, fileName, ploidyFile, write = True):
		fh = ReadFileAtOnceParser(fileName)
		outFh = CsvFileWriter(fileName + '_lohNeutral.txt')
		outFh.write(fh.getSplittedLine())
		ploidyDict = self.__getPloidyDictFromFile(ploidyFile)
		segmentList = []
		for splittedLine in fh:
			sampleName = splittedLine[0].split('.')[0]
			#if sampleName[:2] == 'CS':
				#continue
			if not ploidyDict.has_key(sampleName):
				print 'Passing sample %s' % sampleName
				continue
			ploidy = ploidyDict[sampleName]
			nMin = float(splittedLine[-1])
			nMaj = float(splittedLine[-2])
			if (nMaj, nMin) in [(-ploidy, 0), (ploidy, 0)]:
				if write:
					outFh.write(splittedLine)
				else:
					segmentList.append(splittedLine)
		if not write:
			return segmentList
	
	def _createDendrogram(self, matrixFile, groupDict, colorDict, shapeDict, ploidyFile2, shapeDict2, coeff, keyword = None, defaultEmptyValue = None):
		colorFuncStr = shapeFuncStr = setLabelStr = getShapeFuncStr = legendStr = ''
		colName = 'lab.col'
		adjustY = 'yMax * %f / 3' % coeff
		print groupDict
		#if not adjustShapeSize:
			#adjustY = '1'
		drawShapeFuncStr = '''library(graphics)
library(plotrix)
drawShape <- function (x, shape, color, yMax) {
	downCoeff = %(adjustY)s
	if (shape == "circle") {
		draw.circle(x, -3 * downCoeff, 0.2, border=color, col=color)
	}
	else if (shape == "square") {
		rect(x-.2, -4.3 * downCoeff, x+.2, -1.5 * downCoeff, col=color, border=color)
	}
	else if (shape == "diamond") {
		polygon(c(x, x+.2, x, x-.2), c(-4.3 * downCoeff, -2.9 * downCoeff, -1.5 * downCoeff, -2.9 * downCoeff), col=color, border=color)
	}
	else if (shape == "triangle") {
		polygon(c(x-.2, x+.2, x), c(-4.3 * downCoeff, -4.3 * downCoeff, -1.5 * downCoeff), col=color, border=color)
	}
	else if (shape == "triangle2") {
		polygon(c(x-.2, x+.2, x), c(-1.5 * downCoeff, -1.5 * downCoeff, -4.3 * downCoeff), col=color, border=color)
	}
	#else if (shape == "triangle2") {
	#	plot(x, y, type = "l", ylim = c(-3, 3), main = "rotatexy", col = cols[1], lwd = 2, xlim = c(-1, 7))
	#}
	else {
		stop(paste(c("Shape ", shape, " unhandled")))
	}
}
''' % {'adjustY': adjustY}
		for groupName in groupDict:
			sampleList = groupDict.getall(groupName)
			if self.__useShape:
				color = colorDict[groupName]
				colorFuncStr += 'if (sampleName %%in%% c(%s)) {return ("%s")}\n' % (', '.join(['"%s"' % sampleName.split('_')[0] for sampleName in sampleList]), color)
				shape = shapeDict[groupName]
				shapeFuncStr += 'if (sampleName %%in%% c(%s)) {return ("%s")}\n' % (', '.join(['"%s"' % sampleName.split('_')[0] for sampleName in sampleList]), shape)
		shapeStr = drawShapeStr = ''
		if self.__useShape:
			legendStr = 'legend("topright", "(x,y)", legend=c(%s), col=c(%s), pch=c(%s))' % (', '.join(['"%s"' % groupName for groupName in Utilities.getOrderedKeys(groupDict)]), ', '.join(['"%s"' % colorDict[groupName] for groupName in Utilities.getOrderedKeys(groupDict)]), ', '.join([str(shapeDict2[groupName]) for groupName in Utilities.getOrderedKeys(groupDict)]))
			colName = 'col'
			shapeStr = ', pch = getShapeForSample(label)'
			setLabelStr = 'attr(x, "label") <- ""'
			getShapeFuncStr = 'getShapeForSample <- function(sampleName) {\n	%s\n}' % shapeFuncStr
			drawShapeStr = '''for (sampleName in rownames(m[sample.ord, ])) {
	drawShape(x, getShapeForSample(sampleName), getColorForSample(sampleName), yMax)
	x <- x+1
}'''
		imgFile = matrixFile.replace('.txt', '.png')
		if keyword:
			imgFile = matrixFile.replace('.txt', '_%s.png' % keyword)
		rStr = drawShapeFuncStr + '''a<- read.csv2("%(inputFile)s", h=T, row.names=1, sep = "\t")
b<-as.matrix(a)
c<- dist(b)
d<-hclust(c, method="ward")

getColorForSample <- function(sampleName) {
  %(colorFuncStr)s
}
%(getShapeFuncStr)s
## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
	label <- attr(x, "label")
    %(setLabel)s
    ## set label color to red for A and B, to blue otherwise
    #attr(x, "nodePar") <- list(%(colName)s=getColorForSample(label)%(shapeStr)s)
  }
  return(x)
}

## apply labelCol on all nodes of the dendrogram
d <- dendrapply(as.dendrogram(d), labelCol)
sample.ord <- order.dendrogram(d)
#m <- as.matrix(read.csv2("%(ploidyFile2)s", h=T, sep = "\t", row.names=1))
m <- b

png("%(pngFile)s", width=4000, height=2200, res=300)
plot(d)
xylim <- par("usr")
plotdim <- par("pin")
ymult <- getYmult()
x = 1
yMax = par('usr')[4]
%(drawShapeStr)s
%(legendStr)s
dev.off()''' % {'inputFile': matrixFile, 'pngFile': imgFile, 'colorFuncStr': colorFuncStr, 'getShapeFuncStr': getShapeFuncStr, 'setLabel': setLabelStr, 'colName': colName, 'shapeStr': shapeStr, 'legendStr': legendStr, 'ploidyFile2': ploidyFile2, 'drawShapeStr': drawShapeStr}
		rFileName = FileNameGetter(imgFile).get('R')
		outFh = CsvFileWriter(rFileName)
		outFh.write(rStr)
		outFh.close()
		Utilities.mySystem('xvfb-run -w 5 --auto-servernum  %sRscript --vanilla %s %sout' % (self.__binStr, rFileName, rFileName))
	
	def __createMatrixFile(self, matrixDict, targetDir, keyword = '', ploidyFile = None):
		from stats import getAvgAndStdDevFromList
		outFileName = os.path.join(targetDir, 'matrix_%s%s.txt' % (keyword, self.__suffix))
		outFh = CsvFileWriter(outFileName)
		print '%d in matrixDict' % len(matrixDict)
		outFh.write(['Sample'] + [chrPos for chrPos, ploidy in matrixDict.values()[0]])
		if ploidyFile:
			fh = ReadFileAtOnceParser(ploidyFile)
			header = fh.getSplittedLine()
			groupIdx, sampleIdx = self.__getGroupAndSampleIdxFromHeader(header)
			lineDict = dict([[splittedLine[sampleIdx].split('_')[0], splittedLine] for splittedLine in fh])
		for sampleName in Utilities.getOrderedKeys(matrixDict):
			if sampleName[:2] == 'CS':
				continue
			if ploidyFile:
				splittedLine = lineDict[sampleName]
				if not splittedLine[groupIdx]:
					continue
			dataList = matrixDict[sampleName]
			ploidyLine = [ploidy for chrPos, ploidy in dataList]
			#if len(set(ploidyLine)) == 1:
				#print 'Excluding sample %s' % sampleName
				#continue
			#avg, std = getAvgAndStdDevFromList(ploidyLine)
			#if std < 1:
				#print 'Excluding sample %s based on std' % sampleName
				#continue
			outFh.write([sampleName] + ploidyLine)
		return outFileName
	
	def __getSubsetMatrixForSamples(self, matrixDict, sampleList):
		return dict([[sampleName, matrixDict.get(sampleName.split('_')[0])] for sampleName in sampleList])
	
	def __createMergedHistograms(self, matrixDict, ploidyDict, targetDir, group2Dict, lohMatrixDict, groupDict, centromereDict):
		cluster = ThreadManager(_getNbAvailableCpus())
		matrixDict = self.__getSubsetMatrixForSamples(matrixDict, [sampleName for sampleName in matrixDict if sampleName[:2] != 'CS'])
		paramList = [(None, '')]
		if group2Dict:
			paramList += [(group2Dict.getall((3, 'PET')) + group2Dict.getall((3, 'PET assimilated')), '_PET'),
			     (group2Dict.getall((3, 'SIET')), '_SIET'),
			     (group2Dict.getall((4, 'Insulinoma')), '_Insulinoma'),
			     (group2Dict.getall((4, 'D. Gastrinoma')) + group2Dict.getall((4, 'P. Gastrinoma')), '_Gastrinoma'),
			     (group2Dict.getall((4, 'P. Gastrinoma')), '_P_Gastrinoma'),
			     (group2Dict.getall((4, 'D. Gastrinoma')), '_D_Gastrinoma'),
			     (group2Dict.getall((4, 'G1 NF-PET')) + group2Dict.getall((4, 'G2 NF-PET')), '_NF-PET'),
			     (group2Dict.getall((4, 'G1 NF-PET')), '_NF-PET-G1'),
			     (group2Dict.getall((4, 'G2 NF-PET')), '_NF-PET-G2'),
			     (group2Dict.getall((4, 'G1 SIET')), '_SIET-G1'),
			     (group2Dict.getall((4, 'G2 SIET')), '_SIET-G2')
			]
		else:
			for groupName in groupDict:
				sampleList = groupDict.getall(groupName)
				paramList.append((sampleList, '_%s' % groupName))
		#paramList = [(group2Dict.getall((4, 'Insulinoma')), '_Insulinoma'),]
		#maxValue = None
		#for sampleList, keyword in paramList:
		#	subMatrix = matrixDict
		#	if sampleList:
		#		subMatrix = self.__getSubsetMatrixForSamples(matrixDict, sampleList)
		#	maxValue = max(maxValue, self.__createMergedHistogram(subMatrix, ploidyDict, targetDir, None, keyword, True))
		for sampleList, keyword in paramList:
			subMatrix = matrixDict
			subLohMatrix = lohMatrixDict
			if sampleList:
				subMatrix = self.__getSubsetMatrixForSamples(matrixDict, sampleList)
				subLohMatrix = self.__getSubsetMatrixForSamples(lohMatrixDict, sampleList)
			cluster.submit(self.__createMergedHistogram, subMatrix, ploidyDict, targetDir, subLohMatrix, keyword, False, 95, centromereDict)
		'''
		# all samples
		cluster.submit(self.__createMergedHistogram, matrixDict, ploidyDict, targetDir, lohMatrixDict)
		# PET + PET assimilated col D
		cluster.wait()
		return
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((3, 'PET')) + group2Dict.getall((3, 'PET assimilated'))), ploidyDict, targetDir, '_PET')
		# SIET col D
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((3, 'SIET'))), ploidyDict, targetDir, '_SIET')
		# Insulinoma col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'Insulinoma'))), ploidyDict, targetDir, '_Insulinoma')
		# D. Gastrinoma + P. Gastrinoma col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'D. Gastrinoma')) + group2Dict.getall((4, 'P. Gastrinoma'))), ploidyDict, targetDir, '_Gastrinoma')
		# P. Gastrinoma col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'P. Gastrinoma'))), ploidyDict, targetDir, '_P_Gastrinoma')
		# D. Gastrinoma col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'D. Gastrinoma'))), ploidyDict, targetDir, '_D_Gastrinoma')
		# G1 NF-PET + G2 NF-PET col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'G1 NF-PET')) + group2Dict.getall((4, 'G2 NF-PET'))), ploidyDict, targetDir, '_NF-PET')
		# G1 NF-PET col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'G1 NF-PET'))), ploidyDict, targetDir, '_NF-PET-G1')
		# G2 NF-PET col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'G2 NF-PET'))), ploidyDict, targetDir, '_NF-PET-G2')
		# G1 SIET col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'G1 SIET'))), ploidyDict, targetDir, '_SIET-G1')
		# G2 SIET col E
		cluster.submit(self.__createMergedHistogram, self.__getSubsetMatrixForSamples(matrixDict, group2Dict.getall((4, 'G2 SIET'))), ploidyDict, targetDir, '_SIET-G2')'''
		cluster.wait()
	
	def __createReorderedPhenotypeFileFromMatrixFile(self, matrixFile, ploidyFile, outFileName):
		sampleList = []
		fh = ReadFileAtOnceParser(matrixFile)
		header = fh.getSplittedLine()
		for splittedLine in fh:
			sampleName = splittedLine[0].split('_')[0]
			sampleList.append(sampleName)
		#print len(sampleList), sampleList
		del fh
		fh = ReadFileAtOnceParser(ploidyFile)
		outFh = CsvFileWriter(outFileName)
		outFh.write(fh.getSplittedLine()[2:])
		lineDict = dict([[splittedLine[1].split('_')[0], splittedLine] for splittedLine in fh])
		for sampleName in sampleList:
			outFh.write([sampleName] + lineDict[sampleName][3:])
	
	def __getNextSegmentListInWindow(self, segmentList, sampleName, chrName, windowNb):
		nextSegmentList = []
		start = windowNb * self.__windowSize
		end = (windowNb+1) * self.__windowSize - 1
		windowPos = Position(chrName, start, end)
		while segmentList:
			segment = segmentList.pop(0)
			if segment[:2] != [sampleName, chrName]:
				break
			segmentPos = Position(chrName, segment[2], segment[3])
			overlapPos = windowPos.getOverlapPosition(segmentPos)
			if overlapPos:
				nextSegmentList.append(overlapPos)
				if segmentPos.end > windowPos.end:
					segment[2] = windowPos.end + 1
					segmentList.insert(0, segment)
		return nextSegmentList
	
	def __getCentromereDictFromFile(self, fileName):
		centromereDict = {}
		fh = ReadFileAtOnceParser(fileName)
		for splittedLine in fh:
			chrName, start, end = splittedLine[:3]
			if chrName[:3] == 'chr':
				chrName = chrName[3:]
			start = int(start)
			end = int(end)
			info = centromereDict.get(chrName)
			if info:
				oldStart, oldEnd = info
				start = min(start, oldStart)
				end = max(end, oldEnd)
				if len(set([start, oldStart, end, oldEnd])) != 3:
					print start, oldStart, end, oldEnd
					raise NotImplementedError
			centromereDict[chrName] = start, end
		return centromereDict
	
	def __createFileWithTruncatedCentromereData(self, ascatFile, centromereFile, centromereDict):
		outFileName = FileNameGetter(ascatFile).get('_%s.txt' % os.path.basename(centromereFile).split('.')[0])
		outFh = CsvFileWriter(outFileName)
		fh = ReadFileAtOnceParser(ascatFile)
		outFh.write(fh.getSplittedLine())
		for splittedLine in fh:
			print splittedLine
			pos = OrientedPosition(*(splittedLine[1:4] + ['+']))
			centromerePos = Position(pos.ctgId, *centromereDict[pos.ctgId])
			overlapPos = pos.getOverlapPosition(centromerePos)
			if overlapPos:
				covDiff = Coverage(pos).getCovDiffWithPosition(overlapPos)
				if type(covDiff) == types.TupleType:
					covDiff1, covDiff2 = covDiff
					splittedLine[2] = covDiff1.pos.start
					splittedLine[3] = covDiff1.pos.end
					outFh.write(splittedLine)
					splittedLine[2] = covDiff2.pos.start
					splittedLine[3] = covDiff2.pos.end
				else:
					if not covDiff:
						continue
					splittedLine[2] = covDiff.pos.start
					splittedLine[3] = covDiff.pos.end
			outFh.write(splittedLine)
		return outFileName
	
	def process(self, ascatFile, chrFile, targetDir, ploidyFile, histogram = True, merge = False, dendrogram = False, plotAll = False, centromereFile = None, keyword = None, defaultGroupValue = None, mergeCentromereSegments = None, gcFile = None, platform = None, libDir = None, gw6Dir = None):
		if not os.path.isfile(ascatFile):
			ascatFile = Utilities.getFunctionResultWithCache(os.path.join(ascatFile, 'ascatFile.pyDump'), RunAscat(self.__binDir).process, ascatFile, self.__sampleFile, self.__sampleAliasFile, gcFile, platform, libDir, gw6Dir)
		outFileName, matrixDict, ploidyDict, chrSizeDict, centromereDict = self.__createMatrixFileAndGetDicts(ascatFile, chrFile, targetDir, ploidyFile, centromereFile, mergeCentromereSegments)
		if not ploidyFile:
			ploidyFile = CNView(None, 10, self.__binDir, self.__useShape, self.__sampleFile, self.__sampleAliasFile, self.__groupColumnName)._createPloidyFile(ascatFile, chrFile, targetDir, None, centromereFile, mergeCentromereSegments)
			print 'Using ploidyFile = %s' % ploidyFile
		#sys.exit(1)
		groupDict, group2Dict, colorDict, shapeDict, shapeDict2 = self._getGroupAndColorCodeDictFromFile(ploidyFile, defaultGroupValue)
		if plotAll or histogram:
			lohFileName = ascatFile + '_lohNeutral.txt'
			if not os.path.isfile(lohFileName):
				self._createLOHneutralSegments(ascatFile, ploidyFile)
			self.__sampleName = None
			self.__chrName = None
			lohOutFileName, lohMatrixDict = Utilities.getFunctionResultWithCache(lohFileName + '_%s.pyDump5' % self.__suffix, self.__createSegmentPloidyFile, lohFileName, chrSizeDict, targetDir, dict([[sampleName, 0] for sampleName in ploidyDict]), centromereDict)
			if plotAll or merge:
				self.__createMergedHistograms(matrixDict, ploidyDict, targetDir, group2Dict, lohMatrixDict, groupDict, centromereDict)
			if plotAll or not merge:
				self.__createHistogramBySample(matrixDict, ploidyDict, targetDir, lohMatrixDict, centromereDict)
		if plotAll or dendrogram:
			cluster = ThreadManager(_getNbAvailableCpus())
			if not self.__groupColumnName:
				colorDict['SIET'] = 'black'
			ploidyFile2 = ploidyFile + '_reordered.txt'
			ploidyFile3 = ploidyFile2 + '2'
			Utilities()._runFunc(self.__createReorderedPhenotypeFile, [ploidyFile, ploidyFile2, colorDict], ploidyFile2)
			matrixFile = self.__createMatrixFile(matrixDict, targetDir, None, ploidyFile)
			cluster.submit(self._createDendrogram, matrixFile, groupDict, colorDict, shapeDict, ploidyFile3, shapeDict2, 0.01269454, keyword)
			if not self.__groupColumnName:
				coeff = 0.02296309
				currentMatrixDict = self.__getSubsetMatrixForSamples(matrixDict, groupDict.getall('SIET'))
				matrixFile = self.__createMatrixFile(currentMatrixDict, targetDir, 'SIET')
				currentPhenotypeFile = ploidyFile + '_SIET.txt'
				Utilities()._runFunc(self.__createReorderedPhenotypeFileFromMatrixFile, [matrixFile, ploidyFile, currentPhenotypeFile], currentPhenotypeFile)
				cluster.submit(self._createDendrogram, matrixFile, groupDict, colorDict, shapeDict, currentPhenotypeFile, shapeDict2, coeff)
				
				currentMatrixDict = dict([[sampleName, matrixDict[sampleName.split('_')[0]]] for groupName in set(groupDict.keys()) - set(['SIET']) for sampleName in groupDict.getall(groupName)])
				matrixFile = self.__createMatrixFile(currentMatrixDict, targetDir, 'PET')
				currentPhenotypeFile = ploidyFile + '_PET.txt'
				Utilities()._runFunc(self.__createReorderedPhenotypeFileFromMatrixFile, [matrixFile, ploidyFile, currentPhenotypeFile], currentPhenotypeFile)
				cluster.submit(self._createDendrogram, matrixFile, groupDict, colorDict, shapeDict, currentPhenotypeFile, shapeDict2, coeff)
				cluster.wait()
	
	
def run(options, args):
	if options.all:
		CNView(options.windowSize, options.percentage, options.binDir, options.useShape, options.sampleFile, options.sampleAliasFile, options.groupColumnName).processAll(options.fileName, options.chrFile, options.targetDir, options.ploidyFile, options.percentList, options.baseList, options.histogram, options.merge, options.dendrogram, options.plotAll, options.centromereFile)
	elif options.progName == 'ASCAT':
		RunAscat(options.binDir).process(options.fileName, options.sampleFile, options.sampleAliasFile, options.gcFile, options.platform, options.libDir, options.gw6Dir)
	elif options.progName == 'createFileWithUpdatedPositions':
		RunAscat(options.binDir)._createFileWithUpdatedPositions(options.fileName, options.probeFile, options.targetBuild)
	elif options.progName == 'dendroFeatures':
		RunAscat(options.binDir)._createDendrogramForEachFeature(options.fileName, options.targetDir, options.windowSize, options.percentage, options.fileName2, options.chrFile)
	elif options.progName == 'liftOver':
		RunAscat(options.binDir)._liftOverRawProbeFile(options.fileName, options.targetBuild)
	elif options.progName == 'merge':
		RunAscat(options.binDir)._mergePloidyFileWithSampleInfoFile(options.fileName, options.sampleFile, options.sampleAliasFile)
	elif options.progName == 'ploidy2':
		CNView(options.windowSize, options.percentage, options.binDir, options.useShape, options.sampleFile, options.sampleAliasFile, options.groupColumnName)._createPloidyFile(options.fileName, options.chrFile, options.targetDir, None, options.centromereFile, options.mergeCentromereSegments, options.ploidyFile)
		#CNView(options.windowSize, options.percentage, options.binDir, options.useShape)._createPloidyFile2(options.fileName, options.ploidyFile)
	else:
		CNView(options.windowSize, options.percentage, options.binDir, options.useShape, options.sampleFile, options.sampleAliasFile, options.groupColumnName).process(options.fileName, options.chrFile, options.targetDir, options.ploidyFile, options.histogram, options.merge, options.dendrogram, options.plotAll, options.centromereFile, mergeCentromereSegments = options.mergeCentromereSegments, gcFile = options.gcFile, platform = options.platform, libDir = options.libDir, gw6Dir = options.gw6Dir)

runFromTerminal(__name__, [CommandParameter('a', 'all', CommandParameterType.BOOLEAN),
                           CommandParameter('b', 'binDir', 'string'),
                           CommandParameter('B', 'baseList', CommandParameterType.COMMA_SEP_ID, defaultValue = [int(nb * 1000000) for nb in [0.1, 0.5, 1, 2, 5, 10, 20]]),
                           CommandParameter('c', 'chrFile', 'string'),
                           CommandParameter('C', 'centromereFile', 'string'),
                           CommandParameter('d', 'dendrogram', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('f', 'fileName', 'string'),
                           CommandParameter('F', 'fragmentSize', 'int'),
                           CommandParameter('g', 'gcFile', 'string'),
                           CommandParameter('G', 'groupColumnName', 'string'),
                           CommandParameter('gw6Dir', 'string'),
                           CommandParameter('histogram', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('l', 'libDir', 'string'),
                           CommandParameter('m', 'merge', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('M', 'mergeCentromereSegments', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('fileName2', 'string'),
                           CommandParameter('p', 'percentage', 'float'),
                           CommandParameter('P', 'progName', 'string'),
                           CommandParameter('percentList', CommandParameterType.COMMA_SEP_FLOAT),# defaultValue = [0.5, 1, 2, 5, 10]),
                           CommandParameter('platform', 'string'),
                           CommandParameter('ploidyFile', 'string'),
                           CommandParameter('plotAll', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('probeFile', 'string'),
                           CommandParameter('sampleAliasFile', 'string'),
                           CommandParameter('sampleFile', 'string'),
                           CommandParameter('t', 'targetDir', 'string'),
                           CommandParameter('T', 'targetBuild', 'string'),
                           CommandParameter('u', 'useShape', CommandParameterType.BOOLEAN, defaultValue = False),
                           CommandParameter('w', 'windowSize', 'int')], run)
