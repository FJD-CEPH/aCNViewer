import os

from Utilities import Utilities
from ValueParser import ValueParser

class R:
	def __init__(self, binDir = None, waitTime = 5, libDir = None):
		binDir2 = ''
		if binDir:
			binDir2 = binDir.rstrip(os.path.sep) + os.path.sep
		self.__binDir = binDir2
		self.__waitTime = waitTime
		self.__libDir = libDir
		
	def __isPackageInstalled(self, packageName, libDir = None):
		libStr = ''
		if libDir:
			libStr = "'%s'" % libDir
		rStr = "is.element('%s', installed.packages(%s)[,1])" % (packageName, libStr)
		cmd = '%sRscript --vanilla -e "%s"' % (self.__binDir, rStr)
		fh = os.popen(cmd)
		content = fh.read()
		fh.close()
		result = content.split()[1]
		if result not in ['TRUE', 'FALSE']:
			raise NotImplementedError('Unexpected results "%s" from rExpression [%s]' % (result, rStr))
		if result == 'TRUE':
			return True
		
	def isPackageInstalled(self, packageName):
		isPackagedInstalled = self.__isPackageInstalled(packageName)
		if isPackagedInstalled:
			return True
		if self.__libDir:
			return self.__isPackageInstalled(packageName, self.__libDir)
		
	def installPackage(self, fileName):
		libStr = ''
		if self.__libDir:
			if not os.path.isdir(self.__libDir):
				Utilities.mySystem('mkdir -p %s' % self.__libDir)
			libStr = ", lib = '%s'" % self.__libDir
		rStr = "install.packages('%s'%s)" % (fileName, libStr)
		cmd = '%sRscript --vanilla -e "%s"' % (self.__binDir, rStr)
		Utilities.mySystem(cmd)
		
	def _getStrFromList(self, valueList):
		if valueList and ValueParser().isFloat(valueList[0]):
			valueList = [str(value) for value in valueList]
		else:
			valueList = ['"%s"' % value for value in valueList]
		myStr = ', '.join(valueList)
		return 'c(%s)' % myStr
	
	def _getCmdForFile(self, rFileName):
		waitStr = ''
		if self.__waitTime:
			waitStr = '-w %d' % self.__waitTime
		#cmd = 'cd %s && xvfb-run %s --auto-servernum  %sR CMD BATCH --vanilla %s %sout' % (os.path.dirname(rFileName), waitStr, self.__binDir, rFileName, rFileName)
		cmd = 'cd %s && %sRscript --vanilla %s %sout' % (os.path.dirname(rFileName), self.__binDir, rFileName, rFileName)
		return cmd
		
	def runScript(self, rFileName):
		Utilities.mySystem(self._getCmdForFile(rFileName))
	
	