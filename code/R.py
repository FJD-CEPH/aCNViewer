import os

from Utilities import Utilities
from ValueParser import ValueParser


class R:

    def __init__(self, binDir=None, waitTime=5, libDir=None, repos=None):
        binDir2 = ''
        if binDir and os.path.isfile(os.path.join(binDir, 'R')):
            binDir2 = binDir.rstrip(os.path.sep) + os.path.sep
        self.__binDir = binDir2
        self.__waitTime = waitTime
        self.__libDir = libDir
        if not repos:
            repos = "http://cran.r-project.org"
        self.__repos = repos

    def getLibraryPathList(self):
        # cmd = '%sRscript --vanilla -e ".libPaths()"' % self.__binDir
        return [dirName.strip().split()[-1].strip('"') for dirName in
                self._execString('.libPaths()', True)]

    def getLibStr(self):
        libStr = ''
        if self.__libDir:
            libStr = ', lib.loc = "%s"' % self.__libDir
        return libStr

    def __isPackageInstalled(self, packageName, libDir=None):
        libStr = ''
        if libDir:
            libStr = "'%s'" % libDir
        rStr = "is.element('%s', installed.packages(%s)[,1])" % (
            packageName, libStr)
        # cmd = '%sRscript --vanilla -e "%s"' % (self.__binDir, rStr)
        # fh = os.popen(cmd)
        fh = self._execString(rStr, True)
        content = fh.read()
        fh.close()
        result = content.split()[1]
        if result not in ['TRUE', 'FALSE']:
            raise NotImplementedError(
                'Unexpected results "%s" from rExpression [%s]' % (result,
                                                                   rStr))
        if result == 'TRUE':
            return True

    def isPackageInstalled(self, packageName):
        isPackagedInstalled = self.__isPackageInstalled(packageName)
        if isPackagedInstalled:
            return True
        if self.__libDir:
            return self.__isPackageInstalled(packageName, self.__libDir)

    def installPackageFromUrl(self, url):
        from WebExtractor import WebExtractor
        content = WebExtractor()._getPageContentForUrl(url)
        fileName = WebExtractor()._getStrIncludedInTag(
            content, '<td> Package&nbsp;source: </td>\n<td> <a href="', '"')
        cmd = 'wget %s' % os.path.join(os.path.dirname(url), fileName)
        Utilities.mySystem(cmd)
        self.installPackage(os.path.basename(fileName))
        Utilities.mySystem('rm %s' % os.path.basename(fileName))

    def installPackage(self, fileName):
        libStr = ''
        if self.__libDir:
            if not os.path.isdir(self.__libDir):
                Utilities.mySystem('mkdir -p %s' % self.__libDir)
            libStr = ", lib = '%s'" % self.__libDir
        if '.tar.gz' not in fileName:
            libStr += ", repos='%s'" % self.__repos
        rStr = "install.packages('%s'%s)" % (fileName, libStr)
        # cmd = '%sRscript --vanilla -e "%s"' % (self.__binDir, rStr)
        # Utilities.mySystem(cmd)
        self._execString(rStr)

    def __getValueStrForList(self, value):
        if 'expression' not in str(value):
            value = '"%s"' % value
        return value

    def _getStrFromList(self, valueList):
        isFloatSet = set([ValueParser().isFloat(value) for value in valueList])
        if valueList and False not in isFloatSet:
            valueList = [str(value) for value in valueList]
        else:
            valueList = [self.__getValueStrForList(
                value) for value in valueList]
        myStr = ', '.join(valueList)
        return 'c(%s)' % myStr

    def _getCmdForFile(self, rFileName):
        waitStr = ''
        if self.__waitTime:
            waitStr = '-w %d' % self.__waitTime
        if self.__binDir:
            cmd = 'cd %s && xvfb-run %s --auto-servernum  %sR CMD BATCH \
--vanilla %s %sout' % (os.path.dirname(rFileName), waitStr, self.__binDir,
                       rFileName, rFileName)
        else:
            cmd = 'cd %s && %sRscript --vanilla %s %sout' % (
                os.path.dirname(rFileName), self.__binDir, rFileName,
                rFileName)
        return cmd

    def _execString(self, rStr, returnStdout=False):
        cmd = '%sRscript --vanilla -e "%s"' % (self.__binDir, rStr)
        if returnStdout:
            return os.popen(cmd)
        Utilities.mySystem(cmd)

    def runCmd(self, cmd, rFileName):
        rFileName = os.path.abspath(rFileName)
        fh = open(rFileName, 'w')
        fh.write(cmd)
        fh.close()
        self.runScript(rFileName)

    def runScript(self, rFileName):
        Utilities.mySystem(self._getCmdForFile(rFileName))
