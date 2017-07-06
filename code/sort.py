# import pp
import glob
import sys
import os

from Utilities import Utilities
from FileHandler import FileNameGetter


class Sort:
    def __init__(self, maxMem=None):
        self.__maxMem = maxMem

    def _sort(self, fileName, optionStr, targetFileName=None, scriptName=None):
        optionStr2 = ''
        if self.__maxMem:
            optionStr2 = '-S %s' % self.__maxMem
        fileExt = None
        if len(fileName.split()) > 1:
            fileName = '<(%s)' % fileName
            if not scriptName:
                scriptName = FileNameGetter(targetFileName).get('sh')
        else:
            fileExt = Utilities.getFileExtension(fileName)
        if fileExt == 'gz':
            cmd = 'zcat %s | sort %s %s' % (fileName, optionStr2, optionStr)
        else:
            cmd = 'sort %s %s %s' % (optionStr2, optionStr, fileName)
        cmd = 'export LC_ALL=C && %s' % cmd
        print cmd
        if targetFileName:
            cmd += ' > %s' % targetFileName
            Utilities.mySystem(cmd, scriptName)
        else:
            if scriptName:
                fh = open(scriptName, 'w')
                fh.write('#!/bin/bash\n%s' % cmd)
                del fh
                Utilities.mySystem('chmod +x %s' % scriptName)
                cmd = scriptName
            return os.popen(cmd)
