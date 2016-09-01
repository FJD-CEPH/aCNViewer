''' Team: - Prapat Suriyaphol
	  - Masao Yamaguchi
	  - Victor Renault
				
	Research project: Trans-ethnic Genomics Study of Multigenic Disorders
	                  by Japan/France International Collaboration
						
	Project: Set up a database based system to help the study
				
	Research director: Dr. Fumihiko Matsuda '''

"""
ObjectBase.py
	- 'ObjectBase': Base class for all basic structure class like Position, Sequence, etc.
	- 'CommandKey': the basic class which store all the information about one HTML variable: name and all possible values
last modified : Tuesday 30 November 2004
"""
'test'
import copy
import types
import re
import time

try:
	import Log
except ImportError:
	Log = None

#from Bio.EUtils.MultiDict import *
from commandLineUtils import *
from MultiDict import *

from Utilities import *


def trueEq(obj1, obj2):
	eq = obj1.__eq__
	eq2 = obj2.__eq__
	obj1.__eq__ = None
	obj2.__eq__ = None
	#del obj1.__eq__
	#del obj2.__eq__
	value = obj1 == obj2
	obj1.__eq__ = eq
	obj2.__eq__ = eq2
	return value


class SeqInconsistency:
	def __init__(self, seq):
		self.seq = seq
			

def getFunctionArgList(func):
	if hasattr(func, 'im_func'):
		func = func.im_func
	try:
		code = func.func_code
	except:
		print  func
		raise
	callargs = code.co_argcount
	args = list(code.co_varnames[:callargs])
	return args


def myStr(self):
	if self is None:
		return 'None'
	stringValue = getClassName(self) + ' ('
	if not hasattr(self, '__dict__'):
		if type(self) in [types.ListType, types.TupleType]:
			stringValue = ''
			for elt in self:
				stringValue += myStr(elt) + ', '
			return stringValue[:-2]
		return str(self)
	attrNameList = self.__dict__.keys()
	attrNameList.sort()
	# we remove 'self' from the field list to avoid maximum recursion exception
	if 'self' in attrNameList:
		attrNameList.remove('self')
	for attrName in attrNameList:
		stringValue += attrName + ' = '
		currentAttr = self.__dict__[attrName]
		isList = type(currentAttr) == types.ListType
		isDict = type(currentAttr) == types.DictType
		if not isList and not isDict:
			stringValue += str(currentAttr)
		elif isList:
			stringValue += ObjectString.getStringFromObjectList(currentAttr)
		elif isDict:
			stringValue += ObjectString.getStringFromObjectDict(currentAttr)
		stringValue += ', '
	stringValue = stringValue[0:len(stringValue)-2]
	stringValue += ')'
	return stringValue


class ObjectString:
	@staticmethod
	def str(self):
		return myStr(self)
	
	@staticmethod
	def getStringFromObjectList(list1):
		if list1 is None:
			return str(None)
		stringValue = '['
		for elt in list1:
			stringValue += str(elt) + ',\n'
		stringValue += ']'
		return stringValue
	
	@staticmethod
	def getStringFromObjectDict(dict1):
		if dict1 is None:
			return str(None)
		stringValue = '{'
		for elt in dict1:
			stringValue += str(elt) + ': ' + str(dict1[elt]) + ',\n'
		stringValue += '}'
		return stringValue
	

def init(obj, *args):
	argList = getFunctionArgList(obj.__init__)
	localDict = args[0]
	for arg in argList:
		# case where we have one class 'A' which inherits ObjectBase and another one 'B' which inherits 'A' and which has more
		# parameters than 'A' and which constructor calls the one from 'A': in that case, 'argList' will contain some variables
		# that 'localDict' won't have
		if not localDict.has_key(arg):
			continue
		if arg[-4:] == 'List' and localDict[arg] is None:
			localDict[arg] = []
		elif arg[-4:] == 'Dict' and localDict[arg] is None:
			localDict[arg] = {}
		setattr(obj, arg, localDict[arg])

#class MyObjBase(object):
	#def __eq__(self)

class SimpleObjectBase(object):
	listOfAttributesNotToCompare = None
	
	def __init__(self, *args):
		init(self, *args)
	
	@classmethod
	def getObjWithNoneParam(self):
		argList = getFunctionArgList(self.__init__)
		return self(*[None] * (len(argList)-1))
		
	@classmethod
	def getObjFromParamDict(self, attrDict, obj):
		argList = getFunctionArgList(self.__init__)
		eltNameDict = Utilities.reverseDict(self._eltNameToVarNameDict)
		argList = [obj.getattr(attrDict, eltNameDict[argName]) for argName in argList[1:]]
		obj = self(*argList)
		#print obj
		return obj
		
	def getId(self):
		if hasattr(self, 'id'):
			return self.id
		print self
		raise NotImplementedError
	
	def __checkObjectIsSameClass(self, other):
		if not self.getClassName() == other.getClassName():
			if config1.getboolean('TEST'):
				print 'can not compare %s with %s: ' % (self.__class__, other.__class__)
			raise MemoryError
	
	def __getListOfAttributesToCompare(self):
		attributeListToCompare = set(self.__dict__.keys()) - set(['self'])
		if self.listOfAttributesNotToCompare is not None:
			attributeListToCompare -=  set(self.listOfAttributesNotToCompare)
		return attributeListToCompare
	
	def __isList(self, obj):
		return type(obj) == types.ListType
	
	def __isDict(self, obj):
		return type(obj) == types.DictType or isinstance(obj, MultiDict)
	
	def __cmpAttribute(self, attrName, other, printDiffField = True):
		currentAttr = self.__dict__[attrName]
		otherAttr = other.__dict__[attrName]
		if type(currentAttr) != type(otherAttr):
			if printDiffField:
				print 'koko ', [type(currentAttr), type(otherAttr)], attrName
				print currentAttr
				print otherAttr
			return currentAttr == otherAttr
		if self.__isList(currentAttr):
			if not Utilities.areListsEqual(currentAttr, otherAttr):
				if printDiffField:
					print 'koko ', attrName
				return False
			return True
		if not currentAttr == otherAttr:
			if printDiffField:
				print '-' * 100
				print 'koko:', attrName
				print type(currentAttr), str(currentAttr)
				print type(otherAttr), str(otherAttr)
				print '-' * 100
			return False
		return True
	
	def __eq__(self, other, printDiffField = False):
		''' compare all instance variables '''
		if config1.TRUE_COMPARE == 'True':
			if isinstance(other, SimpleObjectBase):
				other = super(SimpleObjectBase, other)
			if isinstance(self, SimpleObjectBase):
				self = super(SimpleObjectBase, self)
			return self == other
		if other is None or not isinstance(other, SimpleObjectBase):
			if printDiffField:
				print 'other is None'
			return False
		try:
			self.__checkObjectIsSameClass(other)
		except MemoryError:
			return
		attributeListToCompare = self.__getListOfAttributesToCompare()
		for attrName in attributeListToCompare:
			if not self.__cmpAttribute(attrName, other, printDiffField):
				return
		return True
	
	def __ne__(self, other):
		return not self.__eq__(other)
	
	def __str__(self, printType = False):
		stringValue = self.getClassName() + ' ('
		attrNameList = self.__dict__.keys()
		attrNameList.sort()
		# we remove 'self' from the field list to avoid maximum recursion exception
		if 'self' in attrNameList:
			attrNameList.remove('self')
		for attrName in attrNameList:
			stringValue += attrName + ' = '
			currentAttr = self.__dict__[attrName]
			isList = self.__isList(currentAttr)
			isDict = self.__isDict(currentAttr)
			currentVal = None
			if not isList and not isDict:
				if type(currentAttr) == types.StringType:
					currentAttr = "'%s'" % currentAttr
				currentVal = str(currentAttr)
			elif isList:
				currentVal = ObjectString.getStringFromObjectList(currentAttr)
			elif isDict:
				currentVal = ObjectString.getStringFromObjectDict(currentAttr)
			if printType:
				currentVal = '(' + currentVal + ',' + str(type(currentAttr)) + ')'
			stringValue += currentVal + ', '
		stringValue = stringValue[0:len(stringValue)-2]
		stringValue += ')'
		return stringValue

	def getClassName(self):
		return getClassName(self)
	
	def clone(self):
		return copy.deepcopy(self)


class ObjectBase(SimpleObjectBase):
	'''
	define rich comparison methods
	the inherited class defines only 3 methods
		1. __eq__
		2. __gt__
		3. __lt__
	as necessary
	'''
	
	__dictPattern = re.compile('(\w+)Dict$')
	
	@staticmethod
	def getNoneObject(self):
		argList = getFunctionArgList(self.__init__)
		# remove the 'self'
		argList.pop(0)
		argList = [None for arg in argList]
		return self(*argList)
	
	def __gt__(self, other):
		raise NotImplementedError
	
	def __le__(self, other):
		return not self.__gt__(other)
	
	def __lt__(self, other):
		raise NotImplementedError
	
	def __ge__(self, other):
		return not self.__lt__(other)
	
	def runFunctionFromMap(self, map, key, *functionArgs):
		function = map.get(key)
		if not function:
			return
		function(self, *functionArgs)	


	
