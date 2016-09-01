import types
import re

try:
	import Log
except ImportError:
	Log = None

class ValueParser:
	_nbPattern = re.compile('^\d+$')
	
	def getIdxListFromIntervalStrOrList(self, intervalStrOrList, allowFloat = False):
		idxList = []
		if type(intervalStrOrList) in (types.StringType, types.UnicodeType):
			intervalStrOrList = intervalStrOrList.split(',')
		for val in intervalStrOrList:
			if allowFloat and self.isFloat(val):
				val = float(val)
				idxList.append(val)
			elif self.isNb(val):
				val = int(val)
				idxList.append(val)
			elif '-' in val:
				start, end = val.split('-')
				idxList += range(int(start), int(end)+1)
			else:
				print value, val
				raise NotImplementedError
		return idxList
	
	@staticmethod
	def extractIntFromStr(myStr):
		intStr  = ''
		for char in myStr:
			if ValueParser.isNb(char):
				intStr += char
			else:
				break
		if intStr:
			return int(intStr)
	
	@staticmethod
	def isNb(nb):
		if type(nb) == types.IntType:
			return True
		if type(nb) in (types.StringType, types.UnicodeType):
			if nb and nb[0] == '-':
				nb = nb[1:]
			return re.match(ValueParser._nbPattern, nb) is not None

	@staticmethod
	def isFloat(nb):
		if type(nb) in (types.FloatType, types.IntType):
			return True
		if type(nb) in (types.StringType, types.UnicodeType):
			if nb and nb[0] == '-':
				nb = nb[1:]
			partList = nb.split('.')
			if len(partList) not in [1, 2]:
				return
			return set([ValueParser.isNb(part) for part in partList]) == set([True])

	@staticmethod
	def __isNone(value):
		return value is None or (type(value) in (types.StringType, types.UnicodeType) and (value == str(None) or value == ''))

	@staticmethod
	def getFloatValue(value):
		if ValueParser.__isNone(value):
			return value
		return float(value)

	@staticmethod
	def assertValueIsStringAndGetValue(value):
		if type(value) not in (types.StringType, types.NoneType, types.UnicodeType):
			print 'Type of [%s (%s)] should be string' % (value, type(value))
			raise ValueError
		return value

	@staticmethod
	def getIntValue(value):
		if ValueParser.__isNone(value):
			return value
		return int(value)

	@staticmethod
	def getBoolValue(value):
		if ValueParser.__isNone(value):
			return value
		if type(value) in (types.StringType, types.UnicodeType):
			value = value.lower()
			boolDict = {'false': False, 'true': True}
			if boolDict.has_key(value):
				value = boolDict[value]
			else:
				value = int(value)
		return bool(value)

	@staticmethod
	def checkValue(value, commandKeyClass, allowNullValue = True):
		allowedValueList = commandKeyClass.getKeys()
		if allowNullValue:
			allowedValueList += [None]
		if value not in allowedValueList:
			Log.error('value = "' + str(value) + '" not in ' + str(commandKeyClass.getKeys()))
			raise NotImplementedError
		return value
