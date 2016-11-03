import types
import re

try:
	import Log
except ImportError:
	Log = None

from optparse import OptionParser
from ValueParser import ValueParser
from Utilities import getLastErrorMessage

''' This class is designed to contain the different values associated to one variable. Here the variable name is 'name' and to define
its different possible values, it is sufficient to define them at the same level than 'name'. It is assumed that all the variable which
name are upper case are one possible values for the variable '''
class CommandKey:
	name = None
	
	def __init__(self, name):
		self.name = name
	
	# by default, we return all the variables of the class which are in upper case
	def getKeys(self):
		keys = []
		pattern = re.compile('^([A-Z_0-9]+)$')
		for variableName in self.__dict__.keys():
			if re.match(pattern, variableName):
				keys.append(self.__dict__[variableName])
		return keys
	getKeys = classmethod(getKeys)


class CommandParameterType(CommandKey):
	INT = 'int'
	FLOAT = 'float'
	STRING = 'string'
	BOOLEAN = 'bool'
	COMMA_SEP_ID = 'comma_sep_id'
	COMMA_SEP = 'comma_sep'
	COMMA_SEP_FLOAT = 'comma_sep_float'
	DB = 'db'
	DB_LIST = 'db_list'

class CommandParameter:
	_typeDict = {CommandParameterType.BOOLEAN: CommandParameterType.INT, CommandParameterType.DB: CommandParameterType.STRING,
				 CommandParameterType.COMMA_SEP: CommandParameterType.STRING, CommandParameterType.COMMA_SEP_ID: CommandParameterType.STRING,
				 CommandParameterType.DB_LIST: CommandParameterType.STRING, CommandParameterType.COMMA_SEP_FLOAT: CommandParameterType.STRING}
	
	def __init__(self, optionName, definition, commandType = None, varName = None, defaultValue = None, helpString = '',
	             action = 'store'):
		if not commandType:
			commandType = definition
			definition = optionName
			optionName = None
		if varName is None:
			varName = definition
		if type(varName) != types.StringType:
			defaultValue = varName
			varName = definition
		ValueParser.checkValue(commandType, CommandParameterType)
		self.originalType = commandType
		commandType = self._typeDict.get(commandType, commandType)
		self.optionName = optionName
		self.definition = definition
		self.commandType = commandType
		self.varName = varName
		self.defaultValue = defaultValue
		self.helpString = helpString
		self.action = action
	
	def addOption(self, optionParser):
		optList = ["--%s" % self.definition]
		if self.optionName:
			optList.insert(0, "-%s" % self.optionName)
		optionParser.add_option(*optList, action = self.action, type = self.commandType, dest = self.varName,
		                        default = self.defaultValue, help = self.helpString)
	

def __getDbFromStr(value):
	try:
		from Database import Database
	except ImportError:
		class Database:
			def __init__(self, dbName, dbUser, dbPassword, host = None):
				self.dbName = dbName
				self.dbUser = dbUser
				self.dbPassword = dbPassword
				self.host = host

	from config import config1
	if value is None:
		value = config1.DBNAME
	host = config1.DBHOST
	if value and '@' in value:
		value, host = value.split('@')
	if value:
		value = Database(value, config1.DBUSER, config1.DBPASSWORD, host = host)
	return value

def runFromTerminal(name, commandParameterList, mainFunctionToRun):
	if name != '__main__':
		return
	try:
		import psyco
		psyco.full()
	except ImportError:
		pass
	parser = OptionParser()
	for commandParameter in commandParameterList:
		commandParameter.addOption(parser)
	(options, args) = parser.parse_args()
	for commandParameter in commandParameterList:
		value = getattr(options, commandParameter.varName)
		if commandParameter.originalType == CommandParameterType.BOOLEAN:
			if value is not None:
				value = bool(value)
		elif commandParameter.originalType == CommandParameterType.DB:
			value = __getDbFromStr(value)
		elif commandParameter.originalType == CommandParameterType.DB_LIST:
			valueList = []
			for val in value.split(','):
				valueList.append(__getDbFromStr(val))
			value = valueList
		elif commandParameter.originalType in (CommandParameterType.COMMA_SEP_ID, CommandParameterType.COMMA_SEP, CommandParameterType.COMMA_SEP_FLOAT):
			if value is None:
				value = []
			else:
				if type(value) != types.ListType:
					value = value.split(',')
				if commandParameter.originalType == CommandParameterType.COMMA_SEP_ID:
					value = ValueParser().getIdxListFromIntervalStrOrList(value)
				elif commandParameter.originalType == CommandParameterType.COMMA_SEP_FLOAT:
					value = ValueParser().getIdxListFromIntervalStrOrList(value, True)
		setattr(options, commandParameter.varName, value)
	try:
		mainFunctionToRun(options, args)
	except:
		if Log:
			Log.error(getLastErrorMessage())
		else:
			print getLastErrorMessage()
		raise