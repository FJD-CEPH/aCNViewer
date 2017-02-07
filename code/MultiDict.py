import types
import re


def getClassName(self):
    if isinstance(self, types.TypeType):
        className = str(self)
    else:
        className = str(self.__class__)
    if className[:7] == '<class ':
        pattern = "^\<class '([\w|\_]+\.)*([\w|\_]+)'\>$"
    elif className[:6] == '<type ':
        pattern = "^\<type '(\w+)'\>$"
    else:
        pattern = '(\w+\.)*([\w|\.]+)'
        # raise NotImplementedError('className [%s] not handled' % className)
    matchObj = re.match(pattern, className)
    if matchObj is None:
        print 'Warning: Matching error [%s]' % className
        return str(type(self))
    return matchObj.groups()[-1]

try:
    from collections import defaultdict
    # raise ImportError
except ImportError:
    class defaultdict(dict):

        def __init__(self, eltType):
            self.__eltType = eltType

        def __getitem__(self, key):
            value = self.get(key)
            if value is None:
                value = self.__eltType()
                self[key] = value
            return value


class SimpleMultiDict(object):

    def __init__(self, valueList=None):
        if valueList is None:
            valueList = []
        self.__valDict = {}
        for key, value in valueList:
            self[key] = value

    def __delitem__(self, key):
        del self.__valDict[key]

    def __iter__(self):
        return self.__valDict.__iter__()

    def __len__(self):
        return len(list(self.allvalues()))

    def iteritems(self):
        for key, value in self.__valDict.iteritems():
            yield key, value[-1]

    def allvalues(self):
        for valueList in self.__valDict.values():
            for value in valueList:
                yield value

    def allitems(self):
        for key, valueList in self.__valDict.items():
            for value in valueList:
                yield key, value

    def items(self):
        return [(key, valueList[-1]) for key, valueList in
                self.__valDict.items()]

    def values(self):
        return [valueList[-1] for valueList in self.__valDict.values()]

    def __str__(self):
        return getClassName(self) + '(' + str(self.__valDict.items()) + ')'

    def __getitem__(self, key):
        valueList = self.__valDict.get(key, [None])
        if len(valueList):
            return valueList[-1]

    def __setitem__(self, key, value):
        currentValue = self.get(key)
        if currentValue is None:
            self.__valDict[key] = []
        self.__valDict[key].append(value)

    def get(self, key, defaultValue=None):
        if key in self:
            return self[key]
        return defaultValue

    def keys(self):
        return self.__valDict.keys()

    allkeys = keys

    def has_key(self, key):
        # return key in self.keys()
        return key in self.__valDict

    def getall(self, key, default=None):
        if default is None:
            default = []
        return self.__valDict.get(key, default)

    def update(self, otherMultiDict):
        iterAllItems = otherMultiDict.items
        if not isinstance(otherMultiDict, types.DictType):
            iterAllItems = otherMultiDict.allitems
        for key, value in iterAllItems():
            self.__setitem__(key, value)


class MultiDict(SimpleMultiDict):

    def __init__(self, objList=None):
        self.__nbObj = 0
        self.orderedKeyList = []  # to keep the original order of the keys
        super(MultiDict, self).__init__(objList)

    def __eq__(self, obj1):
        for key in self.keys():
            if self.getall(key) != obj1.getall(key):
                return False
        return len(self) == len(obj1)

    def removeValueForKey(self, value, key):
        valueList = self.getall(key)
        valueList.remove(value)
        self.order_data.remove((key, value))

    def __setitem__(self, key, value):
        super(MultiDict, self).__setitem__(key, value)
        if key not in self.orderedKeyList:
            self.orderedKeyList.append(key)
        self.__nbObj += 1

    def __delitem__(self, key):
        objLength = len(self.getall(key))
        super(MultiDict, self).__delitem__(key)
        self.orderedKeyList.remove(key)
        self.__nbObj -= objLength

    def __len__(self):
        return self.__nbObj
