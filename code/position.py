import string
import types
import os

from ObjectBase import *


class Position(ObjectBase):

    def __init__(self, ctgId, start, end, allowNegativeValues=False):
        self.ctgId = ctgId
        self.start = ValueParser.getIntValue(start)
        self.end = ValueParser.getIntValue(end)
        if self.end < self.start:
            raise NotImplementedError(
                'Position is inconsistent: start %s is < end %s' %
                (start, end))
        if not allowNegativeValues and (start and start < 0) or \
                (end and end < 0):
            raise NotImplementedError(
                'Position can not have negative values: start = %s, end = %s' %
                (start, end))

    def getNiceString(self):
        return str((self.ctgId, self.start, self.end))

    def __getPosListFromAssayPos(self):
        posList = PosSampleTable(self.db).getAllPos()
        posList.sort()
        return posList

    @staticmethod
    def getMergedPosList(posList):
        if isinstance(posList[0], types.TupleType):
            posList.sort()
        else:
            Utilities.sortByPosition(posList)
        mergedPosList = []
        currentPos = Position(*posList.pop(0))
        nb = 0.0
        for pos in posList:
            nb += 1
            pos = Position(*pos)
            if pos.isMergeable(currentPos):
                pos = Position(pos.ctgId, min(
                    pos.start, currentPos.start), max(pos.end, currentPos.end))
            else:
                mergedPosList.append(currentPos)
            currentPos = pos
            if nb % 100000 == 0:
                print nb / len(posList), Utilities.getTimeString()
        if not currentPos.isMergeable(mergedPosList[-1]):
            mergedPosList.append(currentPos)
        return mergedPosList

    @staticmethod
    def mergePosWithList(position, posList):
        if position is None or position == NonePos:
            return
        for pos in posList:
            mergedPos = position.getOverlapPosition(pos)
            if mergedPos is not None:
                posList.remove(pos)
                start = min(pos.start, position.start)
                end = max(pos.end, position.end)
                posList.append(Position(pos.ctgId, start, end))
                return
        posList.append(position)

    def setCtgId(self, ctgName, locusCtgId, locusId, objectName):
        if self.ctgId != locusCtgId:
            print "Inconsistency in db for locus [%d]: %s's ctgId [%s] and \
locus ctgId [%d] are different!!!" % \
                  (locusId, objectName, str(self.ctgId), locusCtgId)
            raise NotImplementedError
        self.ctgId = ctgName

    def __len__(self):
        return self.end - self.start + 1

    def __hash__(self):
        ''' xor all values '''
        return hash(self.ctgId) ^ hash(self.start) ^ hash(self.end)

    def getPositionDifference(self, other):
        ''' create a list from the (start, end) of each position, sort the \
list and return the difference between the 3rd and the 2nd element '''
        _tmp = [self.start, self.end, other.start, other.end]
        _tmp.sort()
        return _tmp[2] - _tmp[1]

    @staticmethod
    def comparePosition(self, other):
        # cmp result 1 = more, 0 = equal, -1 = less
        if other is None:
            return 1
        myPos = (self.ctgId, self.start, self.end)
        otherPos = (other.ctgId, other.start, other.end)
        return cmp(myPos, otherPos)

    def __gt__(self, other):
        if other is None:
            return True
        if self.ctgId != other.ctgId:
            return False
        return self.comparePosition(self, other) == 1

    def __lt__(self, other):
        return other.__gt__(self)

    def sortPosition(self, other):
        if self.__gt__(other):
            return (other, self)
        return (self, other)

    def touch(self, pos):
        (first, second) = self.sortPosition(pos)
        return (first.end == second.start - 1)

    def doesOverlapWithPosList(self, posList):
        for pos in posList:
            if self.getOverlapPosition(pos):
                return True

    def getOverlapPosition(self, pos):
        if not pos:
            return
        # different ctg id
        try:
            if self.ctgId != pos.ctgId:
                return None
        except:
            print pos
            raise
        (first, second) = self.sortPosition(pos)
        if (first.start <= second.start and first.end >= second.start):
            result = Position(self.ctgId, second.start,
                              min(first.end, second.end))
            return result
        return None

    def overlap(self, pos):
        if self.getOverlapPosition(pos) is not None:
            return True
        return False

    def isMergeable(self, pos):
        ''' 2 positions are mergable when
                1. same contig Id
                2. position is overlap or touched each other '''
        if self.ctgId != pos.ctgId:
            return False
        return self.overlap(pos) or self.touch(pos)

    def isIncludedIn(self, otherPos):
        if otherPos.ctgId != self.ctgId or self.start < otherPos.start or \
                self.end > otherPos.end:
            return False
        return True

NonePos = Position(None, None, None)


class Orientation:
    FORWARD = '+'
    REVERSE = '-'

    __reverseMap = {FORWARD: REVERSE, REVERSE: FORWARD}

    @staticmethod
    def reverse(ori):
        return Orientation.__reverseMap[ori]


class OrientedPosition(Position):

    def __init__(self, ctgId, start, end, ori):
        super(OrientedPosition, self).__init__(ctgId, start, end)
        self.ori = ori

    def getNiceString(self):
        return str((self.ctgId, self.start, self.end, self.ori))

    def getOverlapPosition(self, pos):
        # if self.ori != pos.ori:
            # return
        pos = super(OrientedPosition, self).getOverlapPosition(pos)
        if pos is None:
            return
        return OrientedPosition(pos.ctgId, pos.start, pos.end, self.ori)

noneOrientedPos = OrientedPosition(None, None, None, None)
