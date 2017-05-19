import types

from ObjectBase import SimpleObjectBase
from position import OrientedPosition, Orientation, Position

try:
    from operon.byteReader import ByteReader
    from operon.byteWriter import ByteWriter
except ImportError:
    ByteReader = ByteWriter = None


class Coverage(SimpleObjectBase):
    _separator = '|'

    # split by Orientation and chr
    def __init__(self, pos, nb=1, score=0):
        if isinstance(pos, types.StringType):
            byteReader = ByteReader(pos)
            self.pos = OrientedPosition(
                *byteReader.readCString().split(self._separator) +
                [Orientation.FORWARD])
            self.nb = byteReader.readUInt32()
            self.score = byteReader.readFloat32()
        else:
            self.pos = pos
            self.nb = nb
            self.score = score

    def getNbBases(self):
        return len(self) * self.nb

    def __len__(self):
        return len(self.pos)

    def __gt__(self, other):
        otherPos = None
        if other:
            otherPos = other.pos
        return self.pos > otherPos

    def __lt__(self, other):
        otherPos = None
        if other:
            otherPos = other.pos
        return self.pos < otherPos

    def toEntry(self):
        byteWriter = ByteWriter()
        byteWriter.writeCString(self._separator.join(
            [self.pos.ctgId, str(self.pos.start), str(self.pos.end)]))
        byteWriter.writeUInt32(self.nb)
        byteWriter.writeFloat32(self.score)
        return byteWriter.getBytes()

    def getOverlapPosition(self, other):
        if self.pos is None or other.pos is None:
            return
        return self.pos.getOverlapPosition(other.pos)

    def getCovDiffWithPosition(self, pos):
        if pos is None:
            return
        # if pos.ori != self.pos.ori:
            # raise NotImplementedError('Positions should have the same
            # orientation...')
        if self.pos.ctgId != pos.ctgId or (self.pos.start == pos.start and
                                           self.pos.end == pos.end):
            return
        if self.pos.start != pos.start:
            start = min(self.pos.start, pos.start)
            end = max(self.pos.start, pos.start) - 1
        else:
            start = min(self.pos.end, pos.end) + 1
            end = max(self.pos.end, pos.end)
        if start > end:
            raise NotImplementedError('Inconsistent position (start > end): \
start = %d, end = %d.\nself.pos = %s\npos= %s' % (
                                                   start, end, self.pos.getNiceString(), pos.getNiceString()))
        # print "END=", end, start
        pos1 = OrientedPosition(pos.ctgId, start, end, pos.ori)
        end = max(self.pos.end, pos.end)
        # print "END=", end, pos1.end, self.pos.end, pos.end, end
        if end > pos1.end and self.pos.end != pos.end:
            # print '~' * 20
            # print self.pos, pos
            pos2 = OrientedPosition(pos.ctgId, min(
                self.pos.end, pos.end) + 1, end, pos.ori)
            # print 'KOKO2'
            return Coverage(pos1, self.nb, self.score),\
                   Coverage(pos2, self.nb, self.score)
        # print 'KOKO1'
        return Coverage(pos1, self.nb, self.score)

    def __checkDiffCov1IsTupleAndDiffCov2IsNone(self, diffCov1, diffCov2):
        if isinstance(diffCov1, types.TupleType):
            if diffCov2 is not None:
                raise NotImplementedError('diffCov1 %s is tuple and diffCov2 \
is not None %s...' % (str((diffCov1[0].getNiceString(),
                           diffCov1[1].getNiceString())),
                      diffCov2.getNiceString()))
            return diffCov1

    def __checkDiffCovAndReformatDiffPos(self, diffCov1, diffCov2):
        diffRes = self.__checkDiffCov1IsTupleAndDiffCov2IsNone(
            diffCov1, diffCov2)
        diffRes2 = self.__checkDiffCov1IsTupleAndDiffCov2IsNone(
            diffCov2, diffCov1)
        if diffRes:
            # print '1'
            return diffRes
        if diffRes2:
            # print '2'
            return diffRes2
        # print 'KOKO'
        # print diffCov1
        # print diffCov2
        return diffCov1, diffCov2

    #def __getattr__(self, name):
        #try:
            #return super().__getattr__(self, name)
        #except:
            #return getattr(self.pos, name)

    def merge(self, other):
        # print '-' * 20
        # print 'OTHER', other
        # print 'sle', self
        overlapPos = self.getOverlapPosition(other)
        # print overlapPos
        # print self
        overlapCov = Coverage(overlapPos, self.nb +
                              other.nb, self.score + other.score)
        diffCov1 = self.getCovDiffWithPosition(overlapPos)
        # print '1->', diffCov1
        diffCov2 = other.getCovDiffWithPosition(overlapPos)
        diffCov1, diffCov2 = self.__checkDiffCovAndReformatDiffPos(
            diffCov1, diffCov2)
        # print diffCov1
        # print diffCov2
        # return
        if diffCov1 and diffCov1 > diffCov2:
            tmp = diffCov1
            diffCov1 = diffCov2
            diffCov2 = tmp
        if diffCov2 and overlapCov > diffCov2:
            tmp = diffCov2
            diffCov2 = overlapCov
            overlapCov = tmp
        return [diffCov1, overlapCov, diffCov2]


class CoverageMerger:
    def __getFirstOverlappingCovInList(self, cov, covList):
        for i in range(len(covList) - 1, -1, -1):
            cov2 = covList[i]
            if cov2.pos.ctgId != cov.pos.ctgId:
                #print '==>'
                return i+1
            if cov.pos.start >= cov2.pos.start:
                #print '-->'
                return i
        #if cov.pos.start <= cov2.pos.start and cov.pos.end <= cov2.pos.end:
        return 0
        #return 1

    def __insertCoverageAtIdx(self, cov, idx, covList):
        #print '(' *10
        #print 'Inserting', idx, cov, len(covList)
        #for i in covList:
                #print i
        #print ')' * 10
        if idx >= len(covList):
            covList.append(cov)
            return
        try:
            covInList = covList[idx]
        except:
            print 'Exception:'
            print cov
            print len(covList), idx
            for i in range(1, min(10, len(covList) + 1)):
                print 'Cov', covList[-i]
            raise
        leftCov, middleCov, rightCov = covInList.merge(cov)
        #print '%' * 20
        #print 'Merging', cov, covInList
        #print leftCov, middleCov, rightCov
        if [middleCov.pos, leftCov, rightCov] == [None, None, None]:
            if cov.pos < covInList.pos:
                covList.insert(idx, cov)
            else:
                covList.append(cov)
            return
        if middleCov.pos is None:
            middleCov = cov
        if leftCov:
            #print 'lll', idx, len(covList)
            covList[idx] = leftCov
            #print 'll', len(covList)
            idx += 1
            covList.insert(idx, middleCov)
            idx += 1
            #print 'L', idx, len(covList)
        else:
            covList[idx] = middleCov
            idx += 1
            #print 'M', idx, len(covList)
        if rightCov:
            #print '>' * 20
            #print 'Right', rightCov
            #print idx, len(covList)
            self.__insertCoverageAtIdx(rightCov, idx, covList)
        #print '%%%%%%%%%'
        #for c in covList:
                #print c
        #print '??'

    def _mergeCov(self, pos, covList):
        from ObjectBase import getClassName
        cov = pos
        if getClassName(pos) != 'Coverage':
            cov = Coverage(pos)
        #print Utilities.getTimeString(), 'START'
        idx = self.__getFirstOverlappingCovInList(cov, covList)
        #print Utilities.getTimeString(), 'After idx'
        #print '+' * 30
        #print pos, idx, len(covList)
        #print len(covList), str(covList[-1])
        #print '[' * 20
        #for a in covList:
                #print a
        #print ']' * 20
        self.__insertCoverageAtIdx(cov, idx, covList)
        #print ':' * 100
        #for co in covList:
                #print co
        #print '.' * 40
        #self.__checkMerge(covList)
        #print Utilities.getTimeString(), 'After insertion'
        return idx