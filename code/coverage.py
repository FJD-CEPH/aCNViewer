import types

from ObjectBase import SimpleObjectBase
from structure import OrientedPosition, Orientation, Position

class Coverage(SimpleObjectBase):
	# split by Orientation and chr
	def __init__(self, pos, nb = 1, score = 0):
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
		
	def getOverlapPosition(self, other):
		if self.pos is None or other.pos is None:
			return
		return self.pos.getOverlapPosition(other.pos)
	
	def getCovDiffWithPosition(self, pos):
		if pos is None:
			return
		#if pos.ori != self.pos.ori:
			#raise NotImplementedError('Positions should have the same orientation...')
		if self.pos.ctgId != pos.ctgId or (self.pos.start == pos.start and self.pos.end == pos.end):
			return
		if self.pos.start != pos.start:
			start = min(self.pos.start, pos.start)
			end = max(self.pos.start, pos.start) - 1
		else:
			start = min(self.pos.end, pos.end) + 1
			end = max(self.pos.end, pos.end)
		if start > end:
			raise NotImplementedError('Inconsistent position (start > end): start = %d, end = %d.\nself.pos = %s\npos= %s' % (start, end, self.pos.getNiceString(), pos.getNiceString()))
		#print "END=", end, start
		pos1 = OrientedPosition(pos.ctgId, start, end, pos.ori)
		end = max(self.pos.end, pos.end)
		#print "END=", end, pos1.end, self.pos.end, pos.end, end
		if end > pos1.end and self.pos.end != pos.end:
			#print '~' * 20
			#print self.pos, pos
			pos2 = OrientedPosition(pos.ctgId, min(self.pos.end, pos.end) + 1, end, pos.ori)
			#print 'KOKO2'
			return Coverage(pos1, self.nb, self.score), Coverage(pos2, self.nb, self.score)
		#print 'KOKO1'
		return Coverage(pos1, self.nb, self.score)
	
	def __checkDiffCov1IsTupleAndDiffCov2IsNone(self, diffCov1, diffCov2):
		if type(diffCov1) == types.TupleType:
			if diffCov2 is not None:
				raise NotImplementedError('diffCov1 %s is tuple and diffCov2 is not None %s...' % (str((diffCov1[0].getNiceString(), diffCov1[1].getNiceString())), diffCov2.getNiceString()))
			return diffCov1
	
	def __checkDiffCovAndReformatDiffPos(self, diffCov1, diffCov2):
		diffRes = self.__checkDiffCov1IsTupleAndDiffCov2IsNone(diffCov1, diffCov2)
		diffRes2 = self.__checkDiffCov1IsTupleAndDiffCov2IsNone(diffCov2, diffCov1)
		if diffRes:
			#print '1'
			return diffRes
		if diffRes2:
			#print '2'
			return diffRes2
		#print 'KOKO'
		#print diffCov1
		#print diffCov2
		return diffCov1, diffCov2
	
	def merge(self, other):
		#print '-' * 20
		#print 'OTHER', other
		#print 'sle', self
		overlapPos = self.getOverlapPosition(other)
		#print overlapPos
		#print self
		overlapCov = Coverage(overlapPos, self.nb+other.nb, self.score+other.score)
		diffCov1 = self.getCovDiffWithPosition(overlapPos)
		#print '1->', diffCov1
		diffCov2 = other.getCovDiffWithPosition(overlapPos)
		diffCov1, diffCov2 = self.__checkDiffCovAndReformatDiffPos(diffCov1, diffCov2)
		#print diffCov1
		#print diffCov2
		#return
		if diffCov1 and diffCov1 > diffCov2:
			tmp = diffCov1
			diffCov1 = diffCov2
			diffCov2 = tmp
		if diffCov2 and overlapCov > diffCov2:
			tmp = diffCov2
			diffCov2 = overlapCov
			overlapCov = tmp
		return [diffCov1, overlapCov, diffCov2]
	