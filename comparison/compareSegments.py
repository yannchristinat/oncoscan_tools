#!/usr/bin/python
'''
Author: Yann Christinat
Date: unkown

Small utility to compare two oncoscan results and test if one can be the metastasis of the other.

Note: Most likely buggy... (no tests)

Script to run 100 random comparisons
for (( i=0; i<100; i++ )); do 
randi=$(( RANDOM % 359 )); 
randj=$(( RANDOM % 359 )); 
if [ $randi -ne $randj ]; then 
fileA=`head -$randi ../filelist359.txt | tail -1`; 
fileB=`head -$randj ../filelist359.txt | tail -1`;
python -m oncoscan_tools.comparison.compareSegments A:$fileA B:$fileB | grep '#DATA';
fi
done

'''

import sys

from oncoscan_tools.genome import Segment
from oncoscan_tools.oncoscan import Oncoscan
from numpy.f2py.crackfortran import nameargspattern


class Comparison:
	def __init__(self, nameA, fnA, nameB, fnB):
		'''
		Creates a new instance.
		
		:param nameA: Display name for the first ChAS export file 
		:param fnA: Filename of the first ChAS export file
		:param nameB: Display name for the second ChAS export file 
		:param fnB: Filename of the second ChAS export file
		'''
		self.oncoA = Oncoscan(fnA)
		self.oncoB = Oncoscan(fnB)
		
		self.nameA = nameA
		self.nameB = nameB
		

	def mergeSegments(self, dist):
		'''
		Merge segments within 'dist' base pairs.
		
		:param dist: distance below which segments are merged.
		'''
		self.oncoA.segs = self.oncoA.mergeSegments(dist)
		self.oncoB.segs = self.oncoB.mergeSegments(dist)
	
	
	def compareCN(self, s1cn, s2cn):
		'''
		Compare two copy numbers. 
		
		Returns:
		
		- 'eq' if equal
		- 'gt' if more alterations in s1cn than s2cn (could be more loss or more gain)
		- 'lt' if less alterations in s1cn than s2cn (could be less loss or less gain)
		- 'na' if one is a loss and the other a gain (or vice-versa)
		'''
		if (s1cn>2 and s2cn<2) or (s1cn<2 and s2cn>2): #gain/loss or loss/gain
			return 'na'
		if abs(s1cn-s2cn)<1: #same CN
			return 'eq'
		if abs(s1cn-2) >= abs(s2cn-2)+1:
			return 'gt'
		if abs(s1cn-2) <= abs(s2cn-2)-1:
			return 'lt'
		
	
	def printDiffs(self):
		'''
		Prints to terminal all segments not common to both files or present but with different CN or LOH.
		The main purpose is to compare different versions of ChAS.
		'''
		for chrom in self.oncoA.segs:
			for seg in self.oncoA.segs[chrom]:
				if chrom not in self.oncoB.segs or seg not in self.oncoB.segs[chrom]:
					print str(seg)+' '+str(self.oncoA.segs[chrom][seg])+'\tNot found in '+self.nameB
				elif self.oncoA.segs[chrom][seg] != self.oncoB.segs[chrom][seg]:
					print str(seg)+' '+str(self.oncoA.segs[chrom][seg])+'\tDiff CN/LOH in'+self.nameB+': '+str(self.oncoB.segs[chrom][seg])
		for chrom in self.oncoB.segs:
			for seg in self.oncoB.segs[chrom]:
				if chrom not in self.oncoA.segs or seg not in self.oncoA.segs[chrom]:
					print str(seg)+' '+str(self.oncoB.segs[chrom][seg])+'\tNot found in '+self.nameA
	
	def countDiffs(self):
		'''
		Returns the number of segments not common to both files or present but with different CN or LOH.
		The main purpose is to compare different versions of ChAS.
		'''
		
		ncommon = 0
		for chrom in self.oncoA.segs:
			for seg in self.oncoA.segs[chrom]:
				if chrom in self.oncoB.segs and seg in self.oncoB.segs[chrom] and self.oncoA.segs[chrom][seg] == self.oncoB.segs[chrom][seg]:
					ncommon += 1
		return(ncommon)
	
	
	def getOverlapPercent(self):
		stats={'same':{'eq':0, 'lt':0, 'gt':0,'na':0}, 'smaller':{'eq':0, 'lt':0, 'gt':0,'na':0}, 'unmatched':0}
		for chrom in self.oncoA.segs:
			if chrom in self.oncoB.segs:
				for s1 in self.oncoA.segs[chrom]:
					found_doubleOverlap = dict()
					found_simpleOverlap = dict()
					
					for s2 in self.oncoB.segs[chrom]:
						ab = s1.overlap(s2)
						if ab>=0.9: #same segment
							cndiff = self.compareCN(self.oncoA.segs[chrom][s1]['CN'], self.oncoB.segs[chrom][s2]['CN'])
							ba = s2.overlap(s1)
							if ba>=0.9:
								found_doubleOverlap[cndiff] = 1
							else:
								found_simpleOverlap[cndiff] = 1
					
					if len(found_doubleOverlap)==0 and len(found_simpleOverlap)==0:
						stats['unmatched'] += 1
					else:				
						for d in found_doubleOverlap:
							stats['same'][d] += 1
						for d in found_simpleOverlap:
							stats['smaller'][d] += 1
		
		return stats

	def getOverlapPercent_rev(self):
		tmp = self.oncoA
		self.oncoA = self.oncoB
		self.oncoB = tmp
		
		stats = self.getOverlapPercent()
		
		tmp = self.oncoA
		self.oncoA = self.oncoB
		self.oncoB = tmp
		
		return(stats)
									
###### MAIN #####
if len(sys.argv)!=3 or sys.argv[1]=='-h' or sys.argv[1]=='--help':
	print 'Usage: python compareSegments.py nameA:filepathA nameB:filepathB'
	sys.exit(0)
	
nameA=sys.argv[1].split(':')[0]
fnA=sys.argv[1].split(':')[1]

nameB=sys.argv[2].split(':')[0]
fnB=sys.argv[2].split(':')[1]

comp = Comparison(nameA, fnA, nameB, fnB)
comp.mergeSegments(1000000)


print '\t'.join([str(comp.oncoA.nsegs), str(comp.oncoB.nsegs), str(comp.countDiffs())])
#comp.printDiffs()
stats_AB = comp.getOverlapPercent()
stats_BA = comp.getOverlapPercent_rev()

print stats_AB
print stats_BA

#H1: A -> B
h1 = stats_AB['smaller']['lt'] + stats_AB['smaller']['lt'] + stats_BA['unmatched']
print 'Alterations propre a '+nameB+':\t'+str(h1)+' segments'

#H2: B -> A
h2 = stats_BA['smaller']['lt'] + stats_BA['smaller']['eq'] + stats_AB['unmatched']
print 'Alterations propre a '+nameA+':\t'+str(h2)+' segments'

#H3: A = B
h3 = stats_AB['same']['eq'] + stats_AB['same']['lt'] + stats_AB['same']['gt']
print nameA+' equivalent a '+nameB+'\tSupported by '+str(h3)+' segments'

#H4: A != B
h4 = stats_AB['same']['na'] + stats_AB['smaller']['gt'] + stats_AB['smaller']['na'] + stats_BA['smaller']['gt'] + stats_BA['smaller']['na']
print nameA+' incompatible avec '+nameB+'\tSupported by '+str(h4)+' segments'

#Print data for multisample processing
print '\t'.join(['#DATA', nameA, nameB, str(h2), str(h1), str(h3), str(h4)])

