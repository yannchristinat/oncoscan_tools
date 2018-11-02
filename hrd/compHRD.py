#!/usr/bin/python
'''
Author: Yann Christinat
Date: 12.03.2018

Script to compute the LOH, LST, TD, TDplus, ploidy scores for all Oncoscan ChAS export files within a given directory (first 
argument). Prints the results to the terminal (each line has: filename, DIAMIC, number of segments, LOH, TAI, LST, TD, TDplus, 
ploidy).

Usage: python -m oncoscan_tools.hrd.compHRD BASE_DIR

Note:

- The DIAMIC number is retrieved expected to be the part in the filename before the first '_' character.
- Parses all txt files in the directory. Therefore there should not be any other non-ChAS export files in this directory.

'''
import sys, os, math
from oncoscan_tools.oncoscan import Oncoscan
from oncoscan_tools.genome import Segment

class HRD(Oncoscan):
	'''
	Class to compute different homologous recombination defect (HRD) scores
	'''
	def __init__(self, segfile):
		'''
		Creates a new instance by calling the init method from the Oncoscan class.
		
		:param segfile: filename of the ChAS export file.
		'''
		Oncoscan.__init__(self, segfile)
		#super(HRD, self).__init__(segfile)
		
	def computeLOH(self):
		'''
		Computes the number loss of heterozygosity (LOH) segments. Consider only segments >15Mb but <90% of the chromosome
		length.
		
		Definition by Abkevich et al. (2012): Number of LOH region longer than 15Mb but shorter than the whole chromosome.
		'''
		totLOH = 0
		for chrom in self.segs:
			for s in self.segs[chrom]:
				if len(s)>15*math.pow(10,6) and len(s)<0.9*len(self.regs[chrom]) and self.segs[chrom][s]['LOH']:
					totLOH += 1
		return totLOH
		
	def computeTAI(self):
		'''
		Computes the telomerase allelic imbalance score (TAI). Consider regions >11Mb and whose start or end falls within 1Kb of the 
		start/end of the chromosome. Excludes regions that cross the centromere.
		
		Definition by Birkbak et al. (2012): Number of gain/loss regions that span to the telomere but do not cross the centromere.
		
		Addition by Myriad Inc.: Regions >11Mb.
		'''
		totTAI = 0
		for chrom in self.segs:
			for s in self.segs[chrom]:
				if len(s)<11*math.pow(10,6) or round(self.segs[chrom][s]['CN'])==2:
					continue

				if not self.regs[chrom].hasQ(): #only probes on p arm
					if s.start-self.regs[chrom].p_start < 1000:
						totTAI += 1
				elif not self.regs[chrom].hasP(): #only probes on q arm
					if self.regs[chrom].q_end-s.end < 1000:
						totTAI += 1
				else:	
					#Span to p-telomere
					if s.start-self.regs[chrom].p_start < 1000 and s.end < self.regs[chrom].p_end:
						totTAI += 1
					#Span to q telomere
					elif self.regs[chrom].q_end-s.end < 1000 and s.start > self.regs[chrom].q_start:
						totTAI += 1
		return(totTAI)
		
	def getChromPloidy(self, chrom):
		'''
		Computes the ploidy for a given chromosome. Assumes that there are no overlapping CNV segments!
		
		The function parses through all segments (excluding LOH) and get the minimum copy number that span 90% of the 
		chromosome.
		
		:param chrom: chromosome name
		
		'''
		
		#Returns 2 if no CNV is reported in the chromosome
		if chrom not in self.segs:
			return(2)
		
		#Parses through all deletions and gains within the chromosome and computes the total length of the segments that have a 
		#given number of copies. For instance gains[3] = 123.456Mb means that we have regions with at least 3 copies that sum up 
		#in total to 123.456Mb. 
		dels = dict()
		gains = dict()
		for s in self.segs[chrom]:
			cn = int(round(self.segs[chrom][s]['CN']))
			if cn<2:
				for n in range(cn, 2): #If CN=2, also add the length to CN=1  
					if n not in dels:
						dels[n] = 0
					dels[n] += len(s)
			if cn>2:
				for n in range(3,cn+1): #Also add the length to any copy gain below this CN.
					if n not in gains:
						gains[n] = 0
					gains[n] += len(s)
		
		#Filters out CNs that do not cover more than 90% of the chromosome. I.e. get the minimum copy number that span 90%.  
		validCN = [2]
		for n in dels:
			if dels[n]>0.9*len(self.regs[chrom]):
				validCN.append(n)
		for n in gains:
			if gains[n]>0.9*len(self.regs[chrom]):
				validCN.append(n)
		
		maxGain = max(validCN)
		minDel = min(validCN)
		
		#Small check...
		if maxGain>2 and minDel<2:
			print('ERROR')
			return(None)
		
		if maxGain>2:
			return(maxGain)
		elif minDel<2:
			return(minDel)
		else:
			return(2)
		
		
	def computePloidy(self):
		'''
		Computes the ploidy for the whole genome.
		
		Calls getChromPloidy on all chromosomes covered by Oncoscan except sexual chromosome (hence a normal sample should have 
		44 chromosomes).
		'''
		totChr = 0
		for chrom in self.regs:
			if chrom=='chrX' or chrom=='chrY':
				continue
			
			totChr += self.getChromPloidy(chrom)
			
		return(totChr)
		

	def computeLST(self):
		'''
		Computes the number of large-state transitions (LST). Perform a smoothing with a 1Kb window and prune segments below 3Mb. 
		Returns the number of segments that are >10Mb and within 3Mb of another >10Mb segment with a different CN.
		
		Definition by Popova et al. (2012): Chromosomal break between adjacent regions of at least 10 Mb after removing regions <3Mb.
		'''
		totLST = 0
		for chrom in self.regs:
			if chrom not in self.segs:
				continue
				
			#smooth segments
			segs = list() #keep only segments with an altered copy number (no LOH only segments)
			cnvs = dict()
			for s in self.segs[chrom]:
				if self.segs[chrom][s]['CN']!=2:
					segs.append(s)
					cnvs[s] = self.segs[chrom][s]['CN']
			segs.sort()
			
			i=0
			while i<len(segs)-1:
				si = segs[i]
				sj = segs[i+1]
				if si.end+1000000>sj.start: #The next segment is within 1Mb
					if len(si)<3000000 or len(sj)<3000000 or cnvs[si]==cnvs[sj]: #The next segment is below 3Mb or both segments have the same CN
						#Create new segment
						s = Segment(s.chrom, si.start, sj.end)
						segs[i] = s
						cnvs[s] = cnvs[si]
						
						#debug
						#print 'Merging '+str(si)+' (CN='+str(cnvs[si])+' len='+str(len(si)/1000000)+') and '+str(sj)+' (CN='+str(cnvs[sj])+' len='+str(len(sj)/1000000)+')'
						#print 'Merge results: '+str(s)+' (CN='+str(cnvs[s])+' len='+str(len(s)/1000000)+')'	
						
						del segs[i+1]
						del cnvs[si]
						del cnvs[sj]
					else:
						i += 1
				else:
					i += 1
			
			#Prune segments below 3Mb
			i=0
			while i<len(segs)-1:
				si = segs[i]
				if len(si)<3000000:
					del segs[i]
				else:
					i += 1
					
			#Find breakpoints of regions >=10Mb
			i=0
			while i<len(segs)-1:
				si = segs[i]
				sj = segs[i+1]
				if si.end+1000000>sj.start: #The next segment is within 3Mb
					if cnvs[si]!=cnvs[sj]:
						totLST += 1
				i += 1
			
		return(totLST)
		
	def computeTD(self):
		'''
		Computes the number of tandem repeats (TD). Returns the number of segments <1Mb with copy gains.
		'''
		totTD = 0
		for chrom in self.segs:
			for s in self.segs[chrom]:
				cn = self.segs[chrom][s]['CN']
				if len(s)<1*math.pow(10,6) and cn>2 and cn<5:
					totTD += 1
		return(totTD)
	
	def computeTDplus(self):
		'''
		Computes the "tandem repeat plus" score. Returns the number of segments >1Mb but <10Mb with copy gains.
		'''
		totTDp = 0
		for chrom in self.segs:
			for s in self.segs[chrom]:
				cn = self.segs[chrom][s]['CN']
				if len(s)>=1*math.pow(10,6) and len(s)<10*math.pow(10,6) and cn>2 and cn<5:
					totTDp += 1
		return(totTDp)


if __name__ == '__main__':
	fndir = sys.argv[1]
	print('Filename\tDIAMIC\tNsegs\tLOH\tTAI\tLST\tTD\tTDplus\tPloidy')
	for fn in os.listdir(fndir):
		if fn[-4:]=='.txt':
			onco = HRD(os.path.join(fndir, fn))
			
			bn = os.path.basename(fn)
			diamic = bn.split('_')[0]
			
			#print onco.regs
			loh = str(onco.computeLOH())
			tai = str(onco.computeTAI())
			lst = str(onco.computeLST())
			td = str(onco.computeTD())
			tdplus = str(onco.computeTDplus())
			ploidy = str(onco.computePloidy())
			
			print('\t'.join([bn, diamic, str(onco.nsegs), loh, tai, lst, td, tdplus, ploidy]))






				