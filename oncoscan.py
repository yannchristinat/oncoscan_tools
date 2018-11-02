'''
Author: Yann Christinat
Date: 12.03.2018

Module to represent an Oncoscan ChAS export file and compute several metrics.
 
'''
import os, math, copy
from oncoscan_tools.genome import Segment, Chromosome

class Oncoscan:
    '''
    Class to represent an Oncoscan ChAS export file.
    '''
    def __init__(self, segfile):
        '''
        Creates a new instance. Oncoscan region definition (annotfile and chrfile) are defined here.
        
        :param segfile: filename of the ChAS export file.
        '''
        self.annotfile = 'Q:/1.BIOINFORMATIQUE/GitCentralRepo/CNV-tools/CNV-tools/data/OncoScan.na33.r1.annot.csv'
        self.chrfile = 'Q:/1.BIOINFORMATIQUE/GitCentralRepo/CNV-tools/CNV-tools/data/OncoScan.na33.r1.chromStats.tsv'
        
        self.getSegments(segfile)
        self.loadChromRegions()
        

    def getSegments(self, fn):
        '''
        Parses the export file and stores all alterations as segments with an associated LOH or CN tag in self.segs. 
        Directly sets self.nsegs which contains the total number of segments.
        
        Structure of self.segs: {chrom: {segment: {'CN': n, 'LOH': True/False}}}
        '''
        f=open(fn,'r')
        head = f.readline().strip().split('\t')
        iLoc = 0
        while iLoc<len(head) and head[iLoc]!='Full Location':
            iLoc += 1
        if iLoc==len(head):
            raise Exception('No full location column!!')
        
        
        col_cn = head.index('CN State')
        col_cntype = head.index('Type')
        col_pos = head.index('Full Location')
        
        self.segs = dict()
        self.nsegs = 0
        for l in f:
            toks = l.strip().split('\t')
            cn = toks[col_cn]
            cntype = toks[col_cntype]
            chrom = toks[col_pos].split(':')[0]
            pos = toks[col_pos].split(':')[1].split('-')
            
            s = Segment(chrom, int(pos[0]), int(pos[1]))
            if chrom not in self.segs:
                self.segs[chrom] = dict()
            if s not in self.segs[chrom]:
                self.segs[chrom][s]={'CN':2, 'LOH':False}
                self.nsegs +=1
                
            if cntype=='LOH':
                self.segs[chrom][s]['LOH'] = True
            else:
                self.segs[chrom][s]['CN'] = float(cn)

        f.close()
        
    def loadChromRegions(self):
        '''
        Loads the chromosome definition file (chrfile). Calls compChromRegions if the chrfile does not exists.
        
        Note: adds 'chr' to the chromosome name if not already present.  
        '''
        if not os.path.exists(self.chrfile):
            self.compChromRegions()
        
        self.regs = dict()
        f=open(self.chrfile, 'r')
        f.readline() #skip header
        for l in f:
            toks = l.strip().split('\t')
            chrom = toks[0]
            if chrom[:3]!='chr':
                chrom = 'chr'+chrom
            
            self.regs[chrom] = Chromosome(chrom, toks[1], toks[2], toks[3], toks[4])
        return self.regs
    
    def compChromRegions(self):
        '''
        Parses the Oncoscan annotation file (positions of SNPs covered by Oncoscan) and writes the result into the self.chrfile file.
        
        Writes to the self.chrfile ('../data/OncoScan.na33.r1.chromStats.tsv') one line per chromosome with the following columns:
        Chrom, P_Start, P_End, Q_Start, Q_End.
        
        Expects the following columns from the annotation file:
        
        - column 3: chromosome name
        - column 4: genomic position of the SNP
        - column 9: cytoband
        '''
        f = open(self.annotfile,'r')
        l = f.readline()
        while l.strip()[0]=='#':
            l = f.readline()
        l = f.readline()
        
        stats = dict() #variable to store the start/end position of each arm: {chrom: {pstart: x, pend: x, qstart: x, qend: x}}
        for l in f:
            toks = l.strip().split(',')
            chrom = toks[2].replace('"','')
            pos = int(toks[3].replace('"',''))
            cytoband = toks[8].replace('"','')
            
            if chrom not in stats:
                stats[chrom] = {'pstart':None, 'pend':None, 'qstart':None, 'qend':None}
            
            if cytoband[0]=='p':
                if stats[chrom]['pend'] is None or pos>stats[chrom]['pend']:
                    stats[chrom]['pend'] = pos
                if stats[chrom]['pstart'] is None or pos<stats[chrom]['pstart']:
                    stats[chrom]['pstart'] = pos
                    
            if cytoband[0]=='q':
                if stats[chrom]['qend'] is None or pos>stats[chrom]['qend']:
                    stats[chrom]['qend'] = pos
                if stats[chrom]['qstart'] is None or pos<stats[chrom]['qstart']:
                    stats[chrom]['qstart'] = pos
        
        f.close()
        
        f = open(self.chrfile,'w')
        f.write('\t'.join(['Chrom', 'P_Start', 'P_End', 'Q_Start', 'Q_End'])+'\n')
        for c in stats:
            f.write('\t'.join([c, str(stats[c]['pstart']), str(stats[c]['pend']), str(stats[c]['qstart']), str(stats[c]['qend'])])+'\n')
        f.close()
        
    def mergeSegments(self, dist):
        '''
        Merge segments within  'dist' Mb.
        Returns a new dictionary with the merged segments.
        
        :param dist: Number of base pairs at which two segments are merged.
        '''
        merged = copy.deepcopy(self.segs)
        for chrom in merged:
            merged_segs = sorted(merged[chrom].keys())
            i=0
            while i<len(merged_segs)-1:
                si = merged_segs[i]
                sj = merged_segs[i+1]
                
                #The next segment is within 'dist' and both segments have the same CN/LOH -> Merge.
                if si.end+dist>sj.start and merged[chrom][si]['CN']==merged[chrom][sj]['CN'] and merged[chrom][si]['LOH']==merged[chrom][sj]['LOH']: 
                    #Replace the segment i with the merged one and delete the i+1 segment
                    s = Segment(chrom, si.start, sj.end)
                    merged_segs[i] = s
                    del merged_segs[i+1]
                    merged[chrom][s] = merged[chrom][si]
                else:
                    i += 1
                    

            for s in merged[chrom].keys():
                if s not in merged_segs:
                    del merged[chrom][s]

        return merged
    
