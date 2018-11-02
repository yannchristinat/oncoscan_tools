'''
Created on 1 fevr. 2018

@author: YCHR
'''
class Segment:
    '''
    Small class to represent a genomic segment and associated operations.
    '''
    def __init__(self, chrom, start, end):
        '''
        Creates a new segment. (Does not check if start<=end!)
        
        :param chrom: chromosome.
        :param start: first position of the segment.
        :param end: last position of the segment.
        '''
        self.chrom = chrom
        self.start = start
        self.end = end
        
    def __str__(self):
        '''
        Returns a string representation of the segment ("chrom:start-end")
        '''
        return self.chrom+':'+str(self.start)+'-'+str(self.end)
    
        
    def __hash__(self):
        '''
        Returns a unique of the segment (based on the string representation).
        '''
        return hash(str(self))
        
    def __len__(self):
        '''
        Returns the length of the segment
        '''
        return self.end-self.start+1
        
    def __eq__(self, s):
        '''
        Tests if two segment are identical.
        '''
        return self.chrom==s.chrom and self.start==s.start and self.end==s.end

    def __ne__(self, s):
        '''
        Tests if two segment are not identical.
        '''
        return(not (self==s))
        
    def __gt__(self, s):
        '''
        Tests if the segment is greater than segment s.
        
        Returns true if chromosome number is higher or (if same chromosome) segment start is higher.
        '''
        if self.chrom>s.chrom:
            return(True)
        elif self.chrom<s.chrom:
            return(False)
        else:
            return(self.start>s.start)
        
    def __ge__(self, s):
        '''
        Returns true if the segment is identical of greater than s.
        '''
        return(self>s or self==s)
    
    def __lt__(self, s):
        '''
        Returns true if the segment is identical of smaller than s.
        '''
        return(s>self and self!=s)
        
    def __le__(self, s):
        '''
        Returns true if the segment is smaller than s.
        '''
        return(self<s or self==s)
        
    def overlap(self, s):
        '''
        Returns true if the two segment are on the same chromosome and overlap.
        '''
        if self.chrom!=s.chrom: 
            return 0
        if self.start>s.end or self.end<s.start:
            return 0

        start = self.start
        if self.start<s.start:
            start = s.start
        end = self.end
        if self.end>s.end:
            end = s.end
        return float(end-start+1)/len(self)
    
    def getCytoband(self, chroms):
        '''
        Returns the chromosome arm (p or q) on which is situated the segment. Returns 'pq' if the segment spans both arms.
        
        :param chroms: list of Chromosome objects 
        '''
        if self.chrom not in chroms:
            return(None)
        cyto = ''
        if chroms[self.chrom].hasP() and self.start<chroms[self.chrom].p_end:
            cyto += 'p'
        if chroms[self.chrom].hasQ() and self.end>chroms[self.chrom].q_start:
            cyto += 'q'
        return(cyto)
            

class Chromosome:
    '''
    Small class to represent a chromosome as covered by Oncoscan.
    '''
    
    def __init__(self, name, pstart, pend, qstart, qend):
        '''
        Creates a new instance.
        
        :param name: chromosome name
        :param pstart: first genomic position of the p arm
        :param pend: last genomic position of the p arm
        :param qstart: first genomic position of the q arm
        :param qend: last genomic position of the q arm 
        '''
        
        self.name = name
    
        if pstart=='None':
            self.p_start = None
        else:
            self.p_start = int(pstart)
        
        if pend=='None':
            self.p_end = None
        else:
            self.p_end = int(pend)
        
        if qstart=='None':
            self.q_start = None
        else:
            self.q_start = int(qstart)
        
        if qend=='None':
            self.q_end = None
        else:
            self.q_end = int(qend)
        
    def hasP(self):
        '''
        Returns true if the chromosome has a p arm.
        '''
        return(self.p_start is not None)

    def hasQ(self):
        '''
        Returns true if the chromosome has a q arm.
        '''
        return(self.q_start is not None)
        
    def lenP(self):
        '''
        Return the length of the p arm.
        '''
        lenP = 0
        if self.hasP():
            lenP = self.p_end - self.p_start +1
        return(lenP)
            
    def lenQ(self):
        '''
        Return the length of the p arm.
        '''
        lenQ = 0
        if self.hasQ():
            lenQ = self.q_end - self.q_start +1
        return(lenQ)

    def __len__(self):
        '''
        Returns the length of the chromosome (length of p + q).
        '''
        return(self.lenP()+self.lenQ())
        
    def __str__(self):
        '''
        Returns a string representation of the chromosome as "name: p=start-end, q=start-end"
        '''
        if self.hasP() and self.hasQ():
            return(self.name+': p='+str(self.p_start)+'-'+str(self.p_end)+', q='+str(self.q_start)+'-'+str(self.q_end))
        if not self.hasP():
            return(self.name+': p=not_covered, q='+str(self.q_start)+'-'+str(self.q_end))
        if not self.hasQ():
            return(self.name+': p='+str(self.p_start)+'-'+str(self.p_end)+', q=not_covered')