'''
Author: Yann Christinat
Date: 04.05.2018

Module to test the function in the oncoscan_tools package
 
'''

import unittest
from oncoscan_tools.oncoscan import Oncoscan
from oncoscan_tools.genome import Segment

class DummyOncoscan(Oncoscan):
    def __init__(self):
        self.annotfile = 'Q:/GitCentralRepo/CNV-tools/CNV-tools/data/OncoScan.na33.r1.annot.csv'
        self.chrfile = 'Q:/GitCentralRepo/CNV-tools/CNV-tools/data/OncoScan.na33.r1.chromStats.tsv'
        
        self.loadChromRegions()
        

class TestOncoscanMethods(unittest.TestCase):    

    def test_mergeSegments(self):
        onco = DummyOncoscan()
        
        #Added fake segments
        s1 = Segment('chr1', 100, 200)
        s2 = Segment('chr1', 300, 400)
        s3 = Segment('chr1', 1000, 1200)
        
        ### CASE 1: merge first 2 segments but not third
        onco.segs = {'chr1': {s1: {'CN':3, 'LOH': False}, s2: {'CN':3, 'LOH': False}, s3: {'CN':1, 'LOH': True}}}
        merged = onco.mergeSegments(200)
        
        #Test result
        s1s2 = Segment('chr1', 100, 400)
        expected = {'chr1': {s1s2: {'CN':3, 'LOH': False}, s3: {'CN':1, 'LOH': True}}}
        self.assertEqual(expected, merged, 'FAILED CASE 1')
        
        ### CASE 2: no merge as different LOH
        onco.segs = {'chr1': {s1: {'CN':3, 'LOH': True}, s2: {'CN':3, 'LOH': False}, s3: {'CN':1, 'LOH': True}}}
        merged = onco.mergeSegments(200)
        
        #Test result
        self.assertEqual(onco.segs, merged, 'FAILED CASE 2')
        
        ### CASE 3: merge all segments
        onco.segs = {'chr1': {s1: {'CN':3, 'LOH': False}, s2: {'CN':3, 'LOH': False}, s3: {'CN':3, 'LOH': False}}}
        merged = onco.mergeSegments(1000)
        
        #Test result
        s1s2s3 = Segment('chr1', 100, 1200)
        expected = {'chr1': {s1s2s3: {'CN':3, 'LOH': False}}}
        self.assertEqual(expected, merged, 'FAILED CASE 3')
        
        ### CASE 4: merge all segments but first because of LOH
        onco.segs = {'chr1': {s1: {'CN':3, 'LOH': True}, s2: {'CN':3, 'LOH': False}, s3: {'CN':3, 'LOH': False}}}
        merged = onco.mergeSegments(1000)
        
        #Test result
        s2s3 = Segment('chr1', 300, 1200)
        expected = {'chr1': {s1: {'CN':3, 'LOH': True}, s2s3: {'CN':3, 'LOH': False}}}
        self.assertEqual(expected, merged, 'FAILED CASE 4')
            
        ### CASE 5: merge all segments but first because of CN
        onco.segs = {'chr1': {s1: {'CN':4, 'LOH': False}, s2: {'CN':3, 'LOH': False}, s3: {'CN':3, 'LOH': False}}}
        merged = onco.mergeSegments(1000)
        
        #Test result
        s2s3 = Segment('chr1', 300, 1200)
        expected = {'chr1': {s1: {'CN':4, 'LOH': False}, s2s3: {'CN':3, 'LOH': False}}}
        self.assertEqual(expected, merged, 'FAILED CASE 5')
        
    
### MAIN ###
if __name__ == '__main__':
    unittest.main()
    