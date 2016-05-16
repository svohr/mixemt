"""
Unit Tests for preprocess module.
"""

import unittest
import sys

import phylotree
import preprocess

# Stuff to test:
# class HapVarBaseMatrix(object):
#     def __init__(self, refseq, phylo, mut_wt=0.01, mut_max=0.5):
#     def add_hap_markers(self, hap_var):
#     def _prob(self, hap_pos, pos, base):
#     def prob_for_vars(self, hap, pos_obs):
# def build_em_matrix(refseq, phylo, reads, haplogroups, verbose=False):
# def build_em_input(samfile, refseq, phylo, args):


class TestHapVarBaseMatrix(unittest.TestCase):
    def setUp(self):
        phy_in = ['I, A1G ,,',
                  ',H, A3T A5T ,,',
                  ',,F, A6T ,,',
                  ',,,B, A8T ,,',
                  ',,,C, T5A ,,',
                  ',,G, A7T ,,',
                  ',,,D, A9T ,,',
                  ',,,E, A4T ,,',
                  ',A, A2T A4T ,,']
        self.phy = phylotree.Phylotree(phy_in)
        self.ref = "AAAAAAAAA"

    def test_hapvarbasematrix_init(self):
        mut_wt, mut_max = 0.02, 0.6
        hvb = preprocess.HapVarBaseMatrix(self.ref, self.phy, mut_wt, mut_max)
        self.assertEqual(hvb.refseq, self.ref)
        self.assertEqual(hvb.phylo, self.phy)
        self.assertEqual(hvb.mut_wt, mut_wt)
        self.assertEqual(hvb.mut_max, mut_max)
        markers = {'A':{1:'T', 3:'T', 0:'G'},
                   'B':{0:'G', 2:'T', 4:'T', 5:'T', 7:'T'},
                   'C':{0:'G', 2:'T', 5:'T'},
                   'D':{0:'G', 2:'T', 4:'T', 6:'T', 8:'T'},
                   'E':{0:'G', 2:'T', 3:'T', 4:'T', 6:'T'},
                   'F':{0:'G', 2:'T', 4:'T', 5:'T'},
                   'G':{0:'G', 2:'T', 4:'T', 6:'T'},
                   'H':{0:'G', 2:'T', 4:'T'},
                   'I':{0:'G'}}
        self.assertEqual(hvb.markers, markers)

    def test_probs(self): 
        mut_wt, mut_max = 0.10, 0.10
        hvb = preprocess.HapVarBaseMatrix(self.ref, self.phy, mut_wt, mut_max)
        # hits
        self.assertEqual(hvb._prob(hvb.markers['I'], 0, 'G'), 1.0 - 0.1)
        self.assertEqual(hvb._prob(hvb.markers['I'], 4, 'A'), 1.0 - 0.1)
        self.assertEqual(hvb._prob(hvb.markers['I'], 3, 'A'), 1.0 - 0.1)
        # misses
        self.assertEqual(hvb._prob(hvb.markers['I'], 0, 'A'), 0.1 / 3.0)
        self.assertEqual(hvb._prob(hvb.markers['I'], 4, 'T'), 0.1 / 3.0)
        self.assertEqual(hvb._prob(hvb.markers['I'], 3, 'T'), 0.1 / 3.0)

    def test_probs_mut_wts(self):
        mut_wt, mut_max = 0.10, 0.50
        hvb = preprocess.HapVarBaseMatrix(self.ref, self.phy, mut_wt, mut_max)
        # hits
        self.assertEqual(hvb._prob(hvb.markers['I'], 0, 'G'), 1.0 - 0.1)
        self.assertEqual(hvb._prob(hvb.markers['I'], 4, 'A'), 1.0 - 0.2)
        self.assertEqual(hvb._prob(hvb.markers['I'], 3, 'A'), 1.0 - 0.2)
        # misses
        self.assertEqual(hvb._prob(hvb.markers['I'], 0, 'A'), 0.1 / 3.0)
        self.assertEqual(hvb._prob(hvb.markers['I'], 4, 'T'), 0.2 / 3.0)
        self.assertEqual(hvb._prob(hvb.markers['I'], 3, 'T'), 0.2 / 3.0)

    def test_probs_for_vars(self):
        mut_wt, mut_max = 0.10, 0.10
        hvb = preprocess.HapVarBaseMatrix(self.ref, self.phy, mut_wt, mut_max)
        obs_I = zip(range(9), "GAAAAAAAA")
        self.assertEqual(str(hvb.prob_for_vars('I', obs_I)), str(0.9 ** 9))
        self.assertEqual(str(hvb.prob_for_vars('C', obs_I)), 
                         str((0.9 ** 7) * ((0.1 / 3) ** 2)))
        self.assertEqual(str(hvb.prob_for_vars('D', obs_I)), 
                         str((0.9 ** 5) * ((0.1 / 3) ** 4)))

    def test_probs_for_vars_mut_wts(self):
        mut_wt, mut_max = 0.10, 0.50
        hvb = preprocess.HapVarBaseMatrix(self.ref, self.phy, mut_wt, mut_max)
        obs_I = zip(range(9), "GAAAAAAAA")
        self.assertEqual(str(hvb.prob_for_vars('I', obs_I)), 
                         str((0.9 ** 7) * (0.8 ** 2)))
        self.assertEqual(str(hvb.prob_for_vars('C', obs_I)), 
                         str((0.9 ** 5) * (0.8 ** 2) * ((0.1 / 3) ** 2)))


class TestProcessReads(unittest.TestCase):
    # Need a bam file to test this out on.
    pass


class TestReadSignatures(unittest.TestCase):
    def test_read_signature(self):
        obs_tab = {1:'A', 2:'C', 3:'G', 4:'T'}
        sig = preprocess.read_signature(obs_tab)
        self.assertEqual(sig, "1:A,2:C,3:G,4:T")
        obs_from_sig = preprocess.pos_obs_from_sig(sig)
        self.assertEqual(obs_tab, dict(obs_from_sig))
    
    def test_read_signature_bad_input(self):
        with self.assertRaises(TypeError):
            obs_tab = {'A':'A', 2:'C', 3:'G', 4:'T'}
            preprocess.read_signature(obs_tab)
        with self.assertRaises(TypeError):
            obs_tab = "1:A,2:C,3:G,4:T"
            preprocess.read_signature(obs_tab)

    def test_pos_obs_from_sig(self):
        with self.assertRaises(TypeError):
            sig = "A:A,2:C,3:G,4:T"
            preprocess.read_signature(sig)
        with self.assertRaises(TypeError):
            sig = "1:A,2:C,3:G,4:T,5"
            preprocess.read_signature(sig)

    def test_reduce_reads_simple(self):
        reads = {"read1":{1:'A', 2:'C'}, "read2":{3:'G', 4:'T'}}
        res = preprocess.reduce_reads(reads)
        sigs = {"1:A,2:C":["read1"], "3:G,4:T":["read2"]}
        self.assertEqual(res, sigs)

    def test_reduce_reads_overlap(self):
        reads = {"read1":{1:'A', 2:'C'}, "read2":{3:'G', 4:'T'},
                 "read3":{2:'C', 3:'G'}}
        res = preprocess.reduce_reads(reads)
        sigs = {"1:A,2:C":["read1"], "3:G,4:T":["read2"], "2:C,3:G":["read3"]}
        self.assertEqual(res, sigs)

    def test_reduce_reads_match(self):
        reads = {"read1":{1:'A', 2:'C'}, "read2":{3:'G', 4:'T'},
                 "read3":{2:'C', 1:'A'}}
        res = preprocess.reduce_reads(reads)
        sigs = {"1:A,2:C":["read1", "read3"], "3:G,4:T":["read2"]}
        self.assertEqual(res, sigs)

    def test_reduce_reads_diff_bases(self):
        reads = {"read1":{1:'A', 2:'C'}, "read2":{3:'G', 4:'T'},
                 "read3":{2:'C', 1:'T'}}
        res = preprocess.reduce_reads(reads)
        sigs = {"1:A,2:C":["read1"], "3:G,4:T":["read2"], "1:T,2:C":["read3"]}
        self.assertEqual(res, sigs)


#class TestVariantMethods(unittest.TestCase):
#
#
#class TestPhylotreeSimple(unittest.TestCase):
#    def setUp(self):
#        phy_in = ['I, A1G ,,',
#                  ',H, A3T A5T ,,',
#                  ',,F, A6T ,,',
#                  ',,,B, A8T ,,',
#                  ',,,C, T5A! ,,',
#                  ',,G, A7T ,,',
#                  ',,,D, (A9T) ,,',
#                  ',,,E, A4t ,,',
#                  ',A, A2T a4t ,,']
#        self.phy = phylotree.Phylotree(phy_in)



if __name__ == '__main__':
    unittest.main()
