"""
Unit Tests for preprocess module.
"""

import unittest
import sys
import numpy
import argparse
import pysam

import phylotree
import preprocess

# TODO: Stuff to test:
# def process_reads(samfile, var_pos, min_mq, min_bq):
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
    def setUp(self):
        self.mq = 30
        self.bq = 30
        aln1 = pysam.AlignedSegment()
        aln1.reference_start = 10
        aln1.query_name = 'read1'
        aln1.mapping_quality = 30
        aln1.query_sequence = "AAAAATAAAATAAAAT"
        aln1.query_qualities = [30] * 16
        aln1.cigarstring = '16M'

        aln2 = pysam.AlignedSegment()
        aln2.reference_start = 12
        aln2.query_name = 'read2'
        aln2.mapping_quality = 20
        aln2.query_sequence = "AAAGAAGAAAAG"
        qqual = [33] * 12
        qqual[3] = 20
        aln2.query_qualities = qqual
        aln2.cigarstring = '5M2D7M'

        aln3 = pysam.AlignedSegment()
        aln3.mapping_quality = 0
        aln3.query_name = 'read3'

        self.alns = [aln1, aln2, aln3]

    def test_process_reads_read_obs_simple(self):
        res, _ = preprocess.process_reads(self.alns, [15, 20, 25], 20, 10)
        exp = {'read1':{15:'T', 20:'T', 25:'T'},
               'read2':{15:'G', 20:'G', 25:'G'}}
        self.assertEqual(res, exp)

    def test_process_reads_read_obs_min_map_quality(self):
        res, _ = preprocess.process_reads(self.alns, [15, 20, 25], 25, 10)
        exp = {'read1':{15:'T', 20:'T', 25:'T'}}
        self.assertEqual(res, exp)

    def test_process_reads_read_obs_min_base_quality(self):
        res, _ = preprocess.process_reads(self.alns, [15, 20, 25], 20, 30)
        exp = {'read1':{15:'T', 20:'T', 25:'T'},
               'read2':{20:'G', 25:'G'}}
        self.assertEqual(res, exp)

    def test_process_reads_read_obs_paired_end(self):
        aln1b = pysam.AlignedSegment()
        aln1b.reference_start = 30
        aln1b.query_name = 'read1'
        aln1b.mapping_quality = 30
        aln1b.query_sequence = "AAAAACAAAACAAAAT"
        aln1b.query_qualities = [30] * 16
        aln1b.cigarstring = '16M'
        self.alns.append(aln1b)

        var_pos = [15, 20, 25, 35, 40]

        res, _ = preprocess.process_reads(self.alns, var_pos, 20, 10)
        exp = {'read1':{15:'T', 20:'T', 25:'T', 35:'C', 40:'C'},
               'read2':{15:'G', 20:'G', 25:'G'}}
        self.assertEqual(res, exp)

    def test_process_reads_read_obs_paired_end_overlap(self):
        aln1b = pysam.AlignedSegment()
        aln1b.reference_start = 20
        aln1b.query_name = 'read1'
        aln1b.mapping_quality = 20
        aln1b.query_sequence = "AAAAATAAAACAAAAT"
        aln1b.query_qualities = [30] * 16
        aln1b.cigarstring = '16M'
        self.alns.append(aln1b)

        var_pos = [15, 20, 25, 35]

        res, _ = preprocess.process_reads(self.alns, var_pos, 20, 10)
        exp = {'read1':{15:'T', 25:'T', 35:'T'},
               'read2':{15:'G', 20:'G', 25:'G'}}
        self.assertEqual(res, exp)

    def test_process_reads_read_obs_paired_end_overlap_1bad_base_qual(self):
        aln1b = pysam.AlignedSegment()
        aln1b.reference_start = 20
        aln1b.query_name = 'read1'
        aln1b.mapping_quality = 20
        aln1b.query_sequence = "AAAAATAAAACAAAAC"
        qqual = [30] * 16
        qqual[0] = 5
        aln1b.query_qualities = qqual
        aln1b.cigarstring = '16M'
        self.alns.append(aln1b)

        var_pos = [15, 20, 25, 35]

        res, _ = preprocess.process_reads(self.alns, var_pos, 20, 10)
        exp = {'read1':{15:'T', 20:'T', 25:'T', 35:'C'},
               'read2':{15:'G', 20:'G', 25:'G'}}
        self.assertEqual(res, exp)

    def test_process_reads_base_obs_simple(self):
        _, res = preprocess.process_reads(self.alns, [15, 20, 25], 20, 10)
        exp = {10:{'A':1},
               11:{'A':1},
               12:{'A':2},
               13:{'A':2},
               14:{'A':2},
               15:{'G':1, 'T':1},
               16:{'A':2},
               17:{'A':1},
               18:{'A':1},
               19:{'A':2},
               20:{'G':1, 'T':1},
               21:{'A':2},
               22:{'A':2},
               23:{'A':2},
               24:{'A':2},
               25:{'G':1, 'T':1}}
        self.assertEqual(res, exp)


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


class TestBuildEMMatrix(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args([])
        self.args.verbose = False

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
        self.haps = list("ABCDEFGHI")

    def test_build_em_matrix_simple(self):
        reads = ["1:A,2:C", "1:T,2:C", "3:T,4:T", "2:A,4:T"]
        in_mat = preprocess.build_em_matrix(self.ref, self.phy,
                                            reads, self.haps, self.args)
        r1 = ([(0.01/3) * (0.01/3)] +
             ([0.99 * (0.01/3)] * 8))
        r2 = ([0.99 * (0.01/3)] +
             ([(0.01/3) * (0.01/3)] * 8))
        r3 = ([(0.98) * (0.02/3)] + [(0.02/3) * (0.98)] +
              [(0.02/3) * (0.02/3)] + [(0.02/3) * (0.98)] +
              [0.98 * 0.98] + ([(0.02/3) * (0.98)] * 3) +
              [(0.02/3) * (0.02/3)])
        r4 = ([0.99 * (0.02/3)] + [(0.01/3) * (0.98)] + [(0.01/3) * (0.02/3)] +
              ([(0.01/3) * (0.98)] * 5) + [0.99 * (0.02/3)])
        res_mat = numpy.array([r1, r2, r3, r4])
        self.assertEqual(in_mat.shape, (len(reads), len(self.haps)))
        self.assertTrue(numpy.allclose(in_mat, res_mat))

    def test_build_em_matrix_diff_mut(self):
        # need to implement this first.
        pass


class TestBuildEMInput(unittest.TestCase):
    # Need bam file for this too.
    pass


if __name__ == '__main__':
    unittest.main()
