"""
Unit Tests for assemble module.
"""

import unittest
import sys
import numpy
import argparse
import collections
import pysam

from mixemt import observe
from mixemt import assemble
from mixemt import phylotree
from mixemt import preprocess

# TODO: Write tests for the follwing functions.
# write_haplotypes(samfile, contrib_reads, reads, read_sigs, prefix, verbose)


class TestContributors(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args([])
        self.args.verbose = False
        self.args.min_reads = 1
        self.args.min_var_reads = 1
        self.args.var_fraction = 0.5
        self.args.var_count = None
        self.args.var_check = False
        self.args.contributors = None

        phy_in = ['I, A1G ,,',
                  ',H, A3T A5T ,,',
                  ',,F, A6T ,,',
                  ',,,B, A8T ,,',
                  ',,,C, T5A ,,',
                  ',,G, A7T ,,',
                  ',,,D, A9T ,,',
                  ',,,E, A4T ,,',
                  ',A, A2T A4T ,,']
        self.ref = "AAAAAAAAA"
        self.phy = phylotree.Phylotree(phy_in, refseq=self.ref)

        self.cons = [['A', 0.4], ['E', 0.3]]
        self.obs = observe.ObservedBases()
        self.obs.obs_tab[1]['T'] = 1
        self.obs.obs_tab[3]['T'] = 2
        self.obs.obs_tab[0]['G'] = 1
        self.obs.obs_tab[6]['T'] = 1
        self.obs.obs_tab[2]['T'] = 1
        self.obs.obs_tab[4]['T'] = 1

        self.wts = numpy.array([1, 1, 1])
        self.haps = list('ABCDEFGHI')
        self.props = numpy.array([0.40, 0.01, 0.01, 0.01, 0.3,
                                  0.01, 0.01, 0.01, 0.01])
        self.mix_mat = numpy.array([
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.01, 0.01, 0.01, 0.01, 0.91, 0.01, 0.01, 0.01, 0.01]])
        self.em_results = (self.props, self.mix_mat)

    def test_get_contributors_no_phy_vars(self):
        self.args.var_check = False
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.assertEqual(res, exp)

    def test_get_contributors_no_phy_vars_rm(self):
        self.args.var_check = False
        self.args.min_reads = 2
        self.args.min_var_reads = 2
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'A', 0.40]]
        self.assertEqual(res, exp)

    def test_get_contributors_with_phy_vars(self):
        self.args.var_check = True
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.assertEqual(res, exp)

    def test_get_contributors_with_phy_vars_rm(self):
        self.args.var_check = True
        self.args.min_reads = 2
        self.args.min_var_reads = 2
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = []
        self.assertEqual(res, exp)

    def test_get_contributors_sorted(self):
        self.props[4] = 0.6
        self.args.var_check = False
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'E', 0.60], ['hap2', 'A', 0.4]]
        self.assertEqual(res, exp)

    def test_get_contributors_manual(self):
        self.args.contributors = "A,E"
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.assertEqual(res, exp)

    def test_get_contributors_manual_weird_choice(self):
        self.args.contributors = "E,F"
        res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                        self.wts, self.em_results, self.args)
        exp = [['hap1', 'E', 0.3], ['hap2', 'F', 0.01]]
        self.assertEqual(res, exp)

    def test_get_contributors_manual_bad_choice(self):
        with self.assertRaises(ValueError):
            self.args.contributors = "E,Z"
            res = assemble.get_contributors(self.phy, self.obs, self.haps,
                                            self.wts, self.em_results,
                                            self.args)

    def test_find_contribs_from_reads_wts_all_one(self):
        res = assemble._find_contribs_from_reads(self.mix_mat, self.wts,
                                                 self.args)
        exp = [0, 4]
        self.assertEqual(res, exp)

    def test_find_contribs_from_reads_wts_all_one_min_reads(self):
        self.args.min_reads = 2
        res = assemble._find_contribs_from_reads(self.mix_mat, self.wts,
                                                 self.args)
        exp = [0]
        self.assertEqual(res, exp)

    def test_find_contribs_from_reads_wts_save_min_reads(self):
        self.args.min_reads = 2
        self.wts = [1, 1, 2]
        res = assemble._find_contribs_from_reads(self.mix_mat, self.wts,
                                                 self.args)
        exp = [0, 4]
        self.assertEqual(res, exp)

    def test_check_contrib_phy_vars_no_rm(self):
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons, res)

    def test_check_contrib_phy_vars_rm_none(self):
        self.cons.append(['C', 0.1])
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)

        self.assertNotEqual(self.cons, res)
        self.assertEqual(self.cons[0:2], res)

    def test_check_contrib_phy_vars_no_rm_one_base(self):
        self.cons.append(['C', 0.1])
        self.obs.obs_tab[5]['T'] = 1
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons, res)

    def test_check_contrib_phy_vars_rm_shared(self):
        # Should remove E, A4T already seen in A, does not count for E.
        self.obs.obs_tab[6]['T'] = 0
        self.obs.obs_tab[2]['T'] = 0
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons[0:1], res)

    def test_check_contrib_phy_vars_empty_obs(self):
        # no observations, no keepers.
        self.obs.obs_tab = collections.defaultdict(collections.Counter)
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual([], res)

    def test_check_contrib_phy_vars_high_min_reads(self):
        # required number of reads too high.
        self.args.min_var_reads = 10
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual([], res)

    def test_check_contrib_phy_vars_high_min_fraction_requirement(self):
        self.args.var_fraction = 0.9
        self.cons.append(['C', 0.1])
        self.obs.obs_tab[4]['A'] = 1
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons[0:2], res)

    def test_check_contrib_phy_vars_high_var_count_requirement(self):
        # required number of variants too high.
        self.args.var_count = 1
        self.args.var_fraction = 0.9
        self.cons.append(['C', 0.1])
        self.obs.obs_tab[5]['T'] = 1
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons, res)

    def test_check_contrib_phy_vars_check_ancestral(self):
        # should pass if A isn't present than 4 can be ancestral evidence
        self.cons.append(['C', 0.1])
        self.obs.obs_tab[4]['A'] = 1
        self.obs.obs_tab[1]['T'] = 0
        self.obs.obs_tab[3]['T'] = 0
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons[1:], res)

    def test_check_contrib_phy_vars_ignore_ancestral(self):
        # should not find C, T5A is explained by 'A'
        # self.cons.append(['C', 0.1])
        self.obs.obs_tab[4]['A'] = 1
        res = assemble._check_contrib_phy_vars(self.phy, self.obs,
                                               self.cons, self.args)
        self.assertEqual(self.cons, res)


class TestAssignReads(unittest.TestCase):
    def setUp(self):
        self.haps = list('ABCDEFGHI')
        self.props = numpy.array([0.40, 0.01, 0.01, 0.01, 0.3,
                                  0.01, 0.01, 0.01, 0.01])
        self.mix_mat = numpy.log(numpy.array([
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.30, 0.01, 0.01, 0.01, 0.40, 0.01, 0.01, 0.01, 0.01],
            [0.01, 0.01, 0.01, 0.01, 0.91, 0.01, 0.01, 0.01, 0.01]]))
        self.em_results = (self.props, self.mix_mat)
        self.cons = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.reads = [['A', 'B'], ['C'], ['D'], ['E', 'F', 'G']]

    def test_find_best_n_for_read_n2(self):
        prob = numpy.array([0.1, 0.2, 0.1, 0.3, 0.9, 0.1])
        con_i = [1,3,5]
        res = assemble._find_best_n_for_read(prob, con_i, 2)
        self.assertEqual(res, [3,1])

    def test_find_best_n_for_read_n3(self):
        prob = numpy.array([0.1, 0.2, 0.1, 0.3, 0.9, 0.1])
        con_i = [1,3,5]
        res = assemble._find_best_n_for_read(prob, con_i, 3)
        self.assertEqual(res, [3,1,5])

    def test_assign_reads_simple(self):
        res = assemble.assign_read_indexes(self.cons, self.em_results,
                                           self.haps, self.reads, 2.0)
        exp = {"hap1":{0, 1}, "hap2":{3}, "unassigned":{2}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_low_min_fold(self):
        res = assemble.assign_read_indexes(self.cons, self.em_results,
                                           self.haps, self.reads, 1.5)
        exp = {"hap1":{0, 1}, "hap2":{2, 3}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_high_min_fold(self):
        res = assemble.assign_read_indexes(self.cons, self.em_results,
                                           self.haps, self.reads, 200)
        exp = {"unassigned":{0, 1, 2, 3}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_only_one_con(self):
        res = assemble.assign_read_indexes(self.cons[0:1], self.em_results,
                                           self.haps, self.reads, 2)
        exp = {"hap1":{0, 1, 2, 3}}
        self.assertEqual(res, exp)

    def test_get_contrib_read_ids_simple(self):
        idxs = [0, 2]
        res = assemble.get_contrib_read_ids(idxs, self.reads)
        self.assertEqual(set('ABD'), res)

    def test_get_contrib_read_ids_all(self):
        idxs = range(len(self.reads))
        res = assemble.get_contrib_read_ids(idxs, self.reads)
        self.assertEqual(set('ABCDEFG'), res)

    def test_get_contrib_read_ids_empty(self):
        idxs = []
        res = assemble.get_contrib_read_ids(idxs, self.reads)
        self.assertEqual(set(), res)


class TestExtendAssemblies(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.ref = 'AAAAAAAAAAAAAAAAAAAA'
        self.args = parser.parse_args([])
        self.args.min_mq = 30
        self.args.min_bq = 30
        self.args.cons_cov = 2
        self.aln1 = pysam.AlignedSegment()
        self.aln1.reference_start = 10
        self.aln1.query_name = 'read1'
        self.aln1.mapping_quality = 30
        self.aln1.query_sequence = "AAAAA"
        self.aln1.query_qualities = [30] * 5
        self.aln1.cigarstring = '5M'
        self.aln2 = pysam.AlignedSegment()
        self.aln2.reference_start = 13
        self.aln2.query_name = 'read2'
        self.aln2.mapping_quality = 30
        self.aln2.query_sequence = "AAAAA"
        self.aln2.query_qualities = [30] * 5
        self.aln2.cigarstring = '5M'
        self.aln3 = pysam.AlignedSegment()
        self.aln3.reference_start = 15
        self.aln3.query_name = 'read3'
        self.aln3.mapping_quality = 30
        self.aln3.query_sequence = "TAAAA"
        self.aln3.query_qualities = [30] * 5
        self.aln3.cigarstring = '5M'
        return

    def test_call_consensus_no_aligned_sequences(self):
        res = assemble.call_consensus(self.ref, [], self.args.cons_cov,
                                      self.args, strict=True)
        exp = ''
        self.assertEqual(res, exp)

    def test_call_consensus_simple(self):
        res = assemble.call_consensus(self.ref, [self.aln1, self.aln2],
                                      self.args.cons_cov, self.args,
                                      strict=True)
        exp = 'NNNNNNNNNNNNNAANNNNN'
        self.assertEqual(res, exp)

    def test_call_consensus_high_coverage_requirement(self):
        self.args.cons_cov = 3
        res = assemble.call_consensus(self.ref, [self.aln1, self.aln2],
                                      self.args.cons_cov, self.args,
                                      strict=True)
        exp = 'NNNNNNNNNNNNNNNNNNNN'
        self.assertEqual(res, exp)

    def test_call_consensus_low_coverage_requirement(self):
        self.args.cons_cov = 1
        res = assemble.call_consensus(self.ref, [self.aln1, self.aln2],
                                      self.args.cons_cov, self.args,
                                      strict=True)
        exp = 'NNNNNNNNNNAAAAAAAANN'
        self.assertEqual(res, exp)

    def test_call_consensus_disagreement_low_coverage_requirement(self):
        self.args.cons_cov = 1
        res = assemble.call_consensus(self.ref,
                                      [self.aln1, self.aln2, self.aln3],
                                      self.args.cons_cov, self.args,
                                      strict=True)
        exp = 'NNNNNNNNNNAAAAANAAAA'
        self.assertEqual(res, exp)

    def test_call_consensus_disagreement_coverage_requirement(self):
        res = assemble.call_consensus(self.ref,
                                      [self.aln1, self.aln2, self.aln3],
                                      self.args.cons_cov, self.args,
                                      strict=True)
        exp = 'NNNNNNNNNNNNNAANAANN'
        self.assertEqual(res, exp)

    def test_find_new_variants_empty_input(self):
        res = assemble.find_new_variants(self.ref, {}, self.args)
        exp = {}
        self.assertEqual(res, exp)

    def test_find_new_variants_none_to_find(self):
        res = assemble.find_new_variants(self.ref,
                                         {'A':[self.aln1, self.aln2],
                                          'B':[self.aln1, self.aln1],
                                          'unassigned':[self.aln3]}, self.args)
        exp = {}
        self.assertEqual(res, exp)

    def test_find_new_variants_one_diff_two_contributors(self):
        res = assemble.find_new_variants(self.ref,
                                         {'A':[self.aln2, self.aln2],
                                          'B':[self.aln3, self.aln3],
                                          'unassigned':[self.aln1]}, self.args)
        exp = {(15, 'A'):'A', (15, 'T'):'B'}
        self.assertEqual(res, exp)

    def test_find_new_variants_one_diff_three_contributors(self):
        res = assemble.find_new_variants(self.ref,
                                         {'A':[self.aln2, self.aln2],
                                          'B':[self.aln3, self.aln3],
                                          'C':[self.aln2, self.aln2],
                                          'unassigned':[self.aln1]}, self.args)
        # base 'T' can be assigned to hap 'B', but 'A' cannot be assigned
        # because 'A' and 'C' both have 'A'.
        exp = {(15, 'T'):'B'}
        self.assertEqual(res, exp)

    def test_find_new_variants_two_diff_missing_cov_three_contributors(self):
        self.aln1.query_sequence = "AAGAA"
        res = assemble.find_new_variants(self.ref,
                                         {'A':[self.aln2, self.aln2],
                                          'B':[self.aln3, self.aln3],
                                          'C':[self.aln1, self.aln2,
                                               self.aln2],
                                          'unassigned':[self.aln1]}, self.args)
        # A13G variant is only covered in hap 'C'. Cannot call new variant.
        exp = {(15, 'T'):'B'}
        self.assertEqual(res, exp)

    def test_find_new_variants_one_diff_cant_assign(self):
        res = assemble.find_new_variants(self.ref,
                                         {'A':[self.aln2, self.aln3],
                                          'B':[self.aln2, self.aln3],
                                          'unassigned':[self.aln1]}, self.args)
        exp = {}
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_no_assign(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[self.aln1]}
        new_var = {(15, 'A'):'A', (15, 'T'):'B'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2],
                'B':[self.aln3, self.aln3],
                'unassigned':[self.aln1]}
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_simple_assign(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[self.aln2]}
        new_var = {(15, 'A'):'A', (15, 'T'):'B'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2, self.aln2],
                'B':[self.aln3, self.aln3],
                'unassigned':[]}
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_multi_assign(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[self.aln1, self.aln2, self.aln3]}
        new_var = {(15, 'A'):'A', (15, 'T'):'B'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2, self.aln2],
                'B':[self.aln3, self.aln3, self.aln3],
                'unassigned':[self.aln1]}
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_single_assign(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'C':[self.aln3, self.aln3],
                         'unassigned':[self.aln1, self.aln2, self.aln3]}
        new_var = {(15, 'A'):'A'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2, self.aln2],
                'B':[self.aln3, self.aln3],
                'C':[self.aln3, self.aln3],
                'unassigned':[self.aln1, self.aln3]}
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_disagree(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[self.aln1, self.aln2, self.aln3]}
        new_var = {(13, 'A'):'B', (15, 'A'):'A'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2],
                'B':[self.aln3, self.aln3, self.aln1],
                'unassigned':[self.aln2, self.aln3]}
        self.maxDiff = None
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_empty_unassigned(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[]}
        new_var = {(13, 'A'):'B', (15, 'A'):'A'}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2],
                'B':[self.aln3, self.aln3],
                'unassigned':[]}
        self.maxDiff = None
        self.assertEqual(res, exp)

    def test_assign_reads_from_new_vars_no_vars(self):
        contrib_reads = {'A':[self.aln2, self.aln2],
                         'B':[self.aln3, self.aln3],
                         'unassigned':[self.aln1, self.aln2, self.aln3]}
        new_var = {}
        res = assemble.assign_reads_from_new_vars(contrib_reads, new_var,
                                                  self.args)
        exp =  {'A':[self.aln2, self.aln2],
                'B':[self.aln3, self.aln3],
                'unassigned':[self.aln1, self.aln2, self.aln3]}
        self.maxDiff = None
        self.assertEqual(res, exp)


if __name__ == '__main__':
    unittest.main()
