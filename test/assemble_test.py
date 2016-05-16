"""
Unit Tests for assemble module.
"""

import unittest
import sys
import numpy
import argparse
import collections

import assemble
import phylotree
import preprocess

# TODO: Write tests for the follwing functions.
# _find_best_n_for_read(read_prob, con_indexes, top_n=2)
# get_contrib_read_ids(indexes, reads, read_sigs)
# assign_reads(contribs, haps, reads, read_hap_mat, props, min_fold)
# write_haplotypes(samfile, contrib_reads, reads, read_sigs, prefix, verbose)


class TestContributors(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args([])
        self.args.verbose = False
        self.args.min_reads = 1
        self.args.var_check = False

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

        self.cons = [['A', 0.4], ['E', 0.3]]
        self.obs_tab = collections.defaultdict(collections.Counter)
        self.obs_tab[1]['T'] = 1
        self.obs_tab[3]['T'] = 2
        self.obs_tab[0]['G'] = 1
        self.obs_tab[6]['T'] = 1
        self.obs_tab[2]['T'] = 1
        self.obs_tab[4]['T'] = 1

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
        res = assemble.get_contributors(self.phy, self.obs_tab, self.haps,
                                        self.em_results, self.args)
        exp = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.assertEqual(res, exp)

    def test_get_contributors_no_phy_vars_rm(self):
        self.args.var_check = False
        self.args.min_reads = 2
        res = assemble.get_contributors(self.phy, self.obs_tab, self.haps,
                                        self.em_results, self.args)
        exp = [['hap1', 'A', 0.40]]
        self.assertEqual(res, exp)

    def test_get_contributors_with_phy_vars(self):
        self.args.var_check = True 
        res = assemble.get_contributors(self.phy, self.obs_tab, self.haps,
                                        self.em_results, self.args)
        exp = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]
        self.assertEqual(res, exp)

    def test_get_contributors_with_phy_vars_rm(self):
        self.args.var_check = True
        self.args.min_reads = 2
        res = assemble.get_contributors(self.phy, self.obs_tab, self.haps,
                                        self.em_results, self.args)
        exp = []
        self.assertEqual(res, exp)

    def test_get_contributors_sorted(self):
        self.props[4] = 0.6
        self.args.var_check = False
        res = assemble.get_contributors(self.phy, self.obs_tab, self.haps,
                                        self.em_results, self.args)
        exp = [['hap1', 'E', 0.60], ['hap2', 'A', 0.4]]
        self.assertEqual(res, exp)

    def test_check_contrib_phy_vars_no_rm(self):
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertEqual(self.cons, res)

    def test_check_contrib_phy_vars_rm_none(self):
        self.cons.append(['C', 0.1])
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertNotEqual(self.cons, res)
        self.assertEqual(self.cons[0:2], res)

    def test_check_contrib_phy_vars_no_rm_one_base(self):
        self.cons.append(['C', 0.1])
        self.obs_tab[4]['A'] = 1
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertEqual(self.cons, res)

    def test_check_contrib_phy_vars_rm_shared(self):
        # Should remove E, A4T already seen in A, does not count for E.
        self.obs_tab[6]['T'] = 0
        self.obs_tab[2]['T'] = 0
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertEqual(self.cons[0:1], res)

    def test_check_contrib_phy_vars_empty_obs(self):
        # no observations, no keepers.
        self.obs_tab = collections.defaultdict(collections.Counter)
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertEqual([], res)

    def test_check_contrib_phy_vars_high_min_reads(self):
        # required number of reads too high.
        self.args.min_reads = 10
        res = assemble.check_contrib_phy_vars(self.phy, self.obs_tab,
                                              self.cons, self.args)
        self.assertEqual([], res)


class TestAssignReads(unittest.TestCase):
    def setUp(self):
        self.haps = list('ABCDEFGHI')
        self.props = numpy.array([0.40, 0.01, 0.01, 0.01, 0.3,
                                  0.01, 0.01, 0.01, 0.01])
        self.mix_mat = numpy.array([
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.91, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01],
            [0.30, 0.01, 0.01, 0.01, 0.40, 0.01, 0.01, 0.01, 0.01],
            [0.01, 0.01, 0.01, 0.01, 0.91, 0.01, 0.01, 0.01, 0.01]])
        self.reads = ['read%d' % (i + 1) for i in xrange(4)]
        self.em_results = (self.props, self.mix_mat)
        self.cons = [['hap1', 'A', 0.40], ['hap2', 'E', 0.3]]

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
        res = assemble.assign_reads(self.cons, self.haps, self.reads,
                                    self.mix_mat, self.props, 2.0)
        exp = {"hap1":{0, 1}, "hap2":{3}, "unassigned":{2}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_low_min_fold(self):
        res = assemble.assign_reads(self.cons, self.haps, self.reads,
                                    self.mix_mat, self.props, 1.5)
        exp = {"hap1":{0, 1}, "hap2":{2, 3}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_high_min_fold(self):
        res = assemble.assign_reads(self.cons, self.haps, self.reads,
                                    self.mix_mat, self.props, 200)
        exp = {"unassigned":{0, 1, 2, 3}}
        self.assertEqual(res, exp)

    def test_assign_reads_simple_only_one_con(self):
        res = assemble.assign_reads(self.cons[0:1], self.haps, self.reads,
                                    self.mix_mat, self.props, 2)
        exp = {"hap1":{0, 1, 2, 3}}
        self.assertEqual(res, exp)


if __name__ == '__main__':
    unittest.main()
