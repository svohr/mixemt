"""
Unit Tests for observe module.
"""

import unittest
import sys
import pysam

from mixemt import observe


class TestObservedBases(unittest.TestCase):
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

    def test_observed_bases_init_empty(self):
        res = observe.ObservedBases()
        self.assertEqual(res.obs_tab, {})

    def test_observed_bases_init_from_alns(self):
        res = observe.ObservedBases(self.alns, 20, 10)
        exp = {10:{'A':1},
               11:{'A':1},
               12:{'A':2},
               13:{'A':2},
               14:{'A':2},
               15:{'G':1, 'T':1},
               16:{'A':2},
               17:{'A':1, '-':1},
               18:{'A':1, '-':1},
               19:{'A':2},
               20:{'G':1, 'T':1},
               21:{'A':2},
               22:{'A':2},
               23:{'A':2},
               24:{'A':2},
               25:{'G':1, 'T':1}}
        self.assertEqual(res.obs_tab, exp)

    def test_observed_bases_update_after_init(self):
        res = observe.ObservedBases(self.alns, 20, 10)
        res.update(self.alns)
        exp = {10:{'A':2},
               11:{'A':2},
               12:{'A':4},
               13:{'A':4},
               14:{'A':4},
               15:{'G':2, 'T':2},
               16:{'A':4},
               17:{'A':2, '-':2},
               18:{'A':2, '-':2},
               19:{'A':4},
               20:{'G':2, 'T':2},
               21:{'A':4},
               22:{'A':4},
               23:{'A':4},
               24:{'A':4},
               25:{'G':2, 'T':2}}
        self.assertEqual(res.obs_tab, exp)

    def test_observed_bases_init_reverse_strand(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        exp = {10:{'A':1},
               11:{'A':1},
               12:{'A':1, 'a':1},
               13:{'A':1, 'a':1},
               14:{'A':1, 'a':1},
               15:{'g':1, 'T':1},
               16:{'A':1, 'a':1},
               17:{'A':1, '+':1},
               18:{'A':1, '+':1},
               19:{'A':1, 'a':1},
               20:{'g':1, 'T':1},
               21:{'A':1, 'a':1},
               22:{'A':1, 'a':1},
               23:{'A':1, 'a':1},
               24:{'A':1, 'a':1},
               25:{'g':1, 'T':1}}
        self.assertEqual(res.obs_tab, exp)

    def test_observed_bases_obs_at_basic_no_base(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        self.assertEqual(res.obs_at(10), {'A':1})
        self.assertEqual(res.obs_at(14), {'A':2})
        self.assertEqual(res.obs_at(17), {'A':1, '-':1})

    def test_observed_bases_obs_at_basic_no_base_stranded(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        self.assertEqual(res.obs_at(10, stranded=True), {'A':1})
        self.assertEqual(res.obs_at(14, stranded=True), {'A':1, 'a':1})
        self.assertEqual(res.obs_at(17, stranded=True), {'A':1, '+':1})

    def test_observed_bases_obs_at_with_base(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        self.assertEqual(res.obs_at(14, 'A'), 2)
        self.assertEqual(res.obs_at(14, 'G'), 0)
        self.assertEqual(res.obs_at(15, 'G'), 1)
        self.assertEqual(res.obs_at(17, '-'), 1)
        self.assertEqual(res.obs_at(17, '+'), 1)

    def test_observed_bases_obs_at_with_base_symmetry(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        for base in 'ACGTN':
            self.assertEqual(res.obs_at(15, base),
                             res.obs_at(15, base.lower()))
        self.assertEqual(res.obs_at(17, '-'), res.obs_at(17, '+'))

    def test_observed_bases_obs_at_with_base_stranded(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        self.assertEqual(res.obs_at(14, 'A', True), (1, 1))
        self.assertEqual(res.obs_at(14, 'G', True), (0, 0))
        self.assertEqual(res.obs_at(15, 'G', True), (0, 1))
        self.assertEqual(res.obs_at(17, '-', True), (0, 1))
        self.assertEqual(res.obs_at(17, '+', True), (0, 1))

    def test_observed_bases_obs_at_with_base_stranded_symmetry(self):
        self.alns[1].is_reverse = True
        res = observe.ObservedBases(self.alns, 20, 10)
        for base in 'ACGTN':
            self.assertEqual(res.obs_at(15, base, stranded=True),
                             res.obs_at(15, base.lower(), stranded=True))
        self.assertEqual(res.obs_at(17, '-', True), res.obs_at(17, '+', True))

    def test_observed_bases_obs_at_bad_base(self):
        res = observe.ObservedBases(self.alns, 20, 10)
        with self.assertRaises(ValueError):
            res.obs_at(15, 'Q')


if __name__ == '__main__':
    unittest.main()
