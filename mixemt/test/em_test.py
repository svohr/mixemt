"""
Unit Tests for EM module.
"""

import unittest
import sys
import argparse
import math
import numpy

from mixemt import phylotree
from mixemt import preprocess
from mixemt import em


class TestEMHelpers(unittest.TestCase):
    def test_init_props(self):
        props = em.init_props(10)
        self.assertEqual(len(props), 10)
        self.assertTrue(numpy.abs(numpy.sum(props) - 1.0) < 0.000000001)

    def test_converged(self):
        prev = numpy.array([math.log(1.0)] * 10)
        cur  = numpy.array([math.log(2.0)] * 10)
        self.assertTrue(em.converged(cur, cur))
        self.assertTrue(em.converged(prev, prev))
        self.assertFalse(em.converged(prev, cur))
        self.assertFalse(em.converged(cur, prev))
        close = numpy.array(cur)
        close[3] = math.log(2.0001)
        self.assertTrue(em.converged(cur, prev, 20.0))
        self.assertFalse(em.converged(cur, close))
        self.assertTrue(em.converged(cur, close, 0.001))

    def test_em_step_simple(self):
        inf = float('Inf')
        in_mat = numpy.array([[ 0.0, -inf, -inf],
                              [-inf,  0.0, -inf],
                              [-inf, -inf,  0.0]])
        wts = numpy.array([1, 1, 1])
        props = numpy.log(numpy.array([0.6, 0.2, 0.2]))
        mix_mat = numpy.empty_like(in_mat)

        res_mat, res_props = em.em_step(in_mat, wts, props, mix_mat)
        self.assertTrue(numpy.all(res_mat == mix_mat))

        self.assertTrue(numpy.all(in_mat == mix_mat))
        self.assertTrue(numpy.all(res_props ==
                                  numpy.log(numpy.array([1.0,1.0,1.0]) / 3.0)))

    def test_em_step_weights(self):
        inf = float('Inf')
        in_mat = numpy.array([[ 0.0, -inf, -inf],
                              [-inf,  0.0, -inf],
                              [-inf, -inf,  0.0]])
        wts = numpy.array([2,1,1])
        props = numpy.log(numpy.array([0.6, 0.2, 0.2]))
        mix_mat = numpy.empty_like(in_mat)
        new_props = numpy.empty_like(props)

        res_mat, res_props = em.em_step(in_mat, wts, props, mix_mat)

        self.assertTrue(numpy.all(in_mat == mix_mat))
        self.assertTrue(numpy.all(res_props ==
                                  numpy.log(numpy.array([2.0,1.0,1.0]) / 4.0)))


class TestAllEM(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args([])
        self.args.init_alpha = 1.0
        self.args.tolerance  = 0.0001
        self.args.max_iter   = 1000
        self.args.n_multi    = 1
        self.args.verbose    = False

        phy_in = ['I, A1G ,,',
                  ',H, A3T A5T ,,',
                  ',,F, A6T ,,',
                  ',,,B, A8T ,,',
                  ',,,C, T5A ,,',
                  ',,G, A7T ,,',
                  ',,,D, A9T ,,',
                  ',,,E, A4T ,,',
                  ',A, A2T A4T ,,']
        phy = phylotree.Phylotree(phy_in)
        ref = "AAAAAAAAA"
        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDEFGHI')
        self.input_mat = preprocess.build_em_matrix(ref, phy, reads,
                                                    haps, self.args)
        self.wts = numpy.ones(len(reads))

        self.true_props = numpy.array(
                            [0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0])
        inf = float('Inf')
        self.true_haps = numpy.full_like(self.input_mat, -inf)
        self.true_haps[0:8, 1] = 0.0
        self.true_haps[8:10, 4] = 0.0

    def test_em_1_run(self):
        props, read_mix = em.run_em(self.input_mat, self.wts, self.args)
        self.assertTrue(numpy.allclose(props, self.true_props, atol=0.02))
        self.assertTrue(numpy.allclose(numpy.exp(read_mix),
                                       numpy.exp(self.true_haps), atol=0.05))

    def test_em_10_run(self):
        self.args.n_multi = 10
        true_props = numpy.array([0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0])
        props, read_mix = em.run_em(self.input_mat, self.wts, self.args)
        self.assertTrue(numpy.allclose(props, self.true_props, atol=0.02))
        self.assertTrue(numpy.allclose(numpy.exp(read_mix),
                                       numpy.exp(self.true_haps), atol=0.05))


if __name__ == '__main__':
    unittest.main()
