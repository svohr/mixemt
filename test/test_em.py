"""
Unit Tests for EM module.
"""

import unittest
import sys
import numpy
import argparse

import phylotree
import preprocess
import em 


class TestEMHelpers(unittest.TestCase):
    def test_init_props(self):
        props = em.init_props(10)
        self.assertEqual(len(props), 10)
        self.assertTrue(numpy.abs(numpy.sum(props) - 1.0) < 0.000000001)
    
    def test_converged(self):
        prev = numpy.array([1.0] * 10)
        cur  = numpy.array([2.0] * 10)
        self.assertTrue(em.converged(cur, cur))
        self.assertTrue(em.converged(prev, prev))
        self.assertFalse(em.converged(prev, cur))
        self.assertFalse(em.converged(cur, prev))
        close = numpy.array(cur)
        close[3] = 2.0001
        self.assertTrue(em.converged(cur, prev, 20.0))
        self.assertFalse(em.converged(cur, close))
        self.assertTrue(em.converged(cur, close, 0.001))

    def test_em_step_simple(self):
        in_mat = numpy.array([[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]])
        wts = numpy.array([1,1,1])
        props = numpy.array([0.6, 0.2, 0.2])
        mix_mat = numpy.empty_like(in_mat)
        new_props = numpy.empty_like(props)

        res_mat, res_props = em.em_step(in_mat, wts, props, mix_mat, new_props)
        self.assertTrue(numpy.all(res_mat == mix_mat))
        self.assertTrue(numpy.all(res_props == new_props))
        
        self.assertTrue(numpy.all(in_mat == mix_mat))
        self.assertTrue(
            numpy.all(new_props == numpy.array([1.0,1.0,1.0]) / 3.0))
    
    def test_em_step_weights(self):
        in_mat = numpy.array([[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]])
        wts = numpy.array([2,1,1])
        props = numpy.array([0.6, 0.2, 0.2])
        mix_mat = numpy.empty_like(in_mat)
        new_props = numpy.empty_like(props)

        res_mat, res_props = em.em_step(in_mat, wts, props, mix_mat, new_props)

        self.assertTrue(numpy.all(in_mat == mix_mat))
        self.assertTrue(
            numpy.all(new_props == numpy.array([2.0,1.0,1.0]) / 4.0))


class TestAllEM(unittest.TestCase):
    def setUp(self):
        parser = argparse.ArgumentParser()
        self.args = parser.parse_args()
        self.args.init_alpha = 1.0
        self.args.tolerance  = 0.0001
        self.args.max_iter   = 1000
        self.args.n_multi    = 1
        self.args.verbose    = False

        phy = phylotree.example()
        ref = "AAAAAAAAA"
        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDEFGHI')
        self.input_mat = preprocess.build_em_matrix(ref, phy, reads, haps)
        self.wts = numpy.ones(len(reads))
        
        self.true_props = numpy.array(
                            [0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0])
        self.true_haps = numpy.zeros_like(self.input_mat)
        self.true_haps[0:8,1] = 1.0
        self.true_haps[8:10,4] = 1.0
    
    def test_em_1_run(self):
        props, read_mix = em.run_em(self.input_mat, self.wts, self.args)
        self.assertTrue(numpy.all(numpy.abs(props - self.true_props) < 0.02))
        self.assertTrue(numpy.all(numpy.abs(read_mix - self.true_haps) < 0.05))

    def test_em_10_run(self):
        self.args.n_multi = 10
        true_props = numpy.array([0.0, 0.8, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0, 0.0])
        props, read_mix = em.run_em(self.input_mat, self.wts, self.args)
        self.assertTrue(numpy.all(numpy.abs(props - self.true_props) < 0.02))
        self.assertTrue(numpy.all(numpy.abs(read_mix - self.true_haps) < 0.05))

    

if __name__ == '__main__':
    unittest.main()
