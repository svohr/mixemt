"""

em.py

This file contains functions that implement an Expectation-Maximization
approach to identify contributing haplotypes and estimate the mixture
proportions.

Sam Vohr (svohr@soe.ucsc.edu)
Mon Apr  4 09:38:08 PDT 2016

"""

import sys
import numpy
import pysam

import phylotree
import preprocess


def init_props(nhaps, alpha=1.0):
    """
    Initializes haplogroup proportions by random draw from Dirichlet with
    a uniform alpha parameter.
    """
    return numpy.random.dirichlet([alpha] * nhaps)


def converged(prop, last_prop, tolerance=0.0001):
    """
    Check if our proportion parameters has converged, i.e. the absolute
    difference between the this and the previous iteration is within a 
    fixed tolerance value. 
    """
    return numpy.sum(numpy.abs(prop - last_prop)) < tolerance


def run_em(read_hap_mat, weights, max_iter=1000):
    """ 
    Runs the EM algorithm to estimate haplotype contributions.
    """
    # initialize haplogroup proportions
    props = init_props(read_hap_mat.shape[1], alpha=1)
    
    # arrays for new calculations
    read_mix_mat = numpy.empty_like(read_hap_mat)
    new_props = numpy.empty_like(props)
   
    for iter_round in xrange(max_iter):
        if (iter_round + 1) % 10 == 0:
            sys.stderr.write('.')
        # Set z_j,g - probablilty that read j originates from haplogroup g 
        # given this proportion in the mixture..
        for i in xrange(read_hap_mat.shape[0]):
            prop_read = props * read_hap_mat[i, ]
            read_mix_mat[i, ] = prop_read / numpy.sum(prop_read)
        
        # Set theta_g - contribution of g to the mixture
        for i in xrange(read_hap_mat.shape[1]):
            new_props[i] = numpy.sum(read_mix_mat[:, i] * weights)
        new_props /= numpy.sum(new_props)

        # check for convergence.
        if converged(props, new_props):
            print "\nConverged!"
            break
        else:
            # Use new_prop as prop for the next iteration, use the old one
            # for empty space.
            new_props, props = props, new_props
    return new_props, read_mix_mat


def main():
    """ Simple example for testing """
    if len(sys.argv) > 2:
        phy_fn = sys.argv[1]
        bam_fn = sys.argv[2]
        ref_in = pysam.FastaFile('../ref/RSRS.mtDNA.fa')
        refseq = ref_in.fetch(ref_in.references[0])
        with open(phy_fn, 'r') as phy_in:
            var_pos, hap_var = phylotree.read_phylotree(phy_in,
                                                        False, False, False)
        with pysam.AlignmentFile(bam_fn, 'rb') as samfile:
            input_mat, weights, haps, reads, read_sigs = \
                preprocess.build_em_input(samfile, refseq, var_pos, hap_var)
        print run_em(input_mat, weights, max_iter=10000)
    else:
        ref = "GAAAAAAAA"
        var_pos = range(1, 9)
        hap_var = dict({'A':['A1T','A3T'],
                        'B':['A2T','A4T','A5T','A7T'],
                        'C':['A2T','A5T'],
                        'D':['A2T','A4T','A6T','A8T'],
                        'E':['A2T','A3T','A4T','A6T']})
        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDE')
        input_mat = preprocess.build_em_matrix(ref, hap_var, reads, haps)
        weights = numpy.ones(len(reads))
        props, read_mix = run_em(input_mat, weights, max_iter=100)
        print props
        print read_mix
    return 0


if __name__ == "__main__":
    sys.exit(main())

