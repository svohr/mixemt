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


def _run_em_once(read_hap_mat, weights, init_alpha=1.0, 
                 max_iter=1000, tolerance=0.0001, verbose=True):
    """ 
    Runs the EM algorithm to estimate haplotype contributions.
    """
    # initialize haplogroup proportions
    props = init_props(read_hap_mat.shape[1], alpha=init_alpha)

    # arrays for new calculations
    read_mix_mat = numpy.empty_like(read_hap_mat)
    new_props = numpy.empty_like(props)
   
    for iter_round in xrange(max_iter):
        if verbose and (iter_round + 1) % 10 == 0:
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
        if converged(props, new_props, tolerance):
            if verbose:
                sys.stderr.write("\nConverged! (%d)\n" % (iter_round + 1))
            break
        else:
            # Use new_prop as prop for the next iteration, use the old one
            # for empty space.
            new_props, props = props, new_props
    return new_props, read_mix_mat


def run_em(read_hap_mat, weights, n_runs=1, init_alpha=1.0, 
             max_iter=1000, tolerance=0.0001, verbose=True):
    """
    Runs EM until convergence several times (n_runs) and returns proportions
    and read/haplogroup probabilities averages over the EM runs.
    """
    if n_runs == 1:
        return _run_em_once(read_hap_mat, weights, 
                            init_alpha=init_alpha,
                            max_iter=max_iter,
                            tolerance=tolerance,
                            verbose=verbose)

    # create empty proportion vector, read matrix to keep track of average.
    avg_props = numpy.zeros(read_hap_mat.shape[1])
    avg_read_mix = numpy.zeros_like(read_hap_mat)

    for i in xrange(n_runs):
        if verbose:
            sys.stderr.write("EM iteration %d\n" % (i + 1))
        props, read_mix_mat = _run_em_once(read_hap_mat, weights, 
                                           init_alpha=init_alpha,
                                           max_iter=max_iter,
                                           tolerance=tolerance,
                                           verbose=verbose)
        avg_props += props
        avg_read_mix += read_mix_mat
    return avg_props / n_runs, avg_read_mix / n_runs
 

def main():
    """ Simple example for testing """
    if len(sys.argv) > 0:
        ref = "GAAAAAAAA"
        hap_var = dict({'A':['A2T','A4T'],
                        'B':['A3T','A5T','A6T','A8T'],
                        'C':['A3T','A6T'],
                        'D':['A3T','A5T','A7T','A9T'],
                        'E':['A3T','A4T','A5T','A7T'],
                        'F':['A3T','A5T','A6T'],
                        'G':['A3T','A5T','A7T'],
                        'H':['A3T','A5T'],
                        'I':[]})
        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDEFGHI')
        input_mat = preprocess.build_em_matrix(ref, hap_var, reads, haps)
        weights = numpy.ones(len(reads))
        props, read_mix = run_em(input_mat, weights, max_iter=100)
        print props
        print read_mix
    return 0


if __name__ == "__main__":
    sys.exit(main())

