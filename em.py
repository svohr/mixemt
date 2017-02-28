"""

em.py

This file contains functions that implement an Expectation-Maximization
approach to identify contributing haplotypes and estimate the mixture
proportions.

Sam Vohr (svohr@soe.ucsc.edu)
Mon Apr  4 09:38:08 PDT 2016

"""

import sys
import argparse
import numpy
import scipy.misc

import preprocess
import phylotree


def init_props(nhaps, alpha=1.0):
    """
    Initializes haplogroup proportions by random draw from Dirichlet with
    a uniform alpha parameter.

    Args:
        nhaps: number of haplogroups and length of the array to produce
        alpha: alpha parameter for dirichlet distribution
    Returns:
        a numpy array of random values that sum to 1.0
    """
    return numpy.random.dirichlet([alpha] * nhaps)


def converged(prop, last_prop, tolerance=0.0001):
    """
    Check if our proportion parameters has converged, i.e. the absolute
    difference between the this and the previous iteration is within a
    fixed tolerance value.

    Args:
        prop, last_prop: two proportion results from em steps
        tolerance: an absolute tolerance value to use when checking for
                   convergence.
    Returns:
        True if sum of absolute differences between proportions is within
        tolerance value, False otherwise.
    """
    return numpy.sum(numpy.abs(numpy.exp(prop)
                               - numpy.exp(last_prop))) < tolerance


def em_step(read_hap_mat, weights, ln_props, read_mix_mat):
    """
    Runs a single iteration of the EM algorithm given:
    1. The original input read-hap probability matrix (read_hap_mat)
    2. The weights, i.e. the number times we observe each sub-haplotype.
    3. Current haplogroup proportion estimates (props)
    4. A matrix of the same size as 'read_hap_mat' that will be written over
       with the conditional probabililies (read_mix_mat)
    5. A vector of the same size as 'props' that will be filled by the new
       haplogroup proportion estimates.

    Args:
        read_hap_mat: starting read/haplogroup conditional prob. matrix
        weights: arrray containing the number of times each sub-haplotype
                 was found in the reads (same length as read_map_hat height)
        ln_props: starting mixture proportion array, log transformed
        read_mix_mat: destination for updated conditional probability matrix
    Returns:
        Updated conditional probability matrix and mixture proportions array
    """
    # E-Step:
    # Set z_j,g - probablilty that read j originates from haplogroup g
    # given this proportion in the mixture.
    numpy.add(ln_props, read_hap_mat, read_mix_mat)
    numpy.subtract(read_mix_mat,
                   scipy.misc.logsumexp(read_mix_mat, axis=1).reshape((-1, 1)),
                   read_mix_mat)

    # M-Step:
    # Set theta_g - contribution of g to the mixture
    new_props = scipy.misc.logsumexp(read_mix_mat, axis=0,
                                     b=weights.reshape((-1, 1)))
    new_props -= scipy.misc.logsumexp(new_props)

    return read_mix_mat, new_props


def run_em(read_hap_mat, weights, args):
    """
    Runs the EM algorithm on the input read-haplogroup probability matrix and
    weights, using the arguments read from the command line.

    Args:
        read_hap_mat: input matrix containing probability of each haplogroup
                      producing each unique fragment.
        weights: arrray containing the number of times each sub-haplotype
                 was found in the reads (same length as read_map_hat height)
        args: arguments namespace from argparse.
    Returns:
        res_props: final estimated mixture proportions (linear)
        res_read_mix: final estimated read/haplogroup conditional probability
                      matrix (log)
    """
    # arrays for new calculations
    read_mix_mat = numpy.empty_like(read_hap_mat)
    new_props = numpy.empty(read_hap_mat.shape[1])

    # results for multiple runs if necessary.
    res_props, res_read_mix = None, None

    for i in xrange(args.n_multi):

        if args.verbose:
            sys.stderr.write("Starting EM run %d...\n" % (i + 1))

        # initialize haplogroup proportions
        props = numpy.log(init_props(read_hap_mat.shape[1],
                                     alpha=args.init_alpha))

        for iter_round in xrange(args.max_iter):
            if args.verbose and (iter_round + 1) % 10 == 0:
                sys.stderr.write('.')
            # Run a single step of EM
            read_mix_mat, new_props = em_step(read_hap_mat, weights,
                                              props, read_mix_mat)
            # Check for convergence.
            if converged(props, new_props, args.tolerance):
                if args.verbose:
                    sys.stderr.write("\nConverged! (%d)\n" % (iter_round + 1))
                break
            else:
                # Use new_prop as prop for the next iteration, use the old one
                # for empty space.
                new_props, props = props, new_props
        else:
            # Never converged, flip props and new_props back
            new_props, props = props, new_props

        if res_props is None and res_read_mix is None:
            # First iteration, store the results
            res_props = new_props
            res_read_mix = read_mix_mat
            if args.n_multi > 1:
                # If we are doing more than 1 iteration, make a new matrix
                read_mix_mat = numpy.empty_like(read_hap_mat)
                new_props = props
        else:
            # Nth iteration, add the new results to the running total.
            res_props += new_props
            numpy.logaddexp(res_read_mix, read_mix_mat, res_read_mix)

    if args.n_multi > 1:
        # Take the average if this is a multi-run.
        res_props /= args.n_multi
        res_read_mix -= numpy.log(args.n_multi)

    res_props = numpy.exp(res_props)

    return res_props, res_read_mix


def main():
    """ Simple example for testing """
    if len(sys.argv) > 0:
        parser = argparse.ArgumentParser()
        args = parser.parse_args()
        args.init_alpha = 1.0
        args.tolerance = 0.0001
        args.max_iter = 1000
        args.n_multi = 1
        args.verbose = True

        phy = phylotree.example()

        ref = "AAAAAAAAA"

        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDEFGHI')
        input_mat = preprocess.build_em_matrix(ref, phy, reads, haps, args)
        weights = numpy.ones(len(reads))
        props, read_mix = run_em(input_mat, weights, args)
        print props
        print read_mix
    return 0


if __name__ == "__main__":
    sys.exit(main())
