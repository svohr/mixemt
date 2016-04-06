"""

em.py

This file contains functions that implement an Expectation-Maximization
approach to identify contributing haplotypes and estimate the mixture
proportions.

Sam Vohr (svohr@soe.ucsc.edu)
Mon Apr  4 09:38:08 PDT 2016

"""

import numpy

import phylotree


def init_props(haplogroups, alpha=1.0):
    """
    Initializes haplogroup proportions by random draw from Dirichlet with
    a uniform alpha parameter.
    """
    return numpy.random.dirichlet([alpha] * len(haplogroups))


def converged(prop, last_prop, tolerance=0.0001):
    """
    Check if our proportion parameters has converged, i.e. the absolute
    difference between the this and the previous iteration is within a 
    fixed tolerance value. 
    """
    return numpy.sum(numpy.abs(prop - last_prop)) < tolerance


def run_em(haplogroups, max_iter=10000):
    """ 
    Runs the EM algorithm to estimate haplotype contributions.
    """
    # initialize haplogroup proportions
    props = init_props(haplogroups)
    for i in xrange(max_iter):
        # Set z_j,g - probablilty that read j originates from haplogroup g.
        # Set theta_g - contribution of g to the mixture
        # check for convergence.
        if converged(prop, last_prop):
            break
    return props
