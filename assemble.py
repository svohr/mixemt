"""

assemble.py

This module contains functions for interpretting the output from our EM
algorithm: the vector of haplogroup contributions and the matrix of
read-haplogroup assignments. From these two data sets, we want to extract the
number of contributors and the fraction of the sample that is attributible to
each and individual assemblies for each contributor.

Sam Vohr (svohr@soe.ucsc.edu)

Wed Apr 13 10:57:39 PDT 2016

"""

import numpy


def get_contributors(haplogroups, props, read_hap_mat, min_prob):
    """ 
    Takes a list of haplogroup IDs, a vector of their relative contributions to
    the sample and a matrix of read haplogroup probability assignments and
    returns a list of haplogroup, proportion pairs for our putative haplogroup
    contributors. For a haplogroup to be considered as a contributor, there
    must be at least 1 read that identifies it as the most likely contributor
    according to the matrix of reads/haplogroup probabilities given mixture
    proportions.
    """
    contributors = numpy.unique(numpy.argmax(read_hap_mat, 1))
    return [(haplogroups[con], props[con]) for con in contributors]


def assign_reads(contribs, haps, read_hap_mat, min_prob):
    """
    Takes the list of identified contributors, the list of haplotype ids and
    the read-haplogroup probability under mixture proportions matrix, and
    returns a table mapping the identified contributors to the indexes of all
    reads (rows in the matrix) that have been assigned to that haplotype. A
    read is assigned to a haplogroup with the probability it originated from
    that haplogroup is greater or equal to the minimum cutoff provided.
    """
    contrib_reads = dict()
    for hap, prop in contribs:
        hap_col = haps.index(hap)
        contrib_reads[hap] = numpy.nonzero(read_hap_mat[:, hap_col] 
                                             >= min_prob)
    return contrib_reads

