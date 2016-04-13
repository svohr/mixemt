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


def get_contributors(haplogroups, props, read_hap_mat):
    """ 
    Takes a list of haplogroup IDs, a vector of their relative contributions to
    the sample and a matrix of read haplogroup probability assignments and
    returns a list of haplogroup, proportion pairs for our putative haplogroup
    contributors. For a haplogroup to be considered as a contributor, there
    must be at least 1 read that identifies it as the most likely contributor
    according to the matrix of reads/haplogroup probabilities given mixture
    proportions.
    """
    contributors = numpy.unique(numpy.argmax(read_hap_mat))
    return [(haplogroups[con], props[con]) for con in contributors]


