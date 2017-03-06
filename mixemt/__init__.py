"""

mixemt

This package contains modules used by the mixemt binary to read in mapped
sequences, store the haplogroup tree, perform mixture analysis using
expectation maximization (EM) and report the results.

modules:
    phylotree - Class and functions for Phylotree haplogroups and variants
    preprocess - Functions for building EM input from mapped reads
    em - Functions for running the expectation maximization algorithm
    observe - Class for storing the bases observed at each reference position.
    assemble - Functions for interpretting EM results.
    stats - Functions for summarizing and writing out results.

"""

from mixemt import assemble
from mixemt import em
from mixemt import observe
from mixemt import phylotree
from mixemt import preprocess
from mixemt import stats
