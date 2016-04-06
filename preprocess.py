"""

preprocess.py

This directory contains functions and classes for preprocessing raw data into
matrices and vectors we can run our EM algorithm on. A big part of this is to
reduce the overall size of our input matrixes to make our algorithm run quickly
and without using too much memory.

Sam Vohr (svohr@soe.ucsc.edu)
Wed Apr  6 10:48:21 PDT 2016

"""

import numpy
import pysam
import collections
import bisect

import phylotree


class HapVarBaseMatrix(object):
    """
    This class is designed to simplify a pretty large 3 dimensional matrix of
    indicator variables. Since our matrix is sparse and only need to look up
    single cells, we can store this information in a simplied way. The first
    part is that there are only two possible values that we store, a value if
    the haplogroup, variant site, observed base is expected (1) or unexpected
    (0). In practice, we will use fractional values (e.g. 0.991 and 0.003) to
    allow for sequencing error or de novo mutations.  Second, we keep a
    dictionary of haplogroup, variant site, base tuples that mark where a
    non-reference base occurs. If an entry is not present in the dictionary, it
    must be the case that the reference base is the "expected" base and all
    others are "unexpected".
    """
    def __init__(self, refseq, hap_var=None, exp=0.991, unexp=0.003):
        """ Initialzes the HapVarBaseMatrix. """
        self.refseq = refseq
        self.exp = exp 
        self.unexp = unexp
        self.markers = set()
        
        if hap_var is not None:
            add_hap_markers(hap_var)
        return

    def add_hap_markers(hap_var):
        """ 
        Populates the table with bases that are associated with haplogroups.
        """
        for hap in hap_var:
            for var in hap_var[hap]:
                pos = phylotree.pos_from_var(var)
                der = phylotree.der_allele(var)
                if der != self.refseq[pos]:
                    self.markers.add((hap, pos, der))
        return self

    def prob(self, hap, pos, base):
        """ 
        Returns the "probability" of observing this base, at this position
        in this haplogroup.
        """
        if (hap, pos, base) in self.markers:
            return self.exp
        else:
            if self.refseq[pos] == base:
                return self.exp
        return self.unexp 


def process_reads(samfile, var_pos, min_mq, min_bq):
    """
    Reads in a set of reads from a SAM/BAM file and makes observations for
    known variant positions. Each read is simplified into a signature of
    observed base per variant site.
    """
    read_sigs = collections.defaultdict(list)
    for aln in samfile:
        if aln.mapping_quality >= min_mq:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                qpos = int(qpos)
                rpos = int(rpos)
                if qpos in var_pos:
                    if aln.query_qualities[qpos] >= min_bq:
                        # Add this to the list
                        aln.query_sequence[qpos].upper()
    return read_sigs


def build_em_input(refseq, var_pos, hap_tab, aln_fn):
    """
    Builds the matrix that describes the "probability" that a read originated
    from a specific haplogroup.
    """
    return
