"""

variants.py

This file contains functions for reading a set of mapped reads in SAM/BAM 
format and discovering variant single nucleotide polymorphism (SNP) sites.
This includes positions that are polymorphic with the sample and positions
where a non-reference base is fixed. We'll assume that the reference is the
Reconstructed Sapiens Reference Sequence (RSRS).

Sam Vohr (svohr@soe.ucsc.edu)
Sun Apr  3 16:12:59 PDT 2016

"""

import pysam
import numpy

# TODO:
# Figure out if how to store this table.
# 1) Basic python list of Counters?
# 2) numpy array?


def init_obs_tab(ref_seq):
    """
    Makes an empty table for counting observed bases at reference positions.
    """
    return pos_obs


def count_read_matches(pos_obs, read, min_bq):
    """
    For all matched positions in the read, increments the count for the 
    observed base in the position table if base quality passes filter.
    """
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        qpos = int(qpos)
        rpos = int(rpos)
        if aln.query_qualities[qpos] >= min_bq:
            pos_obs[qpos][read.query_sequence[qpos]] += 1
        else:
            pos_obs[qpos]['N'] += 1
    return pos_obs


def build_obs_tab(pos_obs, reads, min_bq):
    """
    Returns a table of observation counts per base and position.
    """
    for read in reads:
        count_read_matches(pos_obs, read, min_bq)
    return 


def _find_variants(pos_obs, ref_seq):
    """ 
    Returns a list of (position, alternative base(s)) tuples that describes
    where non-reference bases are observed consistently.
    """
    return


def find_variants(reads, ref_seq, min_bq):
    """
    Takes in a list of pysam.AlignedSegments, a reference sequence and
    additional filtering parameters and returns a list of variant positions
    and alternative base(s) for positions where non-reference bases occur.
    """
    pos_obs = init_obs_tab(ref_seq)
    build_obs_tab(pos_obs, reads, min_bq)
    return _find_variants(pos_obs, ref_seq)



