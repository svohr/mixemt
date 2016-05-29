"""
stats.py

This file contains functions for reporting results from EM and assembly steps
and raw numbers for visualization with another tool.

Sam Vohr (svohr@soe.ucsc.edu)
Fri May 20 14:04:02 PDT 2016
"""


import sys
import numpy
import collections

import phylotree

def report_top_props(haplogroups, props, top_n=10):
    """
    Prints to stderr the names and fractions of the n haplogroups with the
    highest estimated proportions.
    """
    order = numpy.argsort(props)[::-1]
    sys.stderr.write('\nTop %d haplogroups by proportion...\n' % (top_n))
    for i in xrange(top_n):
        sys.stderr.write("%d\t%0.6f\t%s\n" % (i + 1, props[order[i]],
                                              haplogroups[order[i]]))
    sys.stderr.write('\n')
    return


def report_read_votes(haplogroups, read_hap_mat, top_n=10):
    """
    Each read "votes" for a haplogroup; the haplogroup with the highest
    probability. Report the vote counts for the top N.
    """
    votes = numpy.argmax(read_hap_mat, 1)
    vote_count = collections.Counter(votes)
    for hap_i, count in vote_count.most_common(top_n):
        sys.stderr.write("%s\t%d\n" % (haplogroups[hap_i], count))
    return


def report_contributors(out, contribs, contrib_reads):
    """
    Prints a table that summarizes the contributors, their proportions and
    number of reads assigned to each. Formats the output nicely if out is
    a TTY, otherwise prints a tab-delimited table.
    """
    if out.isatty():
        out.write("hap#   Haplogroup      Contribution   Reads\n")
        out.write("-------------------------------------------\n")
    for hap_id, haplogroup, prop in contribs:
        total_reads = len(contrib_reads[hap_id])
        if out.isatty():
            prop_str = '%.4f' % (prop)
            read_str = '%d' % (total_reads)
            out.write('%s %s %s %s\n' % (hap_id.ljust(6),
                                         haplogroup.ljust(15),
                                         prop_str.rjust(12),
                                         read_str.rjust(7)))
        else:
            out.write('%s\t%s\t%.4f\t%d\n' % (hap_id, haplogroup,
                                              prop, total_reads))
    return


def write_coverage(out, obs_tab):
    """
    Write the coverage at each position to the file handle.

    Args:
        out: File handle to write output to.
        obs_tab: Table of base observations for positions in the reference.
    Returns: nothing
    """
    for ref_pos in xrange(max(obs_tab)):
        cov = sum(obs_tab[ref_pos])
        out.write('%d\t%d\n' % (ref_pos, cov))
    return


def write_base_obs(out, obs_tab):
    """
    Write the counts of observed bases for each position to the file handle.

    Args:
        out: File handle to write output to.
        obs_tab: Table of base observations for positions in the reference.
    Returns: nothing
    """
    for ref_pos in xrange(max(obs_tab)):
        out.write("%d\t%s\n" % (ref_pos, ['\t'.join(obs_tab[ref_pos][base])
                                          for base in 'ACGT']))
    return


def write_variants(out, phylo, ref, contribs):
    """
    Write a table of the variants used in this analysis and note whether the
    position is expected to be polymorphic in the sample given the set of
    identified contributors.

    Args:
        out: File handle to write output to.
        phylo: The Phylotree object used in EM analysis
        contribs: Table of identified contributors with fields  hap#,
                  haplogroup, fraction
    Returns: nothing
    """
    for var in phylo.get_variant_pos():
        out.write("%d\n" % (var))
    return

