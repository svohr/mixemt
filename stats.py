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

import preprocess


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


def write_base_obs(out, obs_tab, prefix=''):
    """
    Write the counts of observed bases for each position to the file handle.

    Args:
        out: File handle to write output to.
        obs_tab: Table of base observations for positions in the reference.
    Returns: nothing
    """
    if prefix:
        prefix += '\t'
    for ref_pos in xrange(max(obs_tab)):
        out.write("%s%d\t%s\t%d\n" % (prefix, ref_pos,
                                      '\t'.join([str(obs_tab[ref_pos][base])
                                       for base in 'ACGT']),
                                      sum(obs_tab[ref_pos].values())))
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
    haplogroups = [con[1] for con in contribs]
    polymorphic = set(phylo.polymorphic_sites(haplogroups, ref))
    for var in phylo.get_variant_pos():
        status = "fixed"
        if var in polymorphic:
            status = "polymorphic"
        out.write("%d\t%s\n" % (var, status))
    return


def write_statistics(phylo, ref, contribs, contrib_reads, args):
    """
    Write a bunch of files to use for plotting the results of our EM and
    assembly steps. These will include 1) base observations for each
    contributor and 2) sites from phylotree that were used to estimate mixture
    contributions and whether or not we think these should be polymorphic
    or not.

    Args:
        contribs: The contributor table returned by assembly.get_contributors,
                  a list of (hap#, haplogroup, proportion) tuples.
        contrib_reads: a dictionary mapping hap#s to list of pysam
                       AlignedSegments
        obs_tab: The table of single base per reference position for the
                 entire sample.
        ref: The reference sequence.
        args: The argparse namespace, used for the stats_prefix filename prefix
    Returns: nothing
    """
    haplogroups = {con[0]:con[1] for con in contribs}
    with open("%s.pos.tab" % (args.stats_prefix), 'w') as var_out:
        write_variants(var_out, phylo, ref, contribs)
    for con in contrib_reads:
        with open("%s.%s.obs.tab" % (args.stats_prefix, con), 'w') as obs_out:
            _, obs_tab = preprocess.process_reads(contrib_reads[con], [],
                                                  args.min_mq, args.min_bq)
            haplogroup = ""
            if con in haplogroups:
                haplogroup = haplogroups[con]
            write_base_obs(obs_out, obs_tab, "%s\t%s" % (con, haplogroup))
    return

