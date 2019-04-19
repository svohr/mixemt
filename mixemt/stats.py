"""
stats.py

This file contains functions for reporting results from EM and assembly steps
and raw numbers for visualization with another tool.

Sam Vohr (svohr@soe.ucsc.edu)
Fri May 20 14:04:02 PDT 2016
"""


import sys
import collections
import numpy

from mixemt import observe
from mixemt import phylotree


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
    sys.stderr.write("\nTop 10 haplogroups by read probabilities...\n")
    for hap_i, count in vote_count.most_common(top_n):
        sys.stderr.write("%s\t%d\n" % (haplogroups[hap_i], count))
    sys.stderr.write('\n')
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
    for hap_id, haplogroup, prop  in contribs:
        total_reads = len(contrib_reads[hap_id])
        if out.isatty():
            prop_str = '%.4f' % (prop)
            read_str = '%d' % (total_reads)
            out.write('%s %s %s %s\n' % (hap_id.ljust(6),
                                         haplogroup.ljust(15),
                                         prop_str.rjust(12),
                                         read_str.rjust(7)))
        else:
            out.write('%s\t%s\t%.4f\t%d\n'
                      % (hap_id, haplogroup, prop, total_reads))
    return


def write_base_obs(out, obs_tab, ref, prefix=''):
    """
    Write the counts of observed bases for each position to the file handle.

    Args:
        out: File handle to write output to.
        obs_tab: ObservedBases object of observations per reference position.
        ref: The reference sequence. Used to finding the number of positions
             we must write.
        prefix: string to write before each entry (i.e. an ID followed by a
                tab character)
    Returns: nothing
    """
    if prefix:
        prefix += '\t'
    for ref_pos in xrange(len(ref)):
        out.write(
            "%s%d\t%s\t%d\n" % (prefix, ref_pos,
                                '\t'.join([str(obs_tab.obs_at(ref_pos, base))
                                           for base in 'ACGT']),
                                sum(obs_tab.obs_tab[ref_pos].values())))
    return


def write_variants(out, phylo, contribs, obs_tab, args):
    """
    Write a table of the variants used in this analysis and note whether the
    position is expected to be polymorphic in the sample given the set of
    identified contributors.

    Args:
        out: File handle to write output to.
        phylo: The Phylotree object used in EM analysis
        contribs: Table of identified contributors with fields  hap#,
                  haplogroup, fraction
        args: The argparse namespace
    Returns: nothing
    """
    haplogroups = [con[1] for con in contribs]
    variants = collections.defaultdict(list)
    for hap in haplogroups:
        for var in phylo.hap_var[hap]:
            pos = phylotree.pos_from_var(var)
            variants[pos].append("%s:%s" % (hap, var))

    polymorphic = set(phylo.polymorphic_sites(haplogroups))
    for ref_pos in xrange(len(phylo.refseq)):
        obs = obs_tab.obs_at(ref_pos)

        samp_status = "sample_fixed"
        threshold = max(args.min_var_reads,
                        obs_tab.total_obs(pos) * args.perc_var_reads / 100)
        if sum(obs[base] >= threshold for base in 'ACGT') > 1:
            samp_status = "variant"

        phy_status = "fixed"
        if ref_pos in polymorphic:
            phy_status = "polymorphic"

        out.write("%d\t%s\t%s\t%s\t%s\n" % (ref_pos + 1,
                                            '\t'.join([str(obs[base])
                                                       for base in 'ACGT']),
                                            phy_status, samp_status,
                                            ','.join(variants[ref_pos])))
    return


def write_statistics(phylo, all_obs, contribs, contrib_reads, args):
    """
    Write a bunch of files to use for plotting the results of our EM and
    assembly steps. These will include 1) base observations for each
    contributor and 2) sites from phylotree that were used to estimate mixture
    contributions and whether or not we think these should be polymorphic
    or not.

    Args:
        phylo: The phylotree object these assignments are based on.
        ref: The reference sequence.
        all_obs: ObservedBases object of observations per reference position.
        contribs: The contributor table returned by assembly.get_contributors,
                  a list of (hap#, haplogroup, proportion) tuples.
        contrib_reads: a dictionary mapping hap#s to list of pysam
                       AlignedSegments
        args: The argparse namespace, used for the stats_prefix filename prefix
    Returns: nothing
    """
    haplogroups = {con[0]:con[1] for con in contribs}
    with open("%s.pos.tab" % (args.stats_prefix), 'w') as var_out:
        write_variants(var_out, phylo, contribs, all_obs, args)
    with open("%s.obs.tab" % (args.stats_prefix), 'w') as obs_out:
        for con in sorted(contrib_reads):
            obs_tab = observe.ObservedBases(contrib_reads[con],
                                            args.min_mq, args.min_bq)
            haplogroup = "unassigned"
            if con in haplogroups:
                haplogroup = haplogroups[con]
            write_base_obs(obs_out, obs_tab, phylo.refseq,
                           "%s\t%s" % (con, haplogroup))
        if len(contrib_reads) > 1:
            write_base_obs(obs_out, all_obs, phylo.refseq, "all\tmix")
    return
