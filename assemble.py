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
import operator
import pysam


def get_contributors(haplogroups, props, read_hap_mat, min_prob, min_reads):
    """ 
    Takes a list of haplogroup IDs, a vector of their relative contributions to
    the sample and a matrix of read haplogroup probability assignments and
    returns a list of haplogroup, proportion pairs for our putative haplogroup
    contributors. For a haplogroup to be considered as a contributor, we
    require that there must be minimum number of reads that have a minimum
    probability of originating from haplogroup.
    """
    contributors = list()
    for hap in xrange(len(haplogroups)):
        total_reads = numpy.sum(read_hap_mat[:, hap] >= min_prob)
        if total_reads >= min_reads:
            contributors.append(hap)
    contrib_prop = [[haplogroups[con], props[con]] for con in contributors]
    contrib_prop.sort(key=operator.itemgetter(1), reverse=True)

    # Add a friendly name, not associated with a haplogroup
    name_fmt = "hap%%0%dd" % (len(str(len(contrib_prop) + 1)))
    hap_num = 1
    for con in contrib_prop:
        con.insert(0, name_fmt % (hap_num))
        hap_num += 1

    return contrib_prop


def assign_reads(contribs, haps, reads, read_hap_mat, min_prob):
    """
    Takes the list of identified contributors, the list of haplotype ids and
    the read-haplogroup probability under mixture proportions matrix, and
    returns a table mapping the identified contributors to the indexes of all
    reads (rows in the matrix) that have been assigned to that haplotype and a
    set of read indexes that were not assigned to a contributor. A read is
    assigned to a haplogroup with the probability it originated from that
    haplogroup is greater or equal to the minimum cutoff provided.
    """
    contrib_reads = dict()
    unassigned = set(range(len(reads)))

    for hap, group, prop in contribs:
        hap_col = haps.index(group)
        contrib_reads[hap] = set(numpy.nonzero(read_hap_mat[:, hap_col]
                                                >= min_prob)[0])
        unassigned -= contrib_reads[hap]
    contrib_reads['unassigned'] = unassigned
    return contrib_reads


def report_contributors(out, contribs, contrib_reads, wts, header=True):
    """
    Prints a table that summarizes the contributors, their proportions and 
    number of reads assigned to each.
    """
    for hap_id, haplogroup, prop in contribs:
        total_reads = 0
        for read_index in contrib_reads[hap_id]:
            total_reads += wts[read_index]
        out.write('%s\t%s\t%.4f\t%d\n' % (hap_id, haplogroup,
                                          prop, total_reads))
    return


def get_contrib_read_ids(indexes, reads, reads_sigs):
    """
    Takes a set of indexes from assign_reads and the list of read signatures
    plus the dictionary mapping signatures to aligned read IDs and returns
    the set of corresponding aligned read IDs (SAM/BAM query IDs).
    """
    hap_read_ids = set()
    for read_idx in indexes:
        for read_id in read_sigs[reads]:
            hap_read_ids.add(read_id)
    return hap_read_ids


def write_haplotype_sam(samfile, hap_sam, hap_reads):
    """
    Takes in a SAM/BAM AlignmentFile and read haplotype assignments and
    makes new SAM/BAM files containing the subset of reads assigned to each
    haplogroup.
    """
    for aln in samfile:
        if aln in hap_reads:
            hap_sam.write(aln)
    return
