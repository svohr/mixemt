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
import pysam
import operator
import collections
import sys


def get_contributors(phylo, obs_tab, haplogroups, em_results, args):
    """
    Takes the haplogroup tree, table of observed bases by reference position,
    a list of haplogroup IDs and a tuple representing the results from EM
    (haplogroup proportions and read/haplogroup probabilities) and identifies
    contributors that pass our filtering steps. 

    For a haplogroups to be considered a contributor:
    1) There must exist k reads where that haplogroup represent the most
       probable source of that read. (Candidate contributor)
    2) A majority of the variant bases that are associated with the candidate
       contributor, but no other candidate contributor are observed in the
       sample.
    
    This function returns a list of hap numbers, haplogroup IDs and estimated
    contributions.
    """
    props, read_hap_mat = em_results

    contributors = list()
    vote_count = collections.Counter(numpy.argmax(read_hap_mat, 1))
    contributors = [con for con in vote_count 
                    if vote_count[con] > args.min_reads]
    contrib_prop = [[haplogroups[con], props[con]] for con in contributors]
    contrib_prop.sort(key=operator.itemgetter(1), reverse=True)

    # Add a friendly name, not associated with a haplogroup
    name_fmt = "hap%%0%dd" % (len(str(len(contrib_prop) + 1)))
    hap_num = 1
    for con in contrib_prop:
        con.insert(0, name_fmt % (hap_num))
        hap_num += 1

    return contrib_prop


def check_contrib_phy_vars(phylo, obs_tab, contribs, args):
    """
    Checks if each candidate contributor from contribs passes our variant base 
    check. The strategy for this is to start with the highest estimated
    contributors and an empty list of variant positions. For each contributor,
    we identify the variant bases that are unique from the previous candidates.
    We check the observation table to verify that those bases are observed in
    the sample.
    """
    return


def _find_best_n_for_read(read_prob, con_indexes, top_n=2):
    """
    Takes a vector of haplogroup probabilities for a single read (read_prob)
    and a list of indexes of identified contributors and returns a list of
    N (top_n) indexes in con_indexes that have the highest probabilities.
    """
    order = numpy.argsort(read_prob)[::-1]
    results = list()
    i = 0

    while len(results) < top_n:
        if order[i] in con_indexes:
            results.append(order[i])
        i += 1
    return results
    

def assign_reads(contribs, haps, reads, read_hap_mat, props, min_fold):
    """
    Takes the list of identified contributors, the list of haplotype ids and
    the read-haplogroup probability under mixture proportions matrix, and
    returns a table mapping the identified contributors to the indexes of all
    reads (rows in the matrix) that have been assigned to that haplotype and a
    set of read indexes that were not assigned to a contributor. A read is
    assigned to a haplogroup if its highest probability out of all identified
    contributors is at least min_fold times greater than the probability of
    the next contributor _after_ normalizing out the contribution proportion
    of each. This way, reads that could be assigned to more than 1 
    contributor are not simply assigned to the contributor that represents
    the larger proportion of the mixture.
    """
    contrib_reads = collections.defaultdict(set)
    
    index_to_hap = dict([(haps.index(group), hap_n)
                        for hap_n, group, _ in contribs])
    con_indexes = set(index_to_hap.keys()) 
    for read_i in xrange(len(reads)):
        if len(contribs) > 1:
            read_probs = read_hap_mat[read_i, ]
            best_hap, next_hap = _find_best_n_for_read(read_probs,
                                                       con_indexes,
                                                       top_n=2)
            rel_prob = ((read_probs[best_hap] / props[best_hap]) / 
                        (read_probs[next_hap] / props[next_hap]))
            if rel_prob >= min_fold:
                contrib_reads[index_to_hap[best_hap]].add(read_i)
            else:
                contrib_reads['unassigned'].add(read_i)
        else:
            # Only 1 identified contributor, use the first field of the first
            # entry in the contributors table.
            contrib_reads[contribs[0][0]].add(read_i)
    return contrib_reads


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


def report_contributors(out, contribs, contrib_reads, wts):
    """
    Prints a table that summarizes the contributors, their proportions and 
    number of reads assigned to each. Formats the output nicely if out is
    a TTY, otherwise prints a tab-delimited table.
    """
    if out.isatty():
        out.write("hap#   Haplogroup      Contribution   Reads\n")
        out.write("-------------------------------------------\n")
    for hap_id, haplogroup, prop in contribs:
        total_reads = 0
        for read_index in contrib_reads[hap_id]:
            total_reads += wts[read_index]
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


def get_contrib_read_ids(indexes, reads, read_sigs):
    """
    Takes a set of indexes from assign_reads and the list of read signatures
    plus the dictionary mapping signatures to aligned read IDs and returns
    the set of corresponding aligned read IDs (SAM/BAM query IDs).
    """
    hap_read_ids = set()
    for read_idx in indexes:
        for read_id in read_sigs[reads[read_idx]]:
            hap_read_ids.add(read_id)
    return hap_read_ids


def write_haplotypes(samfile, contrib_reads, reads, read_sigs, prefix, verbose):
    """
    For each contributing haplotype in contrib_reads, creates a new SAM/BAM
    file based on the one provided and the prefix string and writes the reads
    from the original SAM/BAM input that have been assigned to the contributor.
    A separate file is also written for unassigned reads.
    """
    ext = samfile.filename[-3:]
    mode = 'wb'
    if ext != 'bam':
        mode = 'w'
    if verbose:
        sys.stderr.write('\n')
    for contrib in contrib_reads:
        if len(contrib_reads[contrib]) == 0:
            continue
        hap_read_ids = get_contrib_read_ids(contrib_reads[contrib], 
                                            reads, read_sigs)
        hap_fn = "%s.%s.%s" % (prefix, contrib, ext)
        try:
            hap_samfile = pysam.AlignmentFile(hap_fn, mode, template=samfile)
            written = 0
            for aln in samfile.fetch():
                if aln.query_name in hap_read_ids:
                    hap_samfile.write(aln)
                    written += 1
            if verbose:
                sys.stderr.write('Wrote %d aligned segments to %s\n' 
                                 % (written, hap_fn))
            hap_samfile.close()
        except (ValueError, IOError) as inst:
            sys.stderr.write("Error writing '%s': %s\n" % (hap_fn, inst))
            return 1
    return 0
