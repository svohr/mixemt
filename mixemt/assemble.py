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

import sys
import operator
import collections
import itertools
import numpy
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mixemt import phylotree
from mixemt import observe


def get_contributors(phylo, obs_tab, haplogroups, wts, em_results, args):
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

    Args:
        phylo: The Phylotree object used for this analysis
        obs_tab: Table of base observations for positions in the reference.
        haplogroups: List of haplogroups IDs in the order they appear as
                     columns in the em_results matrix
        wts: Vector of interger counts of the number of fragments where each
             sub-haplotype in the matrix
        em_results tuple: results from EM algorithm; includes...
            props: Array of stimatated proportions of haplogroups
            read_hap_mat: Conditional read/haplogroup probability matrix
        args: arguments namespace from argparse.

    Returns:
        A list of lists, one for each contributor, containing:
          - A friendly hap# ID assigned based on the proportion estimate.
          - The haplogroup ID for the contributor.
          - The proprotion estimate from EM.

    Raises:
        ValueError if user-specified haplogroup ID is not found.
    """
    props, read_hap_mat = em_results

    contributors = list()
    if not args.contributors:
        contributors = _find_contribs_from_reads(read_hap_mat, wts, args)
    else:
        # override contributor identification steps.
        con_haps = args.contributors.split(',')
        for con in con_haps:
            try:
                index = haplogroups.index(con)
                contributors.append(index)
            except ValueError:
                raise ValueError("Unknown haplogroup '%s'" % (con))

    # Extract contributor proportions
    contrib_prop = [[haplogroups[con], props[con]] for con in contributors]
    contrib_prop.sort(key=operator.itemgetter(1), reverse=True)

    if args.var_check and not args.contributors:
        # Remove haplogroups with minimal variant support.
        contrib_prop = _check_contrib_phy_vars(phylo, obs_tab,
                                               contrib_prop, args)

    # Add a friendly name, not associated with a haplogroup
    name_fmt = "hap%%0%dd" % (len(str(len(contrib_prop) + 1)))
    hap_num = 1
    for con in contrib_prop:
        con.insert(0, name_fmt % (hap_num))
        hap_num += 1

    return contrib_prop


def _find_contribs_from_reads(read_hap_mat, wts, args):
    """
    Returns a list of indexes to haplogroups that pass the minimum required
    number of reads filters

    Args:
        read_hap_mat: The final conditional probability matrix from EM
        wts: The vector of read counts for each distinct subhaplotype
        args: arguments namespace from argparse.
    Returns:
        A list of column (haplogroup) indexes in read_hap_mat
    """
    best_haps = numpy.argmax(read_hap_mat, 1)
    vote_count = collections.defaultdict(int)

    for hap, read_count in itertools.izip(best_haps, wts):
        vote_count[hap] += read_count

    contributors = [con for con in vote_count
                    if vote_count[con] >= args.min_reads]
    return contributors


def _check_contrib_phy_vars(phylo, obs_tab, contrib_prop, args):
    """
    Checks if each candidate contributor from contribs passes our variant base
    check. The strategy for this is to start with the highest estimated
    contributors and an empty list of variant positions. For each contributor,
    we identify the variant bases that are unique from the previous candidates.
    We check the observation table to verify that those bases are observed in
    the sample.

    Args:
        phylo: the phylotree object that holds the variant information for
               these haplogroups
        obs_tab: Table of base observations for positions in the reference.
        contrib_prop: A list of lists, one for each contributor, containing:
          - The haplogroup ID for the contributor.
          - The proprotion estimate from EM
        args: argparse Namespace with user specified values for:
            min_var_reads: The minimum number of observations required to call
                           a base as present in the mixture sample (int)
            var_fraction: The minimum fraction of defining variants required
                          to be observed to call a haplogroup a contributor
                          (float)
            var_count: Call a haplogroup a contributor if there the number of
                       observed variants is equal or greater than min_count
                       (int)
    Returns:
        contrib_prop, with haplogroups that do not pass filters removed.
    """
    used_vars = set()
    ignore_haps = set()

    if args.verbose:
        sys.stderr.write("Checking diagnostic variants:\n")

    for hap, _ in contrib_prop:
        # get variant for this haplogroup
        uniq_vars = set([(phylotree.pos_from_var(var),
                          phylotree.der_allele(var))
                         for var in phylo.hap_var[hap]])
        uniq_vars -= used_vars
        if args.verbose:
            sys.stderr.write("%s (%d unique variants)\n"
                             % (hap, len(uniq_vars)))

        found_vars = set()
        for pos, der in sorted(uniq_vars):
            if args.verbose:
                var = "%d%s" % (pos + 1, der)
                sys.stderr.write("  %s: %d/%d\n" % (var.rjust(6),
                                                    obs_tab.obs_at(pos, der),
                                                    obs_tab.total_obs(pos)))
            if obs_tab.obs_at(pos, der) >= args.min_var_reads:
                found_vars.add((pos, der))
        if ((len(uniq_vars) == 0)
                or (args.var_count is not None
                    and len(found_vars) >= args.var_count)
                or (float(len(found_vars)) / len(uniq_vars)
                    >= args.var_fraction)):
            if args.verbose:
                sys.stderr.write("Keeping '%s': "
                                 "%d/%d unique variant bases observed at "
                                 "least %d times.\n"
                                 % (hap, len(found_vars), len(uniq_vars),
                                    args.min_var_reads))
            # Looks good, these variants can't be used again.
            used_vars.update(found_vars)
            # Also add the ancestral bases for this haplogroup so we do not
            # mistake backmutations in another haplogroup as a novel allele.
            used_vars.update(phylo.get_ancestral(hap))
        else:
            if args.verbose:
                sys.stderr.write("Ignoring '%s': "
                                 "only %d/%d unique variant bases observed.\n"
                                 % (hap, len(found_vars), len(uniq_vars)))
            ignore_haps.add(hap)

    pass_contribs = [con for con in contrib_prop if con[0] not in ignore_haps]

    return pass_contribs


def update_contribs(contribs, em_results, haps):
    """
    Takes an existing contributor table and appends the refined proportion
    estimates to the end of each entry.

    Args:
        contribs: The table of contributors as returned by get_contributors()
        em_results: A tuple containing the final mixture proportion vector and
                    the read/haplogroup conditional probabilities from the
                    second round of EM.
        haps: A list haplogroup IDs/label for each column in the matrix.
    Returns:
        The update contributor table
    """
    props, _ = em_results
    props_by_hap = {haps[i]:props[i] for i in xrange(len(haps))}
    for con in contribs:
        haplogroup = con[1]
        con[2] = props_by_hap[haplogroup]
    return contribs


def assign_reads(bamfile, contribs, em_results, haps, reads, args):
    """
    Assigns reads from the BAM file to contributors using the results from EM.

    Args:
        bamfile: The original bamfile containing the reads to be assigned to
                 contributors.
        contribs: The table of contributors as returned by get_contributors()
        em_results: A tuple containing the final mixture proportion vector and
                    the read/haplogroup conditional probabilities from EM.
        haps: A list haplogroup IDs/label for each column in the matrix.
        reads: A list of lists where each row index in the matrix is associated
               with the read ID(s) it represents.
        args: From argparse.

    Returns: A table that maps contributor names to a list of AlignedSegments
             reprenting the reads that have been assigned to that contributor.
    """
    contrib_reads = collections.defaultdict(list)
    read_to_con = dict()
    con_read_indexes = assign_read_indexes(contribs, em_results,
                                           haps, reads, args.min_fold)
    for con in con_read_indexes:
        con_read_ids = get_contrib_read_ids(con_read_indexes[con], reads)
        for read_id in con_read_ids:
            read_to_con[read_id] = contrib_reads[con]

    for aln in bamfile.fetch():
        if aln.query_name in read_to_con:
            read_to_con[aln.query_name].append(aln)

    return contrib_reads


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


def assign_read_indexes(contribs, em_results, haps, reads, min_fold):
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

    Args:
        contribs: The table of contributors as returned by get_contributors()
        em_results: A tuple containing the final mixture proportion vector and
                    the read/haplogroup conditional probabilities from EM.
        haps: A list haplogroup IDs/label for each column in the matrix.
        reads: The sub-haplotype identifier (read signature) for each row
               in the matrix.
        min_fold: The minimum odds ratio between top two contributors to assign
                  a read to a contributor.
    Returns: A dictionary mapping contributor names to a list of indexes to
             entries in 'reads'
    """
    props, read_hap_mat = em_results
    contrib_reads = collections.defaultdict(set)

    # move mixture proportions and odds ratio cutoff into log, just once.
    min_fold = numpy.log(min_fold)
    log_props = numpy.log(props)

    index_to_hap = {haps.index(group):hap_n for hap_n, group, _ in contribs}
    con_indexes = set(index_to_hap.keys())
    for read_i in xrange(len(reads)):
        if len(contribs) > 1:
            read_probs = read_hap_mat[read_i, ] - log_props
            best_hap, next_hap = _find_best_n_for_read(read_probs,
                                                       con_indexes,
                                                       top_n=2)[0:2]

            if read_probs[best_hap] - read_probs[next_hap] >= min_fold:
                contrib_reads[index_to_hap[best_hap]].add(read_i)
            else:
                contrib_reads['unassigned'].add(read_i)
        else:
            # Only 1 identified contributor, use the first field of the first
            # entry in the contributors table.
            contrib_reads[contribs[0][0]].add(read_i)
    return contrib_reads


def get_contrib_read_ids(indexes, reads):
    """
    Takes a set of indexes from assign_reads and the list of read signatures
    plus the dictionary mapping signatures to aligned read IDs and returns
    the set of corresponding aligned read IDs (BAM query IDs).
    """
    hap_read_ids = set()
    for read_idx in indexes:
        for read_id in reads[read_idx]:
            hap_read_ids.add(read_id)
    return hap_read_ids


def write_haplotypes(bamfile, contrib_reads, args):
    """
    Writes a new BAM file based on the original 'bamfile' for each contributor
    described in the table 'contrib_reads' that maps a haplotype id to a list
    of pysam AlignedSegment objects.

    Args:
        bamfile: a pysam AlignmentFile for the original BAM file.

        contrib_reads: a dictionary mapping our haplotype label to a list of
            pysam AlignedSegments for the reads that have been assigned to this
            contributor.

        args: arguments namespace from argparse. Used to get filename prefix
            for new BAM files and whether verbose mode is enabled.

    Returns:
        0 if all files written successfully, 1 otherwise.
    """
    # Set up for opening new bam files.
    ext = 'bam'
    mode = 'wb'

    if args.verbose:
        sys.stderr.write('\nWriting haplotype alignment files...\n')

    for contrib in contrib_reads:
        if not contrib_reads[contrib]:
            # No reads assigned to this contributor
            continue
        hap_fn = "%s.%s.%s" % (args.out_prefix, contrib, ext)
        try:
            hap_bamfile = pysam.AlignmentFile(hap_fn, mode, template=bamfile)
            written = 0
            for aln in contrib_reads[contrib]:
                hap_bamfile.write(aln)
                written += 1
            if args.verbose:
                sys.stderr.write('  Wrote %d aligned segments to %s\n'
                                 % (written, hap_fn))
            hap_bamfile.close()
        except (ValueError, IOError) as inst:
            sys.stderr.write("Error writing '%s': %s\n" % (hap_fn, inst))
            return 1
    return 0


def write_consensus_seqs(refseq, contrib_props, contrib_reads, args):
    """
    Generates consensus sequences for each contributor from the assigned reads
    for output in FASTA format and writes them out.

    Args:
        refseq: The reference sequence to which the fragments were aligned.
        contrib_props: A list of lists containing for each contributor
                       - contributor ID (hap#)
                       - haplogroup
                       - proportion in mixture (not used).
        contrib_reads: A table mapping hap# IDs to lists of pysam
                       AlignedSegments + an entry of unassigned.
        args: The argument values from mixemt's argparse results.
    Returns:
        nothing
    """
    with open("%s.fa" % (args.cons_prefix), 'w') as fa_out:
        seqs_to_write = list()
        for con, hap, _ in contrib_props:
            seq = call_consensus(refseq, contrib_reads[con],
                                 1, args, strict=False)
            rec = SeqRecord(Seq(seq), id=con, description=hap)
            seqs_to_write.append(rec)
        if 'unassigned' in contrib_reads:
            seq = call_consensus(refseq, contrib_reads['unassigned'],
                                 1, args, strict=False)
            rec = SeqRecord(Seq(seq), id='unassigned', description='')
            seqs_to_write.append(rec)
        SeqIO.write(seqs_to_write, fa_out, 'fasta')
    return


def call_consensus(refseq, alns, min_cov, args, strict=True):
    """
    Generates a consensus sequence based on the list of AlignedSegments.

    Args:
        refseq: The reference sequence to which the fragments were aligned.
        alns: A list of pysam AlignedSegments
        min_cov: minimum coverage required to call a base.
        args: The argument values from mixemt's argparse results.
        strict: Whether to call a strict consensus or a majority base
    Returns:
        A string representing the consensus of the alignments in alns
    """
    def consensus_base(base_counts):
        """
        Given a Counter for a reference position, return the base that
        represents the consensus:
        N if coverage requirement is not met or if bases disagree
        or the observed base if all observations agree.
        """
        base_counts['N'] = 0 # first, ignore an missing observations
        total_obs = sum(base_counts.values())
        if total_obs < min_cov:
            return 'N'
        base, count = base_counts.most_common(1)[0]
        if strict and count != total_obs:
            return 'N'
        else:
            return base
    if not alns:
        # Sometimes alns can be empty.
        return ""
    obs_tab = observe.ObservedBases(alns, args.min_mq, args.min_bq)
    cons_bases = [consensus_base(obs_tab.obs_at(pos))
                  for pos in xrange(len(refseq))]
    return str(''.join(cons_bases))


def find_new_variants(refseq, contrib_reads, args):
    """
    Produces a dictionary that maps reference positions and base tuples to the
    contributing haplotype with which it is associated.

    Args:
        refseq: The reference sequence to which the fragments were aligned.
        contrib_reads: dictionary mapping hapN ids to lists of pysam
                       AlignedSegments
        args: The argument values from mixemt's argparse results.
    Returns:
        A dictionary mapping (ref pos., base) to a contributor hapN id.
    """
    new_vars = dict()
    contrib_cons = {con:call_consensus(refseq, contrib_reads[con],
                                       args.cons_cov, args, strict=True)
                    for con in contrib_reads if con != 'unassigned'}
    if not contrib_cons:
        return {} # no contributors, no new variants.

    min_cons_len = min([len(contrib_cons[cons]) for cons in contrib_cons])
    for pos in xrange(min_cons_len):
        if any([cons[pos] in 'N-' for cons in contrib_cons.values()]):
            # skip Ns or reference position skips
            continue
        # base observed in all consensus sequences
        for hap in contrib_cons:
            base = contrib_cons[hap][pos]
            base_is_uniq = all([contrib_cons[con][pos] != base
                                for con in contrib_cons if con != hap])
            if base_is_uniq:
                new_vars[pos, base] = hap
    return new_vars


def assign_reads_from_new_vars(contrib_reads, new_variants, args):
    """
    Assigns reads from the 'unassigned' list in contrib reads to contributors
    using a dictionary mapping new variant positions and base observations to
    contributors.

    Args:
        contrib_reads: a dictionary mapping hapN ids to lists of pysam
                       AlignedSegments, including an entry for "unassigned"
        new_variants: a dictionary mapping (ref. position, base) tuples to
                       hapN IDs.
        args: The argument values from mixemt's argparse results.
    Returns:
        Updated version of contrib_reads.
    """
    temp_assigns = collections.defaultdict(set)

    for aln in contrib_reads['unassigned']:
        if aln.mapping_quality >= args.min_mq:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                qpos = int(qpos)
                rpos = int(rpos)
                if (aln.query_qualities is None or
                        aln.query_qualities[qpos] >= args.min_bq):
                    base = aln.query_sequence[qpos].upper()
                    if (rpos, base) in new_variants:
                        contrib = new_variants[(rpos, base)]
                        temp_assigns[aln.query_name].add(contrib)

    # Find reads that have been assigned to exactly 1 contributor
    uniq_assigns = {qname:list(contrib)[0]
                    for qname, contrib in temp_assigns.iteritems()
                    if len(contrib) == 1}

    # Move reads that can be assigned to new contributor assignment
    for aln in contrib_reads['unassigned']:
        if aln.query_name in uniq_assigns:
            contrib_reads[uniq_assigns[aln.query_name]].append(aln)

    # Remove the reads that have been moved from "unassigned"
    contrib_reads['unassigned'] = [aln for aln in contrib_reads['unassigned']
                                   if aln.query_name not in uniq_assigns]
    return contrib_reads


def extend_assemblies(refseq, contrib_reads, args):
    """
    Iteratively tries to extend current assemblies by discovering novel
    variants from haplotype assembly consensus sequences, assigning reads to
    contributors based on new variants and repeating until convergence.

    Args:
        contrib_reads: A table mapping hap# IDs to lists of pysam
                       AlignedSegments + an entry of unassigned.
        args: The argument values from mixemt's argparse results.
    Returns:
        Updated version of contrib_reads.
    """
    last_unassigned = None
    unassigned = len(contrib_reads['unassigned'])
    run = 1

    if args.verbose:
        sys.stderr.write('\nAssembly extension step...\n')

    while last_unassigned != unassigned:
        new_variants = find_new_variants(refseq, contrib_reads, args)
        contrib_reads = assign_reads_from_new_vars(contrib_reads,
                                                   new_variants, args)
        last_unassigned = unassigned
        unassigned = len(contrib_reads['unassigned'])
        if args.verbose:
            sys.stderr.write(
                "  %d: %d/%d reads assigned using %d variants\n"
                % (run, last_unassigned - unassigned, last_unassigned,
                   len(new_variants)))
        run += 1

    if args.verbose:
        sys.stderr.write('\n')

    return contrib_reads
