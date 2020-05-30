"""

preprocess.py

This directory contains functions and classes for preprocessing raw data into
matrices and vectors we can run our EM algorithm on. A big part of this is to
reduce the overall size of our input matrixes to make our algorithm run quickly
and without using too much memory.

Sam Vohr (svohr@soe.ucsc.edu)
Wed Apr  6 10:48:21 PDT 2016

"""

import sys
import math
import collections
import numpy

from mixemt import phylotree


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
    others are "unexpected". Finally, we will now make use of mutation counts
    available from Phylotree to adjust the probabilities based on how "mutable"
    a site appears to be.
    """
    def __init__(self, refseq, phylo, mut_wt=0.01, mut_max=0.5):
        """ Initialzes the HapVarBaseMatrix. """
        self.refseq = refseq
        self.phylo = phylo
        self.mut_wt = mut_wt
        self.mut_max = mut_max
        self.markers = dict()

        self.mut_prob = dict()
        for pos in self.phylo.variants:
            self.mut_prob[pos] = min(self.mut_max,
                                     self.mut_wt
                                     * sum(self.phylo.variants[pos].values()))

        self.add_hap_markers(phylo.hap_var)
        return

    def add_hap_markers(self, hap_var):
        """
        Populates the table with bases that are associated with haplogroups.
        """
        for hap in hap_var:
            self.markers[hap] = dict()
            for var in hap_var[hap]:
                pos = phylotree.pos_from_var(var)
                der = phylotree.der_allele(var)
                if der != self.refseq[pos]:
                    self.markers[hap][pos] = der
        return

    def _prob(self, hap_pos, pos, base):
        """
        Returns the "probability" of observing this base, at this position in
        this haplogroup, given the dictionary that describes the derived
        alleles of this haplogroup.
        """
        if pos in hap_pos:
            # Does this haplogroup carry a derived base?
            if hap_pos[pos] == base:
                # Is it the one we observed?
                return 1.0 - self.mut_prob[pos]
        else:
            # No derived base at this position, does our observed match ref.
            if self.refseq[pos] == base:
                return 1.0 - self.mut_prob[pos]
        return self.mut_prob[pos] / 3.0

    def prob_for_vars(self, hap, pos_obs):
        """
        Finds the probability of observing the read signature from the given
        haplotype by calculating the product of the probability of all
        hap, pos, base tuples it represents.
        """
        total = 0
        hap_pos = self.markers[hap]
        for pos, obs in pos_obs:
            total += math.log(self._prob(hap_pos, pos, obs))
        return total


def process_reads(alns, var_pos, min_mq, min_bq):
    """
    Read alignments from an iterator/iterable and makes observations for known
    variant positions. Each read is simplified into a signature of observed
    base per variant site.

    Args:
        alns: An iterator/iterable of pysam AlignedSegments
        var_pos: A sorted list of all known variant positions from Phylotree.
        min_mq: The minimum mapping quality score required for a read to be
                considered.
        min_bq: The minimum base quality score required for a base observation.

    Returns:
        read_obs: A dictionary that maps AlignedSegment query_names to
                  dictionaries mapping reference positions to base
                  observations.
    """
    read_obs = collections.defaultdict(dict)
    var_pos = set(var_pos)
    for aln in alns:
        if aln.mapping_quality >= min_mq:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                qpos = int(qpos)
                rpos = int(rpos)
                if (aln.query_qualities is None
                        or aln.query_qualities[qpos] >= min_bq):
                    if rpos in var_pos:
                        obs = aln.query_sequence[qpos].upper()
                        # Add this to the list, if this is a known var site.
                        if rpos in read_obs[aln.query_name]:
                            # check if this position has been observed before.
                            if read_obs[aln.query_name][rpos] != obs:
                                read_obs[aln.query_name][rpos] = "N"
                        else:
                            read_obs[aln.query_name][rpos] = obs
    # Finished, do one pass to remove Ns
    for aln_id in read_obs:
        read_obs[aln_id] = {pos:base for pos, base in
                            list(read_obs[aln_id].items()) if base != 'N'}
    return read_obs


def read_signature(obs_by_pos):
    """
    Returns a string made from the positions where bases were able to be
    observed from the read.
    """
    return ','.join(["%d:%s" % (pos, obs_by_pos[pos])
                     for pos in sorted(obs_by_pos)])


def pos_obs_from_sig(read_sig):
    """
    Returns a list of position, observation pairs described in the read
    signature string.
    """
    def pos_obs(var):
        """ Splits var into a int position and a string base. """
        pos, obs = var.split(':')
        return int(pos), obs
    return [pos_obs(var) for var in read_sig.split(',')]


def reduce_reads(read_obs):
    """
    Takes a dictionary of read IDs and base observations at variable positions
    and returns a dictionary observation signatures (a string describing the
    base observations made at variant sites) that maps to a list of read IDs
    that carry that signature.
    """
    read_sigs = collections.defaultdict(list)
    for read_id in read_obs:
        sig = read_signature(read_obs[read_id])
        read_sigs[sig].append(read_id)
    return read_sigs


def build_em_matrix(refseq, phylo, reads, haplogroups, args):
    """
    Returns the matrix that describes the probabiliy of each read
    originating in each haplotype.
    """
    hvb_mat = HapVarBaseMatrix(refseq, phylo)
    read_hap_mat = numpy.empty((len(reads), len(haplogroups)))

    if args.verbose:
        sys.stderr.write('Building EM input matrix...\n')

    for i in range(len(reads)):
        pos_obs = pos_obs_from_sig(reads[i])
        for j in range(len(haplogroups)):
            read_hap_mat[i, j] = hvb_mat.prob_for_vars(haplogroups[j], pos_obs)
        if args.verbose and (i + 1) % 500 == 0:
            sys.stderr.write('  processed %d fragments...\n' % (i + 1))

    if args.verbose:
        sys.stderr.write('Done.\n\n')

    return read_hap_mat


def build_em_input(bamfile, refseq, phylo, args):
    """
    Builds the matrix that describes the "probability" that a read originated
    from a specific haplogroup.
    """
    var_pos = phylo.get_variant_pos()

    input_alns = bamfile.fetch()
    read_obs = process_reads(input_alns, var_pos, args.min_mq, args.min_bq)
    read_sigs = reduce_reads(read_obs)

    if args.verbose:
        sys.stderr.write('Using %d aligned fragments (MQ>=%d) '
                         '(%d distinct sub-haplotypes)\n\n'
                         % (len(read_obs), args.min_mq, len(read_sigs)))

    # This is now the order we will be using for the matrix.
    haplogroups = sorted(phylo.hap_var)
    reads = sorted(read_sigs)
    weights = numpy.array([len(read_sigs[r]) for r in reads])

    em_matrix = build_em_matrix(refseq, phylo, reads, haplogroups, args)

    # make a list mapping matrix indexes to read IDs from bam.
    read_ids = [read_sigs[reads[i]] for i in range(len(reads))]

    return em_matrix, weights, haplogroups, read_ids


def reduce_em_matrix(em_mat, haplogroups, contrib_props):
    """
    Takes the matrix used by the EM algorithm, the column haplogroup labels,
    and the table of identified contributors and returns a new matrix made
    up of only the haplogroups that have passed the contributor filtering
    steps.

    Args:
        em_mat: a numpy matrix
        haplogroups: a list of strings for every column in em_mat
        contrib_props: a list of lists for each identified contributor
                       containing the hap#, haplogroup and initial proportion
                       estimate for each contributor.
    Returns:
        A new matrix made up of only the columns listed in contrib_props and
        a new list of labels for the simplified matrix.
    """
    haps_to_keep = {con[1] for con in contrib_props}
    indexes = [i for i in range(len(haplogroups))
               if haplogroups[i] in haps_to_keep]
    new_haps = [haplogroups[i] for i in indexes]
    return em_mat[:, indexes], new_haps


def main():
    """ Main function for simple testing. """
#   args = {'min_mq':30, 'min_bq':30, 'verbose':True}
#   if len(sys.argv) > 2:
#       phy_fn = sys.argv[1]
#       bam_fn = sys.argv[2]
#       ref_in = pysam.FastaFile('../ref/RSRS.mtDNA.fa')
#       refseq = ref_in.fetch(ref_in.references[0])
#       with open(phy_fn, 'r') as phy_in:
#           var_pos, hap_var = phylotree.read_phylotree(phy_in,
#                                                       False, False, False)
#       with pysam.AlignmentFile(bam_fn, 'rb') as samfile:
#           build_em_input(samfile, refseq, var_pos, hap_var, args)
#   else:
        # Entirely fake data.
#       ref = "GAAAAAAAA"
#       var_pos = range(1, 9)
#       hap_var = dict({'A':['A1T', 'A3T'],
#                       'B':['A2T', 'A4T', 'A5T', 'A7T'],
#                       'C':['A2T', 'A5T'],
#                       'D':['A2T', 'A4T', 'A6T', 'A8T'],
#                       'E':['A2T', 'A3T', 'A4T', 'A6T']})
#       reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
#                     "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
#                     "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
#       haps = list('ABCDE')
#       print build_em_matrix(ref, hap_var, reads, haps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
