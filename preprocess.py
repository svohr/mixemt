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
import numpy
import pysam
import collections

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
        self.markers = dict()

        if hap_var is not None:
            self.add_hap_markers(hap_var)
        return

    def add_hap_markers(self, hap_var):
        """
        Populates the table with bases that are associated with haplogroups.
        """
        for hap in hap_var:
            for var in hap_var[hap]:
                pos = phylotree.pos_from_var(var)
                der = phylotree.der_allele(var)
                if der != self.refseq[pos]:
                    self.markers[(hap, pos)] = der
        return

    def prob(self, hap, pos, base):
        """
        Returns the "probability" of observing this base, at this position
        in this haplogroup.
        """
        if (hap, pos) in self.markers:
            # Does this haplogroup carry a derived base?
            if self.markers[hap, pos] == base:
                # Is it the one we observed?
                return self.exp
        else:
            # No derived base at this position, does our observed match ref.
            if self.refseq[pos] == base:
                return self.exp
        return self.unexp

    def prob_for_vars(self, hap, pos_obs):
        """
        Finds the probability of observing the read signature from the given
        haplotype by calculating the product of the probability of all 
        hap, pos, base tuples it represents.
        """
        total = 1
        for pos, obs in pos_obs: 
            total *= self.prob(hap, pos, obs)
        return total


def build_hap_var_base_mat(refseq, hav_var, exp=0.991, unexp=0.003):
    """
    Build a big matrix of haplogroups, variant sites and bases. 
    """
    return


def process_reads(refseq, samfile, var_pos, min_mq, min_bq):
    """
    Reads in a set of reads from a SAM/BAM file and makes observations for
    known variant positions. Each read is simplified into a signature of
    observed base per variant site.
    """
    read_obs = collections.defaultdict(dict)
    for aln in samfile:
        if aln.mapping_quality >= min_mq:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                qpos = int(qpos)
                rpos = int(rpos)
                if qpos in var_pos:
                    if (aln.query_qualities is None or 
                        aln.query_qualities[qpos] >= min_bq):
                        # Add this to the list
                        obs = aln.query_sequence[qpos].upper()
                        if rpos in read_obs[aln.query_name]:
                            # check if this position has been observed before.
                            if read_obs[aln.query_name][rpos] != obs:
                                read_obs[aln.query_name][rpos] = "N"
                        else:
                            read_obs[aln.query_name][rpos] = obs
    print len(read_obs)
    # Finished, do one pass to remove Ns
    for aln_id in read_obs:
        read_obs[aln_id] = {pos:base for pos, base in 
                            read_obs[aln_id].items() if base != 'N'}
    print len(read_obs)
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


def build_em_matrix(refseq, hap_tab, reads, haplogroups):
    """ 
    Returns the matrix that describes the probabiliy of each read 
    originating in each haplotype. 
    """
    hvb_mat = HapVarBaseMatrix(refseq, hap_tab)
    read_hap_mat = numpy.empty((len(reads), len(haplogroups)))

    for i in xrange(len(reads)):
        pos_obs = pos_obs_from_sig(reads[i]) 
        for j in xrange(len(haplogroups)):
            read_hap_mat[i, j] = hvb_mat.prob_for_vars(haplogroups[j], pos_obs) 
        if i % 1000 == 0:
            print i, read_hap_mat[i, 0]
    return read_hap_mat


def build_em_input(samfile, refseq, var_pos, hap_tab):
    """
    Builds the matrix that describes the "probability" that a read originated
    from a specific haplogroup.
    """
    read_obs = process_reads(refseq, samfile, var_pos, 30, 30)
    read_sigs = reduce_reads(read_obs)

    # This is now the order we will be using for the matrix.
    haplogroups = sorted(hap_tab)
    reads = sorted(read_sigs)
    weights = numpy.array([len(read_sigs[r]) for r in reads])
    
    em_matrix = build_em_matrix(refseq, hap_tab, reads, haplogroups)
    return em_matrix, weights, haplogroups, reads, read_sigs


def main():
    """ Main function for simple testing. """
    if len(sys.argv) > 2:
        phy_fn = sys.argv[1]
        bam_fn = sys.argv[2]
        ref_in = pysam.FastaFile('../ref/RSRS.mtDNA.fa')
        refseq = ref_in.fetch(ref_in.references[0])
        with open(phy_fn, 'r') as phy_in:
            var_pos, hap_var = phylotree.read_phylotree(phy_in,
                                                        False, False, False)
        with pysam.AlignmentFile(bam_fn, 'rb') as samfile:
            build_em_input(samfile, refseq, var_pos, hap_var)
    else:
        # Entirely fake data.
        ref = "GAAAAAAAA"
        var_pos = range(1, 9)
        hap_var = dict({'A':['A1T', 'A3T'],
                        'B':['A2T', 'A4T', 'A5T', 'A7T'],
                        'C':['A2T', 'A5T'],
                        'D':['A2T', 'A4T', 'A6T', 'A8T'],
                        'E':['A2T', 'A3T', 'A4T', 'A6T']})
        reads = list(["1:A,2:T,3:A", "2:T,3:A", "3:A,4:T,5:T", "5:T,6:A",
                      "6:A,7:T", "6:A,7:T,8:A", "7:T,8:A", "4:T,5:T",
                      "1:A,2:T,3:T,4:T", "5:A,6:T,7:A,8:A"])
        haps = list('ABCDE')
        print build_em_matrix(ref, hap_var, reads, haps)
    return 0


if __name__ == "__main__":
    sys.exit(main())
