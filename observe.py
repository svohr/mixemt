"""

This file contains a class for building a table of observations for each
position of a reference genome. This is meant to be a little more robust than
passing around a dictionary of Counters. We will keep track of strand
information as well as gaps.

Sam Vohr (svohr@soe.ucsc.edu)
Mon Jul  4 11:06:25 PDT 2016

"""

import collections


class ObservedBases(object):
    """
    Class for keeping track of observations at reference positions. To do this,
    we use a default dictionary of counters to keep track of reference
    position obervations as single-character strings. The values that these
    may be are 'A', 'C', 'G', 'T' for forward strand observations,
    'a', 'c', 'g', 't' for reverse strand observations, 'N' for unknown base
    and '-' for a gap.

    Atributes:
        obs_tab: A default dict of Counters that contain the observed bases
                 for each observed reference position.
        min_map_qual: Minimum mapping quality to use for reads.
        min_base_qual: Minimum base quality required for observation
    """
    def __init__(self, alns=None, mapq=30, baseq=30):
        """
        Initial ObservedBases as empty, or from a set of alignments
        Args:
            alns: An iterable of pysam AlignedSegments.
        """
        self.obs_tab = collections.defaultdict(collections.Counter)
        self.min_map_qual = mapq
        self.min_base_qual = baseq
        if alns:
            self.update(alns)

    def update(self, alns):
        """
        Update the counts based on these alignments.

        Args:
            alns: An iterable of pysam AlignedSegments.
        Returns:
            nothing
        """
        for aln in alns:
            self._observe_aln(aln)
        return

    def _observe_aln(self, aln):
        """
        Updates obs_tab with the observations from this read.

        Args:
            aln: a pysam AlignedSegment object
        Returns:
            nothing
        """
        if aln.mapping_quality >= self.min_map_qual:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=False):
                if rpos is None:
                    # no matched reference position, skip
                    continue
                rpos = int(rpos)
                obs = ' '
                if qpos is not None:
                    if (aln.query_qualities is None
                      or aln.query_qualities[qpos] >= self.min_base_qual):
                        qpos = int(qpos)
                        obs = aln.query_sequence[qpos].upper()
                    else:
                        obs = 'N'
                    if aln.is_reverse:
                        obs = obs.lower()
                else:
                    obs = '-'
                    if aln.is_reverse:
                        obs = '+'
                self.obs_tab[rpos][obs] += 1
        return

    def obs_at(self, pos, base=None, stranded=False):
        """
        Returns the total number of times the base (on both strands) was
        observed at this position in the reference.

        Args:
            pos: the reference position we are interested in
            base: the DNA base we are asking about, if not set return the
                  counter for the position.
            stranded: A flag for whether the total base count should be
                      returned or a tuple with (for, rev) counts
        Returns:
            A counter if base is None
            An int count if stranded is False
            A tuple of 2 ints if stranded is True
        Raises:
            ValueError if base is not one of 'ACGTN-acgtn+'
        """
        if base is None:
            return self.obs_tab[pos]
        if base.upper() in 'ACGTN':
            for_count = self.obs_tab[pos][base.upper()]
            rev_count = self.obs_tab[pos][base.lower()]
        elif base in '-+':
            for_count = self.obs_tab[pos]['-']
            rev_count = self.obs_tab[pos]['+']
        else:
            raise ValueError("Bad base: %s" % (base))
        if stranded:
            return (for_count, rev_count)
        return for_count + rev_count
