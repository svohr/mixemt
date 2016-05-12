"""
observe.py

This module contains functions for building a table of counts of all bases
observed at all reference positions. This table will be used to associate novel
variants observed in the sample with the haplogroups and variants known from
the tree. We can also use this table to confirm that variants from the tree
are observed within this sample.

Sam Vohr (svohr@soe.ucsc.edu)
Wed May 11 10:47:05 PDT 2016

"""

import collections

def build_obs_table(samfile, min_mq=30, min_bq=30):
    """
    Returns a table of counts for every base observed at every reference
    position in the fragments contained in the provided samfile.
    """
    obs_tab = collections.defaultdict(collections.Counter)

    for aln in samfile.fetch():
        if aln.mapping_quality >= min_mq:
            for qpos, rpos in aln.get_aligned_pairs(matches_only=True):
                qpos = int(qpos)
                rpos = int(rpos)
                if (aln.query_qualities is None or 
                    aln.query_qualities[qpos] >= min_bq):
                    obs = aln.query_sequence[qpos].upper()
                    obs_tab[rpos][obs] += 1
    return obs_tab

