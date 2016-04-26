"""

This file contains functions for reading in phylotree variants and haplogroups
and building a matrix to represent these markers.

Mon Apr  4 12:00:54 PDT 2016

"""

import sys
import collections


class Phylotree:
    """
    Class for representing haplogroups and the variants that define them in
    an explicit tree data structure. The general workflow of this is to load
    in the tree from a file based on the Phylotree mtDNA tree (phylotree.org),
    perform any filtering steps and produce a table of haplogroup IDs and
    the variants that define them.
    """
    class PhyloNode:
        """
        Class that represents a single node in the phylotree.
        """
        def __init__(self, hap='', parent=None):
            """ Initialize new PhyloNode """
            self.hap_id = hap
            self.parent = parent
            self.children = list()
            self.variants = list()

    def __init__(self):
        """ Initialize a blank Phylotree before reading from a file. """
        self.root = None
        self.var_pos = dict()


def pos_from_var(var):
    """ Returns the position of the variant """
    if var.startswith('('):
        var = var[1:-1]
    var = var.rstrip('!')
    return int(var[1:-1]) - 1 # 1-based to 0-based


def der_allele(var):
    """ Returns the derived allele of this variant. """
    var = var.rstrip(')!')
    return var[-1].upper()


def anc_allele(var):
    """ Returns the derived allele of this variant. """
    var = var.lstrip('(')
    return var[0].upper()


def is_snp(var):
    """ Returns true if var is a SNP or False if it is an indel """
    if '.' in var or 'd' in var:
        return False # Variant is an indel
    return True 


def is_unstable(var):
    """ Returns true if var is annotated as unstable. """
    if var[0] == '(':
        return True
    return False


def is_backmutation(var):
    """ Returns true if var is annotated as a backmutation. """
    if '!' in var:
        return True
    return False


def rm_snp_annot(var):
    """ 
    Returns the SNP variant string, nicely formatted with annotation stripped.
    """
    if var.startsiwth('('):
        var = var[1:-1]
    var.rstrip('!')
    return var.upper()


def _read_phy_line(line):
    """
    Reads a single comma-separated line from phylotree and returns the
    indentation level, the haplogroup id, and variants. Importantly, some nodes
    do not have IDs, so we need to check that we do not accidentally grab the
    variant list instead of the blank id.
    """
    items = line.rstrip().split(',')
    level = 0
    while items[level] == '':
        level += 1
    hap_id = items[level]
    while ' ' in hap_id:
        level -= 1
        hap_id = items[level]
    variants = items[level + 1].split()
    return level, hap_id, variants


def _summarize_vars(var_stack):
    """
    Takes a list of variant lists and produces a flat list of variants that
    are associated with the lineage leading to and including this haplogroup.
    Importantly, mutations that occur to the same position are masked by the
    most recent mutation (i.e. C152T will not be included if T152C occurs 
    farther down in the tree).
    """
    summed_vars = dict()
    for hap_id, variants in reversed(var_stack):
        # Go through the stack backwards and only add a variant if we have
        # not seen a variant at the same position.
        for var in variants:
            pos = pos_from_var(var) 
            if pos not in summed_vars:
                summed_vars[pos] = var
    return [summed_vars[pos] for pos in sorted(summed_vars)]


def _anon_hap_name(parent):
    """
    Returns an id for a haplogroup without a name
    """
    _anon_hap_name.counter[parent] += 1
    return "%s[%d]" % (parent, _anon_hap_name.counter[parent])
_anon_hap_name.counter = collections.defaultdict(int)


def _read_phylotree_csv(phy_in, leaves_only=False):
    """
    Reads input from phylotree table that has been converted into
    comma-separated values and produces a table of haplogroups with associated
    SNP variants. The table is constructed so that variants defining a parent
    haplogroup are also associated with all descendent haplogroups.
    """
    hap_tab = dict()
    cur = -1
    var_stack = list()
    for line in phy_in:
        level, hap_id, raw_var = _read_phy_line(line)
        variants = [var for var in raw_var if is_snp(var)]
        if (not leaves_only or level <= cur) and cur >= 0:
            # Store previous entry checking if it was a leaf.
            hap_tab[var_stack[-1][0]] = _summarize_vars(var_stack)
        while cur >= 0 and cur >= level:
            var_stack.pop()
            cur -= 1
        if hap_id == '':
            hap_id = _anon_hap_name(var_stack[-1][0])
        var_stack.append((hap_id, variants))
        cur += 1
    else:
        if (not leaves_only or level <= cur) and cur > 0:
            # Store previous entry checking if it was a leaf.
            hap_tab[var_stack[-1][0]] = _summarize_vars(var_stack)
    return hap_tab


def _flatten_var_pos(hap_tab, rm_unstable=False, rm_backmut=False):
    """ 
    Takes a dictionary of haplogroups and variants and produces a list of
    variant positions. Removes variants and positions that do not pass filters.
    """
    blacklist = set()
    var_pos = set()
    for hap in hap_tab:
        for var in hap_tab[hap]:
            pos = pos_from_var(var)
            if rm_unstable and is_unstable(var):
                blacklist.add(pos)
            elif rm_backmut and is_backmutation(var):
                blacklist.add(pos)
            else:
                var_pos.add(pos)
    if rm_unstable or rm_backmut:
        var_pos -= blacklist
        for hap in hap_tab:
            hap_tab[hap] = [var for var in hap_tab[hap] 
                            if pos_from_var(var) not in blacklist] 
    return list(sorted(var_pos))             


def read_phylotree(phy_in, 
                   leaves_only=False, rm_unstable=False, rm_backmut=False):
    """ 
    Reads input from phylotree table that has been converted into
    comma-separated values and produces a list of variant sites and table of
    haplogroups with associated SNP variants. Optionally removes sites that
    are annotated as unstable or contain backmutations.
    """
    hap_tab = _read_phylotree_csv(phy_in, leaves_only=leaves_only)
    var_pos = _flatten_var_pos(hap_tab, rm_unstable, rm_backmut)
    return var_pos, hap_tab


def main():
    """ Simple test of phylotree functions. """
    if len(sys.argv) > 1:
        phy_fn = sys.argv[1]
        with open(phy_fn, 'r') as phy_in:
            var_pos, hap_var = read_phylotree(phy_in, False, False, False)
        for hap in hap_var:
            print hap, ','.join(hap_var[hap])
            #print len(var_pos)
    return 0


if __name__ == "__main__":
    sys.exit(main())
