"""

This file contains functions for reading in phylotree variants and haplogroups
and building a matrix to represent these markers.

Mon Apr  4 12:00:54 PDT 2016

"""

import sys
import numpy
import bisect


def pos_from_var(var):
    """ Returns the position of the variant """
    if var.startswith('('):
        var = var[1:-1]
    var = var.rstrip('!')
    return int(var[1:-1])


def is_snp(var):
    """ 
    Returns true if var is a SNP or False if it is an indel
    """
    if '.' in var or 'd' in var:
        return False # Variant is an indel
    return True 


def rm_snp_annot(var):
    """ 
    Returns the SNP variant string, nicely formatted with annotation stripped.
    """
    if var.startsiwth('('):
        var = var[1:-1]
    var.rstrip('!')
    return var.upper()


def read_phy_line(line):
    """
    Reads a single comma-separated line from phylotree and returns
    the indentation level, the haplogroup id, and variants.
    """
    items = line.rstrip().split(',')
    level = 0
    while items[level] == '':
        level += 1
    hap_id = items[level]
    variants = items[level + 1].split()
    return level, hap_id, variants


def summarize_vars(var_stack):
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
            if var not in summed_vars:
                summed_vars[pos] = var
    return [summed_vars[pos] for pos in sorted(summed_vars)]


def read_phylotree(phy_in, leaves=False):
    """
    Reads input from phylotree table that has been converted into 
    comma-separated values and produces: 
    1) a list of SNP variant positions
    2) a table of haplogroups with associated SNP variants.
    """
    var_pos = set()
    hap_tab = dict()
    cur = 0
    var_stack = list()
    for line in phy_in:
        level, hap_id, raw_var = read_phy_line(line)
        variants = [var for var in raw_var if is_snp(var)]
        if (not leaves or level <= cur) and cur > 0:
            # Store previous entry checking if it was a leaf.
            hap_tab[var_stack[-1][0]] = summarize_vars(var_stack)
        while cur > 0 and cur >= level:
            var_stack.pop()
            cur -= 1
        var_stack.append((hap_id, variants))
        cur += 1
    return hap_tab


def main():
    """ Simple test of phylotree functions. """
    if len(sys.argv) > 1:
        phy_fn = sys.argv[1]
        with open(phy_fn, 'r') as phy_in:
            hap_var = read_phylotree(phy_in)
            for hap in hap_var:
                print hap, ','.join(hap_var[hap])
    return 0


if __name__ == "__main__":
    sys.exit(main())
