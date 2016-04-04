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


def pass_snp_filter(var, no_unstable=True, no_backmut=False):
    """ 
    Returns true if var is a SNP or False if it is an indel or carries
    annotation that we are filtering against (unstable or backmutation).
    """
    if '.' in var or 'd' in var:
        return False # Variant is an indel
    if no_unstable and var.startswith('('):
        return False # Recurrent/Unstable and we're filtering against that.
    if no_backmut and var.endswith('!'):
        return False # Back-mutation and we're filtering against that.
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
    for variants in reversed(var_stack):
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
        if level <= cur:
            # pop variants off.
        var_stack.append(variants)
        
    return var_pos, hap_tab


def main():
    """ Simple test of phylotree functions. """
    return 0

if __name__ == "__main__":
    sys.exit(main())
