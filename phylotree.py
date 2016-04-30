"""

This file contains functions for reading in phylotree variants and haplogroups
and building a matrix to represent these markers.

Mon Apr  4 12:00:54 PDT 2016

"""

import sys
import collections


class Phylotree(object):
    """
    Class for representing haplogroups and the variants that define them in
    an explicit tree data structure. The general workflow of this is to load
    in the tree from a file based on the Phylotree mtDNA tree (phylotree.org),
    perform any filtering steps and produce a table of haplogroup IDs and
    the variants that define them.
    """
    class PhyloNode(object):
        """
        Class that represents a single node in the phylotree.
        """
        def __init__(self, hap='', parent=None, variants=None):
            """ Initialize new PhyloNode """
            self.hap_id = hap
            self.parent = parent
            self.children = list()
            self.anon_child = 0
            self.anon = (self.hap_id == '')
            if self.parent is not None:
                parent.children.append(self)
                if self.anon:
                    self.hap_id = parent.get_anon_name()
            if variants is None:
                self.variants = list()
            else:
                self.variants = variants

        def all_variants(self):
            """
            Produces a list of variants that define this haplogroup including
            variants inherited from the node's parent. Importantly, mutations
            that occur to the same position are masked by the most recent
            mutation (i.e. C152T will not be included if T152C occurs farther
            down in the tree).
            """
            summed_vars = dict()
            node = self
            while node is not None:
                for var in node.variants:
                    pos = pos_from_var(var)
                    if pos not in summed_vars:
                        summed_vars[pos] = var
                node = node.parent
            return [summed_vars[pos] for pos in sorted(summed_vars)]

        def dump(self, out=sys.stderr, indent=0):
            """
            Dumps the contents of this node and its children recursively.
            """
            prefix = ' ' * indent
            out.write('%sNode: %s\n' % (prefix, self.hap_id))
            out.write('%sVariants: %s\n' % (prefix, ','.join(self.variants)))
            out.write('%sChildren: %d\n' % (prefix, len(self.children)))
            for child in self.children:
                child.dump(out, indent + 2)

        def get_anon_name(self):
            """ 
            Returns a unique haplogroup name based on this one for a child with
            no specified name.
            """
            self.anon_child += 1
            return "%s[%d]" % (self.hap_id, self.anon_child)

    def __init__(self, phy_in=None, rm_unstable=False, rm_backmut=False):
        """ Initialize a blank Phylotree before reading from a file. """
        self.root  = None
        self.nodes = list()
        self.variants = collections.defaultdict(collections.Counter)
        self.ignore = None
        self.hap_var  = None
        if phy_in is not None:
            self.read_csv(phy_in)
            self.process_variants(rm_unstable, rm_backmut)
            self.process_haplotypes()

    def read_csv(self, phy_in):
        """ Reads tree structure and variants from input stream. """
        cur = -1
        node_stack = list([None])
        for line in phy_in:
            level, hap_id, raw_var = _read_phy_line(line)
            variants = [var for var in raw_var if is_snp(var)]
            while cur >= 0 and cur >= level:
                node_stack.pop()
                cur -= 1
            new_node = Phylotree.PhyloNode(hap_id, node_stack[-1], variants)
            self.nodes.append(new_node)
            node_stack.append(new_node)
            cur += 1
        else:
            # Set the root for the tree. First item on node_stack is None.
            self.root = node_stack[1]
        return

    def process_variants(self, rm_unstable, rm_backmut):
        """
        Builds a table of variant sites stored in Phylotree.variants after
        filtering based on phylotree annotation. This table keeps track of
        number of mutations and derived alleles as they occur.
        """
        var_pos = set()
        self.ignore = set()
        for node in self.nodes:
            for var in node.variants:
                pos = pos_from_var(var)
                if rm_unstable and is_unstable(var):
                    self.ignore.add(pos)
                elif rm_backmut and is_backmutation(var):
                    self.ignore.add(pos)
                else:
                    var_pos.add(pos)
                    der = der_allele(var)
                    self.variants[pos][der] += 1
        if rm_unstable or rm_backmut:
            var_pos -= self.ignore
            for pos in self.ignore:
                if pos in self.variants:
                    del self.variants[pos]
            for node in self.nodes:
                node.variants = [var for var in node.variants
                                 if pos_from_var(var) not in self.ignore]
        return

    def process_haplotypes(self):
        """
        Contructs a dictionary that maps haplogroup IDs to a list of
        variants associated with that haplogroup, using the Phylotree
        information and variants after any filtering. Merges haplogroups with
        identical variants.
        """
        haplotypes = collections.defaultdict(list)
        hap_tab = dict()
        for node in self.nodes:
            var_str = ','.join(node.all_variants())
            haplotypes[var_str].append(node)
        for var_str in haplotypes:
            variants = var_str.split(',')
            haplo_id = '/'.join([node.hap_id for node in haplotypes[var_str]])
            hap_tab[haplo_id] = variants
        self.hap_var = hap_tab
        return

    def get_variant_pos(self):
        """
        Returns a list of positions with variants after filtering.
        """
        return sorted(self.variants.keys())

    def add_custom_hap(self, hap_id, variants):
        """
        Adds a custom haplotype to the haplotype variant table. Does not add a
        node in the tree for this haplotype. Ignores positions that have been
        removed from consideration by process_variants().
        """
        if hap_id in self.hap_var:
            raise ValueError("Custom haplogroup name '%d' already in use."
                             % (hap_id))
        ok_vars = [var for var in variants 
                   if pos_from_var(var) not in self.ignore]
        self.hap_var[hap_id] = ok_vars
        return


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


def main():
    """ Simple test of phylotree functions. """
    if len(sys.argv) > 1:
        phy_fn = sys.argv[1]
        with open(phy_fn, 'r') as phy_in:
            phy = Phylotree(phy_in)
            #phy.root.dump(sys.stdout)
            #phy.process_variants(rm_unstable=False, rm_backmut=False)
            for pos in sorted(phy.variants):
                print pos, phy.variants[pos], (sum(phy.variants[pos].values()))
#           for hap in phy.hap_var:
#               print hap, ','.join(phy.hap_var[hap])
    return 0


if __name__ == "__main__":
    sys.exit(main())
