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

    Atributes:
        refseq: The reference sequence upon which this tree is bases (RSRS)
        root: The root node for the internal tree representation
        nodes: tree nodes in a list for easy iteration
        variants: A dictionary of counters tracking the mutations that have
                  occurred according to phylotree for each variant position.
        ignore: A set of reference positions that should be ignored
        hap_var: A dictionary mapping haplogroup ID strings to lists of
                 variants strings associated with each haplogroup
        anon_haps: Flag for whether anonymous haplogroups should be included
                   in the table of haplogroups and variants (hap_var).
        rm_unstable: Flag for whether unstable sites are ignored
        rm_backmut: Flag for whether sites with backmutations are ignored.
    """
    class PhyloNode(object):
        """
        Class that represents a single node in the phylotree.

        Attributes:
            hap_id: Haplogroup ID for this node
            parent: The parent node of this node (None if this is the root)
            children: List of child nodes of this node
            anon_child: Int counting number of anonymous children so far.
            anon: True if this node did not have a Haplogroup ID.
        """
        def __init__(self, hap='', parent=None, variants=None):
            """
            Initialize new PhyloNode

            Args:
                hap: Haplogroup ID (str)
                parent: The parent node for this new node.
                variants: a list of variant string for this node.
            Returns: nothing
            """
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

            Args: (self)
            Returns:
                a list of variant strings, sorted by position
            """
            summed_vars = dict()
            node = self
            while node is not None:
                for var in node.variants:
                    pos = pos_from_var(var)
                    if pos not in summed_vars:
                        summed_vars[pos] = rm_snp_annot(var)
                node = node.parent
            return [summed_vars[pos] for pos in sorted(summed_vars)]

        def dump(self, out=sys.stderr, indent=0):
            """
            Dumps the contents of this node and its children recursively.

            Args:
                out: the destination output stream
                indent: The number of space to indent before writing to
                        indicate depth in the tree.
            Returns: nothing
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

    def __init__(self, phy_in=None, refseq=None, anon_haps=False,
                 rm_unstable=False, rm_backmut=False):
        """
        Initialize a blank Phylotree before reading from a file.

        Args:
            phy_in: input stream of Phylotree in simplified CSV format.
            rm_unstable: Ignore positions that have unstable variants
            rm_backmut: Ignore positions that have backmutations
        Returns: nothing
        """
        self.root = None
        self.nodes = list()
        self.variants = collections.defaultdict(collections.Counter)
        self.ignore = set()
        self.hap_var = None
        self.refseq = refseq
        self.anon_haps = anon_haps
        self.rm_unstable = rm_unstable
        self.rm_backmut = rm_backmut
        if phy_in is not None:
            self.read_csv(phy_in)
            self.process_variants()
            self.process_haplotypes()

    def read_csv(self, phy_in):
        """
        Reads tree structure and variants from input stream. Builds an internal
        representation of the tree.

        Args:
            phy_in: input stream of Phylotree in simplified CSV format.
        Returns: nothing
        """
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
        # Set the root for the tree. First item on node_stack is None.
        self.root = node_stack[1]
        return

    def process_variants(self):
        """
        Builds a table of variant sites stored in Phylotree.variants after
        filtering based on phylotree annotation. This table keeps track of
        number of mutations and derived alleles as they occur.

        Args: (self)
        Returns: nothing
        """
        var_pos = set()
        for node in self.nodes:
            for var in node.variants:
                pos = pos_from_var(var)
                if self.rm_unstable and is_unstable(var):
                    self.ignore.add(pos)
                elif self.rm_backmut and is_backmutation(var):
                    self.ignore.add(pos)
                else:
                    var_pos.add(pos)
                    der = der_allele(var)
                    self.variants[pos][der] += 1
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

        Args: (self)
        Returns: nothing
        """
        haplotypes = collections.defaultdict(list)
        hap_tab = dict()
        for node in self.nodes:
            if self.anon_haps or not node.anon:
                var_str = ','.join(node.all_variants())
                haplotypes[var_str].append(node)
        for var_str in haplotypes:
            variants = [var for var in var_str.split(',') if var != '']
            haplo_id = '/'.join([node.hap_id for node in haplotypes[var_str]])
            hap_tab[haplo_id] = variants
        self.hap_var = hap_tab
        return

    def get_variant_pos(self):
        """
        Get a sorted list of all the variant positions this Phyotree object
        keeps track of.

        Args: (self)
        Returns: a list of positions with variants after filtering.
        """
        return sorted(self.variants.keys())

    def add_custom_hap(self, hap_id, variants):
        """
        Adds a custom haplotype to the haplotype variant table. Does not add a
        node in the tree for this haplotype. Ignores positions that have been
        removed from consideration by process_variants().

        Args:
            hap_id: Haplogroup ID for this custom entry
            variants: A list of variant strings
        Returns:
            nothing
        Raises:
            ValueError: if haplogroup name is already in use.
        """
        if hap_id in self.hap_var:
            raise ValueError("Custom haplogroup name '%s' already in use."
                             % (hap_id))
        ok_vars = [var for var in variants
                   if pos_from_var(var) not in self.ignore]
        self.hap_var[hap_id] = ok_vars
        return

    def ignore_sites(self, sites_str):
        """
        Adds the sites described by the argument into the ignore list and
        rebuilds the variant and haplotype tables. site_str is a string that
        represents a comma-separated list of 1-based positions or 'start-end'
        ranges (inclusive).

        Args:
            site_str: a string representing positions and ranges of positions
                      in the reference in 1-based coords.
        Returns: nothing
        """
        sites = sites_str.split(',')
        for site in sites:
            if '-' in site:
                start, end = site.split('-')
                self.ignore.update(list(range(int(start) - 1, int(end))))
            else:
                self.ignore.add(int(site) - 1)
        self.process_variants()
        self.process_haplotypes()
        return

    def polymorphic_sites(self, haps, ref=None):
        """
        Takes a list of haplogroups and returns the sites that are expected to
        be polymorphic within the sample.

        Args:
            self: this Phylotree object
            haps: a list of haplogroup IDs
            ref: the reference sequence to use when checking for backmutations
        Returns:
            a list of positions in the reference.
        """
        def derived_diff(variants):
            """
            Takes a list of variants and returns True if they all result in the
            same derived base.

            Args:
                variants: A list of variant strings
            Returns:
                True if all result in the same derived base, False otherwise
            """
            der_bases = [der_allele(var) for var in variants]
            return not all([base == der_bases[0] for base in der_bases])

        if ref is None:
            ref = self.refseq
        var_tab = collections.defaultdict(list)
        for hap in haps:
            for var in self.hap_var[hap]:
                pos = pos_from_var(var)
                der = der_allele(var)
                if der != ref[pos]:
                    var_tab[pos].append(var)
        # Position is polymorphic if not all contributors have a non reference
        # allele or if any derived alleles differ.
        poly = [pos for pos in sorted(var_tab)
                if (len(var_tab[pos]) < len(haps)
                    or derived_diff(var_tab[pos]))]
        return poly

    def get_ancestral(self, hap_id):
        """
        Find the variant positions that have not been affected by a mutation
        in this haplogroup and return a list of positions (int, 0-based) and
        base tuples.

        Args:
            self: this Phylotree object
            hap_id: a haplogroup ID
        Returns:
            a list of (position (0-based), base) tuples
        """
        ancestral_bases = {pos:self.refseq[pos]
                           for pos in self.get_variant_pos()}
        for var in self.hap_var[hap_id]:
            pos = pos_from_var(var)
            if pos in ancestral_bases:
                del ancestral_bases[pos]
        return list(ancestral_bases.items())


def pos_from_var(var):
    """
    Extracts the position of the SNP variant, converted to 0-base.

    Args:
        var: A variant string representing a SNP
    Returns:
        The position of the variant as an int, converted from 1- to 0-base
    """
    if var.startswith('('):
        var = var[1:-1]
    var = var.rstrip('!')
    return int(var[1:-1]) - 1 # 1-based to 0-based


def der_allele(var):
    """
    Returns the derived allele of this SNP variant.

    Args:
        var: A variant string representing a SNP
    Returns:
        The derived base as uppercase
    """
    var = var.rstrip(')!')
    return var[-1].upper()


def anc_allele(var):
    """
    Returns the ancestral allele of this SNP variant.

    Args:
        var: A variant string representing a SNP
    Returns:
        The ancestral base as uppercase
    """
    var = var.lstrip('(')
    return var[0].upper()


def is_snp(var):
    """
    Check if var is a SNP (not an indel)

    Args:
        var: A variant string
    Returns:
        True if var represents a SNP, False otherwise.
    """
    if '.' in var or 'd' in var:
        return False # Variant is an indel
    return True


def is_unstable(var):
    """
    Check if SNP var is annotated as unstable.

    Args:
        var: A variant string representing a SNP
    Returns:
        True if marked as unstable, False otherwise.
    """
    if var[0] == '(':
        return True
    return False


def is_backmutation(var):
    """
    Check if SNP var is annotated as a backmutation.

    Args:
        var: A variant string representing a SNP
    Returns:
        True if marked as a backmutation, False otherwise.
    """
    if '!' in var:
        return True
    return False


def rm_snp_annot(var):
    """
    Returns the SNP variant string, nicely formatted with annotation stripped.
    This includes ()'s, !, and lowercase bases.

    Args:
        var: A variant string representing a SNP
    Returns:
        The variant string with annotation removed.
    """
    if var.startswith('('):
        var = var[1:-1]
    var = var.rstrip('!')
    return var.upper()


def example():
    """
    Returns an example tree we can use for testing.
    """
    #            I
    #           / \
    #          /   H
    #         /   / \
    #        /   /   \
    #       /   F     G
    #      /   / \   / \
    #     A   B   C D   E
    phy_in = ['I, A1G ,,',
              ',H, A3T A5T ,,',
              ',,F, A6T ,,',
              ',,,B, A8T ,,',
              ',,,C, T5A ,,',
              ',,G, A7T ,,',
              ',,,D, A9T ,,',
              ',,,E, A4T ,,',
              ',A, A2T A4T ,,']
    return Phylotree(phy_in)


def _read_phy_line(line):
    """
    Reads a single comma-separated line from phylotree and returns the
    indentation level, the haplogroup id, and variants. Importantly, some nodes
    do not have IDs, so we need to check that we do not accidentally grab the
    variant list instead of the blank id.

    Args:
        line: A string representing 1 line in CSV phylotree input
    Returns:
        level: The level of the node this line represents in the tree.
               (the number of blank fields before the haplogroup ID)
        hap_id: The haplogroup ID for this entry (can be '')
        variants: List of variant strings for this node.
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
            phy = Phylotree(phy_in, anon_haps=True)
            for hap in phy.hap_var:
                print(hap, ','.join(phy.hap_var[hap]))
    else:
        phy = example()
        hap_var = dict({'A':['A1G', 'A2T', 'A4T'],
                        'B':['A1G', 'A3T', 'A5T', 'A6T', 'A8T'],
                        'C':['A1G', 'A3T', 'T5A', 'A6T'],
                        'D':['A1G', 'A3T', 'A5T', 'A7T', 'A9T'],
                        'E':['A1G', 'A3T', 'A4T', 'A5T', 'A7T'],
                        'F':['A1G', 'A3T', 'A5T', 'A6T'],
                        'G':['A1G', 'A3T', 'A5T', 'A7T'],
                        'H':['A1G', 'A3T', 'A5T'],
                        'I':['A1G']})
        for hap in sorted(phy.hap_var):
            print(hap, phy.hap_var[hap], hap_var[hap])
        print(phy.variants)
        phy.root.dump()
    return 0


if __name__ == "__main__":
    sys.exit(main())
