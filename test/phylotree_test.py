"""
Unit Tests for Phylotree module.
"""

import unittest
import sys

import phylotree


class TestVariantMethods(unittest.TestCase):
    def test_is_snp(self):
        self.assertTrue(phylotree.is_snp("T14034C"))
        self.assertTrue(phylotree.is_snp("T182C!"))
        self.assertTrue(phylotree.is_snp("(T195C)"))
        self.assertTrue(phylotree.is_snp("(T195C!)"))
        self.assertTrue(phylotree.is_snp("C182T!!"))
        self.assertTrue(phylotree.is_snp("T10454c"))
        self.assertFalse(phylotree.is_snp("745.1T"))
        self.assertFalse(phylotree.is_snp("(745.1T)"))
        self.assertFalse(phylotree.is_snp("T16325d"))
        self.assertFalse(phylotree.is_snp("5899.XCd!"))
        self.assertFalse(phylotree.is_snp("(C965d)"))

    def test_der_allele(self):
        self.assertEqual(phylotree.der_allele("T14034C"), 'C')
        self.assertEqual(phylotree.der_allele("T182C!"), 'C')
        self.assertEqual(phylotree.der_allele("(T195C)"), 'C')
        self.assertEqual(phylotree.der_allele("(T195C!)"), 'C')
        self.assertEqual(phylotree.der_allele("C182T!!"), 'T')
        self.assertEqual(phylotree.der_allele("T10454c"), 'C')

    def test_anc_allele(self):
        self.assertEqual(phylotree.anc_allele("T14034C"), 'T')
        self.assertEqual(phylotree.anc_allele("T182C!"), 'T')
        self.assertEqual(phylotree.anc_allele("(T195C)"), 'T')
        self.assertEqual(phylotree.anc_allele("(T195C!)"), 'T')
        self.assertEqual(phylotree.anc_allele("C182T!!"), 'C')
        self.assertEqual(phylotree.anc_allele("T10454c"), 'T')

    def test_is_unstable(self):
        self.assertTrue(phylotree.is_unstable("(T195C)"))
        self.assertTrue(phylotree.is_unstable("(T195C!)"))
        self.assertTrue(phylotree.is_unstable("(745.1T)"))
        self.assertTrue(phylotree.is_unstable("(C965d)"))
        self.assertFalse(phylotree.is_unstable("C182T!!"))
        self.assertFalse(phylotree.is_unstable("T10454c"))
        self.assertFalse(phylotree.is_unstable("745.1T"))
        self.assertFalse(phylotree.is_unstable("T14034C"))
        self.assertFalse(phylotree.is_unstable("T182C!"))
        self.assertFalse(phylotree.is_unstable("T16325d"))
        self.assertFalse(phylotree.is_unstable("5899.XCd!"))

    def test_is_backmutation(self):
        self.assertTrue(phylotree.is_backmutation("T182C!"))
        self.assertTrue(phylotree.is_backmutation("(T195C!)"))
        self.assertTrue(phylotree.is_backmutation("C182T!!"))
        self.assertTrue(phylotree.is_backmutation("5899.XCd!"))
        self.assertFalse(phylotree.is_backmutation("(745.1T)"))
        self.assertFalse(phylotree.is_backmutation("(C965d)"))
        self.assertFalse(phylotree.is_backmutation("(T195C)"))
        self.assertFalse(phylotree.is_backmutation("T10454c"))
        self.assertFalse(phylotree.is_backmutation("745.1T"))
        self.assertFalse(phylotree.is_backmutation("T14034C"))
        self.assertFalse(phylotree.is_backmutation("T16325d"))

    def test_rm_snp_annot(self):
        self.assertEqual(phylotree.rm_snp_annot("T182C!"), "T182C")
        self.assertEqual(phylotree.rm_snp_annot("(T195C!)"), "T195C")
        self.assertEqual(phylotree.rm_snp_annot("C182T!!"), "C182T")
        self.assertEqual(phylotree.rm_snp_annot("(T195C)"), "T195C")
        self.assertEqual(phylotree.rm_snp_annot("T10454c"), "T10454C")
        self.assertEqual(phylotree.rm_snp_annot("T14034C"), "T14034C")

    def test_pos_from_var(self):
        self.assertEqual(phylotree.pos_from_var("T182C!"), 181)
        self.assertEqual(phylotree.pos_from_var("(T195C!)"), 194)
        self.assertEqual(phylotree.pos_from_var("C182T!!"), 181)
        self.assertEqual(phylotree.pos_from_var("(T195C)"), 194)
        self.assertEqual(phylotree.pos_from_var("T10454c"), 10453)
        self.assertEqual(phylotree.pos_from_var("T14034C"), 14033)


class TestPhylotreeSimple(unittest.TestCase):
    def setUp(self):
        phy_in = ['I, A1G ,,',
                  ',H, A3T A5T ,,',
                  ',,F, A6T ,,',
                  ',,,B, A8T ,,',
                  ',,,C, T5A! ,,',
                  ',,G, A7T ,,',
                  ',,,D, (A9T) ,,',
                  ',,,E, A4t ,,',
                  ',A, A2T a4t ,,']
        self.phy = phylotree.Phylotree(phy_in)
        self.ref = "AAAAAAAAA"

    def test_phylonode_get_anon_name(self):
        node = phylotree.Phylotree.PhyloNode('A')
        for i in xrange(1, 5):
            self.assertEqual('A[%d]' % i, node.get_anon_name())

    def test_get_variant_pos(self):
        self.assertEqual(self.phy.get_variant_pos(), [0,1,2,3,4,5,6,7,8])

    def test_add_custom_hap(self):
        self.assertNotIn('J', self.phy.hap_var)
        self.phy.add_custom_hap('J', ['A1G', 'A5C', 'A9T'])
        self.assertIn('J', self.phy.hap_var)
        self.assertEqual(self.phy.get_variant_pos(), [0,1,2,3,4,5,6,7,8])

    def test_add_custom_hap_same_name(self):
        with self.assertRaises(ValueError):
            self.phy.add_custom_hap('A', ['A1G', 'A5C', 'A9T'])

    def test_phylotree_nodes(self):
        self.assertEquals('I', self.phy.root.hap_id)

    def test_hapvar(self):
        self.assertIn('A', self.phy.hap_var)
        self.assertIn('B', self.phy.hap_var)
        self.assertIn('C', self.phy.hap_var)
        self.assertIn('D', self.phy.hap_var)
        self.assertIn('E', self.phy.hap_var)
        self.assertIn('F', self.phy.hap_var)
        self.assertIn('G', self.phy.hap_var)
        self.assertIn('H', self.phy.hap_var)
        self.assertIn('I', self.phy.hap_var)
        self.assertEqual(self.phy.hap_var['A'], ['A1G', 'A2T', 'A4T'])
        self.assertEqual(self.phy.hap_var['B'], ['A1G', 'A3T', 'A5T',
                                                'A6T', 'A8T'])
        self.assertEqual(self.phy.hap_var['C'], ['A1G', 'A3T', 'T5A', 'A6T'])
        self.assertEqual(self.phy.hap_var['D'], ['A1G', 'A3T', 'A5T',
                                                'A7T', 'A9T'])
        self.assertEqual(self.phy.hap_var['E'], ['A1G', 'A3T', 'A4T',
                                                'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['F'], ['A1G', 'A3T', 'A5T', 'A6T'])
        self.assertEqual(self.phy.hap_var['G'], ['A1G', 'A3T', 'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['H'], ['A1G', 'A3T', 'A5T'])
        self.assertEqual(self.phy.hap_var['I'], ['A1G'])

    def test_variants(self):
        for i in range(9):
            self.assertIn(i, self.phy.variants)
        self.assertNotIn(9, self.phy.variants, "Bad 1-based to 0-based conv.")
        self.assertEqual(self.phy.variants[0], {'G':1})
        self.assertEqual(self.phy.variants[3], {'T':2})
        self.assertEqual(self.phy.variants[4], {'T':1, 'A':1})
        for i in [1, 2, 5, 6, 7, 8]:
            self.assertEqual(self.phy.variants[i], {'T':1})

    def test_ignore_sites(self):
        self.phy.ignore_sites("1")
        self.assertNotIn(0, self.phy.variants)
        self.assertIn(0, self.phy.ignore)
        self.assertEqual(self.phy.hap_var['A'], ['A2T', 'A4T'])
        self.assertEqual(self.phy.hap_var['B'], ['A3T', 'A5T', 'A6T', 'A8T'])
        self.assertEqual(self.phy.hap_var['C'], ['A3T', 'T5A', 'A6T'])
        self.assertEqual(self.phy.hap_var['D'], ['A3T', 'A5T', 'A7T', 'A9T'])
        self.assertEqual(self.phy.hap_var['E'], ['A3T', 'A4T', 'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['F'], ['A3T', 'A5T', 'A6T'])
        self.assertEqual(self.phy.hap_var['G'], ['A3T', 'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['H'], ['A3T', 'A5T'])
        self.assertEqual(self.phy.hap_var['I'], [])

    def test_ignore_sites_range(self):
        self.phy.ignore_sites("1-3")
        for i in xrange(3):
            self.assertNotIn(i, self.phy.variants)
            self.assertIn(i, self.phy.ignore)
        self.assertEqual(self.phy.hap_var['A'], ['A4T'])
        self.assertEqual(self.phy.hap_var['B'], ['A5T', 'A6T', 'A8T'])
        self.assertEqual(self.phy.hap_var['C'], ['T5A', 'A6T'])
        self.assertEqual(self.phy.hap_var['D'], ['A5T', 'A7T', 'A9T'])
        self.assertEqual(self.phy.hap_var['E'], ['A4T', 'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['F'], ['A5T', 'A6T'])
        self.assertEqual(self.phy.hap_var['G'], ['A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['H'], ['A5T'])
        self.assertEqual(self.phy.hap_var['I'], [])

    def test_ignore_sites_mix(self):
        self.phy.ignore_sites("1,2-3")
        for i in xrange(3):
            self.assertNotIn(i, self.phy.variants)
            self.assertIn(i, self.phy.ignore)
        self.assertEqual(self.phy.hap_var['A'], ['A4T'])
        self.assertEqual(self.phy.hap_var['B'], ['A5T', 'A6T', 'A8T'])
        self.assertEqual(self.phy.hap_var['C'], ['T5A', 'A6T'])
        self.assertEqual(self.phy.hap_var['D'], ['A5T', 'A7T', 'A9T'])
        self.assertEqual(self.phy.hap_var['E'], ['A4T', 'A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['F'], ['A5T', 'A6T'])
        self.assertEqual(self.phy.hap_var['G'], ['A5T', 'A7T'])
        self.assertEqual(self.phy.hap_var['H'], ['A5T'])
        self.assertEqual(self.phy.hap_var['I'], [])

    def test_find_polymorphic_one_haplogroup(self):
        self.assertEqual(self.phy.polymorphic_sites(['A'], self.ref), [])

    def test_find_polymorphic_two_haplogroups(self):
        res = self.phy.polymorphic_sites(['F', 'G'], self.ref)
        exp = [5, 6]
        self.assertEqual(res, exp)

    def test_find_polymorphic_two_haplogroups_backmutation(self):
        res = self.phy.polymorphic_sites(['A', 'C'], self.ref)
        exp = [1, 2, 3, 5]
        self.assertEqual(res, exp)

    def test_find_polymorphic_multi_haplogroups(self):
        res = self.phy.polymorphic_sites(['B', 'C', 'D', 'E'], self.ref)
        exp = [3, 4, 5, 6, 7, 8]
        self.assertEqual(res, exp)


if __name__ == '__main__':
    unittest.main()
