import unittest

from format_plasmid_headers import scaffold_type


class ScaffoldType(unittest.TestCase):
    def test_scaffold_titles(self):
        cases = {
            'JACFBF010000002.1 Borreliella bavariensis strain A91S A91S_cp26_NODE_10_length_27256, whole genome': 'cp26',
            'JACFBF010000001.1 Borreliella bavariensis strain A91S A91S_chr_to_PBi_insert, whole genome': 'chromosome',
            'JACFBF010000006.1 Borreliella bavariensis strain A91S A91S_cp32-3lp25_to_PBi, whole genome': 'cp32-3-lp25',
            'JACFBF010000008.1 Borreliella bavariensis strain A91S A91S_lp28-4cp32-1_NODE_4, whole genome': 'lp28-4-cp32-1',
            'X.1 Borreliella bavariensis strain Y Y_lp28-4_cp32-1_NODE_4, whole genome': 'lp28-4-cp32-1',
            'JACFAY010000005.1 Borreliella bavariensis strain PBar PBar_lp17_mapped_to_PNeb, whole genome': 'lp17',
        }
        for header, name in cases.items():
            self.assertEqual(scaffold_type(header), name, header)

    def test_unnamed_scaffold_is_skipped(self):
        header = 'JACFBF010000014.1 Borreliella bavariensis strain A91S NODE_20_length_1563_cov_10.6121, whole genome'
        self.assertIsNone(scaffold_type(header))


if __name__ == '__main__':
    unittest.main()
