import unittest

from format_plasmid_headers import correct_pf32_header, existing_pf32_label, scaffold_type


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


class Pf32Corrections(unittest.TestCase):
    def test_rename_changes_the_id_the_caller_reads(self):
        corrections = {'RS01295_JD1_JD1_cp32-1+5_ParA_X': 'RS01295_JD1_JD1_cp32-5_ParA_X'}
        applied = set()
        header = correct_pf32_header('RS01295_JD1_JD1_cp32-1+5_ParA_X', corrections, applied)
        self.assertEqual(header, 'RS01295_JD1_JD1_cp32-5_ParA_X')
        self.assertEqual(existing_pf32_label(header), 'cp32-5')
        self.assertEqual(applied, set(corrections))

    def test_other_records_pass_through(self):
        applied = set()
        self.assertEqual(correct_pf32_header('X_lp54_ParA extra', {'Y': 'Z'}, applied), 'X_lp54_ParA extra')
        self.assertEqual(applied, set())


if __name__ == '__main__':
    unittest.main()
