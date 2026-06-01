import unittest
import pandas as pd
from virtual_screening import calculate_similarity

class TestVirtualScreening(unittest.TestCase):
    def setUp(self):
        self.target = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
        self.library = [
            {'id': 'CHEMBL941', 'smiles': 'CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5'}, # Imatinib (itself)
            {'id': 'CHEMBL1201138', 'smiles': 'CC1=CN=C(C=N1)C2=C(C=C(C=C2)C)NC3=NC=CC(=N3)C4=CN=CC=C4'}, # Nilotinib
            {'id': 'INVALID', 'smiles': 'NOT_A_SMILES'}, # Invalid
            {'id': 'CHEMBL192', 'smiles': 'CCN(CC)CCNC(=O)C1=CC=C(C=C1)N'} # Low similarity
        ]

    def test_similarity_calculation(self):
        df = calculate_similarity(self.target, self.library)
        # Should contain Imatinib (1.0) and Nilotinib (high similarity)
        self.assertIn('CHEMBL941', df['ChEMBL_ID'].values)
        self.assertIn('CHEMBL1201138', df['ChEMBL_ID'].values)
        # Should not contain INVALID or low similarity molecules (depending on threshold)
        self.assertNotIn('INVALID', df['ChEMBL_ID'].values)
        self.assertTrue(df['Similarity'].iloc[0] == 1.0)

    def test_invalid_target(self):
        with self.assertRaises(ValueError):
            calculate_similarity("INVALID_TARGET", self.library)

if __name__ == '__main__':
    unittest.main()
