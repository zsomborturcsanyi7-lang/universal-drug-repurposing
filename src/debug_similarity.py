from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

target_smiles = "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
nilotinib_smiles = "CC1=CN=C(C=N1)C2=C(C=C(C=C2)C)NC3=NC=CC(=N3)C4=CN=CC=C4"

mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
target_fp = mfpgen.GetFingerprint(Chem.MolFromSmiles(target_smiles))
nilotinib_fp = mfpgen.GetFingerprint(Chem.MolFromSmiles(nilotinib_smiles))

similarity = DataStructs.TanimotoSimilarity(target_fp, nilotinib_fp)
print(f"Similarity: {similarity}")
