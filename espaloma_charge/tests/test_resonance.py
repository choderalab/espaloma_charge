import pytest

def test_cooh():
    from rdkit import Chem, AddHs
    from espaloma_charge.utils import from_rdkit_mol
    mol = Chem.MolFromSmiles("[O-]C=O")
    mol = Chem.AddHs(mol)
    g = from_rdkit_mol(mol)
    print(g)
