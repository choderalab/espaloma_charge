import pytest
# import torch

def test_cooh():
    from rdkit import Chem
    from espaloma_charge.utils import from_rdkit_mol
    import torch
    mol = Chem.MolFromSmiles("[O-]C=O")
    mol = Chem.AddHs(mol)
    g = from_rdkit_mol(mol)
    assert torch.unique(g.ndata["h0"][g.ndata["type"] == 8], dim=0).shape[0] == 1
 
def test_benzamidine():
    from rdkit import Chem
    from espaloma_charge.utils import from_rdkit_mol
    import torch
    mol = Chem.MolFromSmiles("C1=CC=C(C=C1)C(=[NH2+])N")
    mol = Chem.AddHs(mol)
    g = from_rdkit_mol(mol)
    assert torch.unique(g.ndata["h0"][g.ndata["type"] == 7], dim=0).shape[0] == 1

    mol = Chem.MolFromSmiles("C1=CC=C(C=C1)C(=[N])N")
    mol = Chem.AddHs(mol)
    g = from_rdkit_mol(mol)
    assert torch.unique(g.ndata["h0"][g.ndata["type"] == 7], dim=0).shape[0] == 2
        
def test_guanidinium():
    import torch
    from rdkit import Chem
    from espaloma_charge.utils import from_rdkit_mol
    mol = Chem.MolFromSmiles("C(=[NH2+])(N)N")
    mol = Chem.AddHs(mol)
    g = from_rdkit_mol(mol)
    assert torch.unique(g.ndata["h0"][g.ndata["type"] == 7], dim=0).shape[0] == 1
 


