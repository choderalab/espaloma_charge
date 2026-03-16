def test_charge():
    from rdkit import Chem; from espaloma_charge import charge
    molecule = Chem.MolFromSmiles("N#N")
    charge(molecule)

def test_charge_batch():
    from rdkit import Chem
    from espaloma_charge import charge
    molecule = Chem.MolFromSmiles("N#N")
    charge([molecule, molecule])

def test_total_charge_constraint():
    from rdkit import Chem
    from espaloma_charge import charge
    import numpy as np

    """Test that the total_charge argument is respected for a charged molecule."""
    mol = Chem.MolFromSmiles("[NH4+]")
    net_charge = 1.0
    q = charge(mol, total_charge=net_charge)
    
    assert np.isclose(np.sum(q), net_charge, atol=1e-4), \
        f"Expected total charge {net_charge}, but got {np.sum(q):.4f}"