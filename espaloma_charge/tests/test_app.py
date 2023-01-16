def test_charge():
    from rdkit import Chem; from espaloma_charge import charge
    molecule = Chem.MolFromSmiles("N#N")
    charge(molecule)

def test_charge_batch():
    from rdkit import Chem
    from espaloma_charge import charge
    molecule = Chem.MolFromSmiles("N#N")
    charge([molecule, molecule])

