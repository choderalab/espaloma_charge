from espaloma_charge.espaloma_wrapper import EspalomaChargeToolkitWrapper
from openff.toolkit import Molecule


def test_ethanol_molecule_api_w_registry():
    offmol = Molecule.from_smiles('CCO')
    offmol.assign_partial_charges('espaloma-am1bcc',
                                  toolkit_registry=EspalomaChargeToolkitWrapper())
    assert offmol.partial_charges != offmol.partial_charges * 0., offmol.partial_charges

def test_ethanol_molecule_api_wo_registry():
    offmol = Molecule.from_smiles('CCO')
    offmol.assign_partial_charges('espaloma-am1bcc')
    assert offmol.partial_charges != offmol.partial_charges * 0., offmol.partial_charges

def test_ethanol_direct():
    ectkw = EspalomaChargeToolkitWrapper()
    offmol = Molecule.from_smiles('CCO')
    ectkw.assign_partial_charges(offmol, 'espaloma-am1bcc')
    assert offmol.partial_charges != offmol.partial_charges * 0., offmol.partial_charges
