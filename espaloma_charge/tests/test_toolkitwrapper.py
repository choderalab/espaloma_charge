from espaloma_charge.espaloma_wrapper import EspalomaChargeToolkitWrapper


def test_ethanol():
    from openff.toolkit import Molecule
    offmol = Molecule.from_smiles('CCO')
    offmol.assign_partial_charges('espaloma-am1bcc',
                                  toolkit_registry=EspalomaChargeToolkitWrapper())
    print(offmol.partial_charges)
