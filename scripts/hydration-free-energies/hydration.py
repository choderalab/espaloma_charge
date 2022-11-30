#!/bin/env python
"""
Compute hydration free energy for a neutral compound using explicit solvent alchemical replica exchange free energy calculations.

"""

import click

@click.group()
def cli():
    pass

@click.command()
@click.option('--smiles', 
              required=True,
              help='SMILES string of molecule whose hydration free energy is to be computed')
@click.option('--toolkit',
              default='EspalomaCharge',
              type=click.Choice(['EspalomaCharge', 'AmberTools', 'OpenEye', 'RDKit']),
              help='Toolkit to use for assigning charges.')              
@click.option('--method', 
              default='espaloma-am1bcc',
              help='The charge model to use from the toolkit.')
def run(smiles, toolkit='', method='espaloma-am1bcc'):
    """
    Compute the hydration free energy of a specified molecule

    Parameters
    ----------
    smiles : str
        The SMILES string of the compound whose hydration free energy is to be computed.
        The SMILES string can have explicit or implicit protons.
        The compound must be neutral.
    toolkit : str, optional, default='EspalomaCharge'
        The toolkit to use for assigning charges.
        Valid options are ['EspalomaCharge', 'AmberTools', 'OpenEye', 'RDKit']
        'ToolkitWrapper' is appended to the toolkit name.
    method : str, optional, default='espaloma-am1bcc1'
        The method to use for assigning partial charges.
        Valid options depend on the toolkit:
            'espaloma' : ['espaloma-am1bcc']
            'AmberTools' : ['am1bcc', 'am1-mulliken', 'gasteiger']
            'OpenEye' : ['am1bcc', 'am1-mulliken', 'gasteiger', 'mmff94', 'am1bccnosymspt', 'am1elf10', 'am1bccelf10']
            'RDKit' : ['mmff94'] 
    """
    # Create an OpenFF Molecule object
    from openff.toolkit.topology import Molecule
    molecule = Molecule.from_smiles(smiles)

    # Check to make sure molecule is neutral
    if molecule.total_charge.magnitude != 0.0:
        raise ValueError(f'Molecule {smiles} has a net charge of {molecule.total_charge}. Only neutral molecules are supported.')

    # Assign partial charges
    toolkit_wrapper_name = toolkit + 'ToolkitWrapper'
    if toolkit_wrapper_name == 'EspalomaChargeToolkitWrapper':
        # TODO: Eliminate this branch if we integrate this wrapper into the OpenFF Toolkit
        from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        toolkit_wrapper = EspalomaChargeToolkitWrapper()
    else:
        import openff.toolkit.utils.toolkits
        toolkit_wrapper = getattr(openff.toolkit.utils.toolkits, toolkit_wrapper_name)()

    print(toolkit_wrapper)
    molecule.assign_partial_charges(method, toolkit_registry=toolkit_wrapper)

    # Run vacuum phase
    print(f'Assigned partial charges from toolkit {toolkit} : {method}')
    print(molecule.partial_charges)

cli.add_command(run)

if __name__ == '__main__':
    cli()