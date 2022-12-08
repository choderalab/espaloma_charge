#!/bin/env python
"""
Compute hydration free energy for a neutral compound using explicit solvent alchemical replica exchange free energy calculations.

"""

import click

@click.group()
def cli():
    pass

@click.command()
@click.option('--model-url', 
              default=None,
              help='If specified, provide model URL or local filepath to espaloma_charge')
def freesolv(model_url):
    """Run specified molecule index from FreeSolv database
    """
    # Load FreeSolv database
    # TODO: Retrieve automatically if needed
    import json
    with open('freesolv.json', 'rt') as infile:
        freesolv = json.load(infile)

    # Extract info
    print('Creating molecules...')
    molecules = dict()
    for name, entry in freesolv.items():
        smiles = freesolv[name]['smiles']
        from openff.toolkit.topology import Molecule
        from openff.toolkit.utils.exceptions import UndefinedStereochemistryError
        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        molecules[name] = molecule
    print(f'{len(molecules)} molecules read')

    # DEBUG
    molecules = { name : molecules[name] for name in list(molecules.keys())[:10] }
    print(f'{len(molecules)} molecules remaining')

    # Assign charges
    toolkit_methods = {
        'EspalomaCharge' : 'espaloma-am1bcc',
        'AmberTools' : 'am1bcc',
        'OpenEye' : 'am1bcc', 
    }
    toolkits = list(toolkit_methods.keys())
    charged_molecules = dict()
    import time
    for toolkit, method in toolkit_methods.items():
        print(f'Assigning charges for {toolkit}')
        import copy
        charged_molecules[toolkit] = copy.deepcopy(molecules)
        initial_time = time.time()
        assign_charges(charged_molecules[toolkit].values(), toolkit, method, model_url=model_url)
        elapsed_time = time.time() - initial_time
        print(f'Assign {len(molecules)} charges with {toolkit} took {elapsed_time:.3f} s ({elapsed_time/len(molecules):.3f} s / molecule)')

    # Compare charges
    print(f'Charge model comparison RMSE on {len(molecules)} from FreeSolv')
    import itertools
    for toolkit1, toolkit2 in itertools.combinations(toolkits, 2):
        rmse = compute_rmse(charged_molecules[toolkit1].values(), charged_molecules[toolkit2].values())
        print(f'{toolkit1:25s} {toolkit2:25s} {rmse:8.3f}')


def compute_rmse(molecules1, molecules2):
    """Compute RMSE between sets of charges
    
    Parameters
    ----------
    molecules1 : list of Molecules
    molecules2 : list of Molecules
        Molecule sets to compare

    Return
    ------
    rmse : float
        RMSE between sets of charges
    """
    import numpy as np
    mse = 0.0
    assert len(molecules1) == len(molecules2)
    nmolecules = len(molecules1)
    for molecule1, molecule2 in zip(molecules1, molecules2):
        assert molecule1.name == molecule2.name
        mse += np.sum((molecule1.partial_charges - molecule2.partial_charges)**2)
    mse /= nmolecules
    rmse = np.sqrt(mse) 
    return rmse   

def assign_charges(molecules, toolkit, method, model_url=None):
    """Assign charges for all molecules using the specified toolkit and method.

    Parameters
    ----------
    molecules : list of openff.toolkit.topology.Molecule
        The list of molecules to assign charges to.
        The molecules will be modified in place.
    toolkit : str
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
    model_url : str
        Provide optional model_url or filepath to espaloma_charge.app.charge()       
    """
    # Get toolkit
    toolkit_wrapper_name = toolkit + 'ToolkitWrapper'
    if toolkit_wrapper_name == 'EspalomaChargeToolkitWrapper':
        # TODO: Eliminate this branch if we integrate this wrapper into the OpenFF Toolkit
        from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
        toolkit_wrapper = EspalomaChargeToolkitWrapper()
        # TODO: We would still ened to keep this
        if model_url is not None:
            print(f'Setting toolkit_wrapper.model_url = {model_url}')
            toolkit_wrapper.model_url = model_url
            
    else:
        import openff.toolkit.utils.toolkits
        toolkit_wrapper = getattr(openff.toolkit.utils.toolkits, toolkit_wrapper_name)()

    # Assign charges to all molecules
    for molecule in molecules:
        molecule.assign_partial_charges(method, toolkit_registry=toolkit_wrapper)

cli.add_command(freesolv)

if __name__ == '__main__':
    cli()