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
    #molecules = { name : molecules[name] for name in list(molecules.keys())[:5] }
    #print(f'{len(molecules)} molecules remaining')

    # Assign charges
    toolkit_methods = {
        'EspalomaCharge' : 'espaloma-am1bcc',
        'OpenEye' : 'am1bcc', 
        'AmberTools' : 'am1bcc',
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
        # Compute RMSE
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
        # Warning: Using .magnitude is only safe because we know both are in units of 'elementary_charge'
        mse += np.sum((molecule1.partial_charges.magnitude - molecule2.partial_charges.magnitude)**2)
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
    from rich.progress import track
    for molecule in track(molecules, description='Assigning partial charges to molecules...'):
        molecule.assign_partial_charges(method, toolkit_registry=toolkit_wrapper)


@click.command()
@click.option('--model-url', 
              default=None,
              help='If specified, provide model URL or local filepath to espaloma_charge')
def errors(model_url):
    """Write out molecules with largest hydration free energy discrepancies between espaloma and OpenEye am1bccelf10
    also showing charge RMSE.
    """
    # Load FreeSolv data for all pairs of methods.
    print('Loading FreeSolv data...')
    import pandas as pd
    # TODO: Don't hard-code these filenames
    csv_filenames = {
        'OpenEye' : 'openeye.csv',
        'AmberTools' : 'ambertools.csv',
        'EspalomaCharge' : 'espaloma-2023-01-05.csv',
    }
    # Assign charges
    toolkit_methods = {
        'OpenEye' : 'am1bcc', 
        'AmberTools' : 'am1bcc',
        'EspalomaCharge' : 'espaloma-am1bcc',
    }
    toolkits = list(toolkit_methods.keys())
    freesolv_results = dict()
    for toolkit_name, csv_filename in csv_filenames.items():
        freesolv_results[toolkit_name] = pd.read_csv(csv_filename)

    # Merge dataframes
    toolkit_names = csv_filenames.keys()
    df = freesolv_results[toolkits[0]]
    df = pd.merge(df, freesolv_results[toolkits[1]][['name', f'calculated hydration free energy {toolkits[1]} (kcal/mol)']], on='name', how='inner', sort=False)
    df = pd.merge(df, freesolv_results[toolkits[2]][['name', f'calculated hydration free energy {toolkits[2]} (kcal/mol)']], on='name', how='inner', sort=False)


    # Filter by error to only retain molecules with hydration free energy error to OpenEye < 2 kcal/mol
    df = df.assign(error=abs(df[f'calculated hydration free energy {toolkits[0]} (kcal/mol)'] - df[f'calculated hydration free energy {toolkits[2]} (kcal/mol)']))
    df = df[df["error"] > 2.0]     

    # Sort
    df = df.reset_index()  # make sure indexes pair with number of rows
    df.sort_values(by='error', ascending=False)

    # Write molecules with annotations
    from openeye import oechem
    with oechem.oemolostream('charges.oeb') as ofs:
        for index, row in df.iterrows():
            smiles = row['SMILES']
            name = row['name']
            from openff.toolkit.topology import Molecule
            molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
            molecule.name = row['name']
            # Geneate a single conformer to use for alignment
            #molecule.generate_conformers(n_conformers=1)
            # Assign charges
            for toolkit, method in toolkit_methods.items():
                delta_g = row[f'calculated hydration free energy {toolkit} (kcal/mol)']
                print(f'Assigning charges for {toolkit}')
                import copy
                charged_molecule = copy.deepcopy(molecule)
                charged_molecule.name = f"{row['name']} {toolkit} {delta_g:8.2f} kcal/mol"
                assign_charges([charged_molecule], toolkit, method, model_url=model_url)
                # Bugfix: For some reason, we need float64 to convert to OpenEye
                #import pint
                #import numpy as np
                #from openff.units import unit
                #unitless_charges = np.zeros(shape=molecule.n_atoms, dtype=np.float64)
                #for atom_index in range(molecule.n_atoms):
                #    unitless_charges[atom_index] = charged_molecule.partial_charges[atom_index] / unit.elementary_charge
                #charged_molecule.partial_charges = unit.Quantity(unitless_charges, unit.elementary_charge)
                #print(type(charged_molecule.partial_charges))
                #print(type(charged_molecule.partial_charges[0].m))
                # Store as OpenEye
                oemol = charged_molecule.to_openeye()
                oechem.OESetSDData(oemol, 'DeltaG', f'{delta_g:8.2f} kcal/mol')
                tagname = 'PARTIAL_CHARGE'
                tag = oechem.OEGetTag(tagname)
                for oeatom in oemol.GetAtoms():
                    oeatom.SetData(tag, oeatom.GetPartialCharge())

                # Write molecule
                oechem.OEWriteMolecule(ofs, oemol)

            ofs.flush()

    # Write 
    #df[['SMILES', 'name', 'dg_error', 'charge_error']].to_csv('errors.csv', index=False)

@click.command()
def render():
    """Render charges for comparison
    """
    tagname = "PARTIAL_CHARGE"
    minvalue = -1.0
    maxvalue = +1.0

    from openeye import oechem, oedepict, oegrapheme
    molecules = []
    with oechem.oemolistream('charges.oeb') as ifs:        
        oemol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, oemol):
            oedepict.OEPrepareDepiction(oemol, False, False)
            molecules.append(oechem.OEGraphMol(oemol))
    print(f'{len(molecules)} molecules read.')
    n_molecules = len(molecules)

    import math
    rows, cols = int(math.ceil(n_molecules/3)), 3
    width, height = 400, 200
    image = oedepict.OEImage(cols*width, rows*height)

    grid = oedepict.OEImageGrid(image, rows, cols)

    opts = oedepict.OE2DMolDisplayOptions(grid.GetCellWidth(), grid.GetCellHeight(),
                                          oedepict.OEScale_Default)
    opts.SetHydrogenStyle(oedepict.OEHydrogenStyle_ExplicitAll)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    opts.SetTitleLocation(oedepict.OETitleLocation_Top)

    propmap = oegrapheme.OE2DPropMap(opts.GetBackgroundColor())
    propmap.SetNegativeColor(oechem.OEDarkRed)
    propmap.SetPositiveColor(oechem.OEDarkBlue)
    propmap.SetLegendLocation(oegrapheme.OELegendLocation_Left)

    propmap.SetMinValue(minvalue)
    propmap.SetMaxValue(maxvalue)

    def get_min_max_atom_property(mol, tagname):
        minvalue = float("inf")
        maxvalue = float("-inf")
        tag = oechem.OEGetTag(tagname)
        for atom in mol.GetAtoms():
            if atom.HasData(tag):
                val = atom.GetData(tag)
                minvalue = min(minvalue, val)
                maxvalue = max(maxvalue, val)

        return minvalue, maxvalue

    index = 0
    for cell, mol in zip(grid.GetCells(), molecules):
        if index % 3 == 0:
            minvalue, maxvalue = get_min_max_atom_property(mol, tagname)
            propmap.SetMinValue(minvalue*1.2)
            propmap.SetMaxValue(maxvalue*1.2)
        
        disp = oedepict.OE2DMolDisplay(mol, opts)
        propmap.Render(disp, tagname)

        oedepict.OERenderMolecule(cell, disp)
        index += 1

    oedepict.OEWriteImage("charges.pdf", image)


cli.add_command(freesolv)
cli.add_command(errors)
cli.add_command(render)

if __name__ == '__main__':
    cli()