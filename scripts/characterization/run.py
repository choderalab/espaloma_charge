import time
import numpy as np
import torch
import pandas as pd
import scipy
from openff.toolkit.topology import Molecule
from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
from openff.toolkit.utils.toolkit_registry import AmberToolsToolkitWrapper, OpenEyeToolkitWrapper

from typing import Union

import numpy as np

from openff.units import unit
from openff.toolkit.topology import Molecule
from openff.recharge.grids import GridGenerator, MSKGridSettings


def calculate_esp(
    grid_coordinates: unit.Quantity,  # N x 3
    atom_coordinates: unit.Quantity,  # M x 3
    charges: unit.Quantity,  # M
    with_units: bool = False,
) -> Union[float, unit.Quantity]:
    """
    Calculate the electrostatic potential at a set of grid points

    Parameters
    ----------
    grid_coordinates
        The coordinates of the points at which to calculate the ESP
    atom_coordinates
        The coordinates of the atoms in the molecule
    charges
        The charges of the atoms in the molecule.
        Should be in same order as ``atom_coordinates``
    with_units
        Whether to return the ESP as a float in atomic units
        or as a unit.Quantity

    Returns
    -------
    esp: float or unit.Quantity
        The ESP at the grid points
    """
    # make sure charges are in correct unit
    # charges = unify(charges)
    charges = charges * unit.elementary_charge

    AU_ESP = unit.atomic_unit_of_energy / unit.elementary_charge
    ke = 1 / (4 * np.pi * unit.epsilon_0)

    displacement = (
        grid_coordinates[:, None, :]
        - atom_coordinates[None, :, :]  # N x M x 3
    )
    distance = (displacement ** 2).sum(axis=-1) ** 0.5  # N x M
    inv_distance = 1 / distance

    esp = ke * (inv_distance @ charges)  # N

    esp_q = esp.m_as(AU_ESP)
    if not with_units:
        return esp_q
    return esp



def compare_molecule_esp(
    molecule: Molecule,
    charges1: unit.Quantity,
    charges2: unit.Quantity,
) -> float:
    """Calculate the RMSD of the ESPs a molecule with different charges

    Parameters
    ----------
    molecule: openff.toolkit.topology.Molecule
        The molecule to calculate the ESPs for
    charges1: unit.Quantity
        The charges to use for the first ESP
    charges2: unit.Quantity
        The charges to use for the second ESP

    Returns
    -------
    rmsd: float
    """

    assert len(molecule.conformers) > 0, "Molecule has no conformers"

    settings = MSKGridSettings()
    rmses = np.zeros((len(molecule.conformers),), dtype=float)

    for i, conformer in enumerate(molecule.conformers):
        grid = GridGenerator.generate(molecule, conformer, settings)
        esp1 = calculate_esp(grid, conformer, charges1)
        esp2 = calculate_esp(grid, conformer, charges2)
        delta = esp1 - esp2
        rmsd = (delta ** 2).mean() ** 0.5
        rmses[i] = rmsd

    return rmses.mean()

def generate_elf10_conformers(
    molecule: Molecule,
    n_conformer_pool: int = 500,
    n_conformers: int = 10,
    rms_cutoff: Union[float, unit.Quantity] = 0.05,
):
    """Generate ELF conformers for a molecule inplace

    Parameters
    ----------
    molecule: openff.toolkit.topology.Molecule
        The molecule to generate conformers for
    n_conformer_pool: int
        The number of conformers to generate
        to select ELF conformers from
    n_conformers: int
        The number of ELF conformers to select.
        The molecule will have maximum this number
        of conformers
    rms_cutoff: float or unit.Quantity
        The RMSD cutoff to use when selecting conformers.
        If a float is provided, it will be interpreted as
        an Angstrom RMSD cutoff

    Returns
    -------
    molecule: openff.toolkit.topology.Molecule
        The same molecule, with ELF conformers
    """
    if isinstance(rms_cutoff, float):
        rms_cutoff = rms_cutoff * unit.angstrom

    molecule.generate_conformers(
        n_conformers=n_conformer_pool,
        rms_cutoff=rms_cutoff,
        make_carboxylic_acids_cis=True,
    )
    molecule.apply_elf_conformer_selection(limit=n_conformers)
    return molecule

def compute_zap(molecule):
    from openeye import oezap, oegrid, oechem

    mol = molecule.to_openeye()
    oechem.OEAssignBondiVdWRadii(mol)
    ## This is only necessary if the molecule doesn't already have 3D
    ##  coordinates
    oechem.OEGenerate2DCoordinates(mol)
    mol.SetDimension(3)

    # ## Already handled automatically when converting to OE
    # charges = molecule.partial_charges.magnitude
    # for atom in mol.GetAtoms():
    #     atom.SetPartialCharge(charges[atom.GetIdx()])

    zap = oezap.OEZap()
    zap.SetInnerDielectric(1.0)
    zap.SetGridSpacing(0.5)
    zap.SetMolecule(mol)

    grid = oegrid.OEScalarGrid()
    zap.CalcPotentialGrid(grid)
    atom_potentials = oechem.OEFloatArray(mol.GetMaxAtomIdx())
    zap.CalcAtomPotentials(atom_potentials)
    return np.array(atom_potentials)

def run(args):
    smiles = open(args.in_path, "r").readlines()[args.idx].strip()
    molecule = Molecule.from_smiles(smiles)

    toolkit_registry = EspalomaChargeToolkitWrapper()
    molecule.assign_partial_charges('espaloma-am1bcc', toolkit_registry=toolkit_registry)
    time0 = time.time()
    molecule.assign_partial_charges('espaloma-am1bcc', toolkit_registry=toolkit_registry)
    time1 = time.time()
    q_esp = molecule.partial_charges.magnitude.tolist()
    t_esp = time1 - time0
    zap_esp = compute_zap(molecule)

    toolkit_registry = AmberToolsToolkitWrapper()
    molecule.assign_partial_charges('am1bcc', toolkit_registry=toolkit_registry)
    time0 = time.time()
    molecule.assign_partial_charges('am1bcc', toolkit_registry=toolkit_registry)
    time1 = time.time()
    q_amber = molecule.partial_charges.magnitude.tolist()
    t_amber = time1 - time0
    zap_amber = compute_zap(molecule)

    toolkit_registry = OpenEyeToolkitWrapper()
    molecule.assign_partial_charges('am1bcc', toolkit_registry=toolkit_registry)
    time0 = time.time()
    molecule.assign_partial_charges('am1bcc', toolkit_registry=toolkit_registry)
    time1 = time.time()
    q_oe = molecule.partial_charges.magnitude.tolist()
    t_oe = time1 - time0
    zap_oe = compute_zap(molecule)

    generate_elf10_conformers(molecule)

    rmse_esp_amber = compare_molecule_esp(molecule, q_oe, q_amber)
    rmse_esp_espaloma = compare_molecule_esp(molecule, q_oe, q_esp)

    df = pd.DataFrame.from_dict(
        {
            smiles: 
                [
                    q_amber, q_esp, q_oe, 
                    t_amber, t_esp, t_oe,
                    zap_amber, zap_esp, zap_oe,
                    rmse_esp_amber, rmse_esp_espaloma,
                ]
        }, 
        orient="index")

    print(df)
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--idx", type=int)
    parser.add_argument("--in_path", type=str)
    parser.add_argument("--out", type=str)
    args = parser.parse_args()
    df = run(args)
    df.to_csv(args.out, mode="a", header=None)







    
    

