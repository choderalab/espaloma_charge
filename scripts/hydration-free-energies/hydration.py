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

    # DEBUG
    print(f'Assigned partial charges from toolkit {toolkit} : {method}')
    print(molecule.partial_charges)

    #
    # TODO: Refactor the below
    #

    # Generate positions
    molecule.generate_conformers(n_conformers=1)
    positions = molecule.conformers[0].to_openmm()

    # Create vacuum OpenMM topology
    openmm_topology = molecule.to_topology().to_openmm()
    # Create vacuum system
    from openmm import unit
    from openmm import app
    forcefields = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml']
    small_molecule_forcefield = 'gaff-2.11'
    nonperiodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
    periodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu, 'nonbondedMethod' : app.PME }
    cache = 'db.json'
    from openmmforcefields.generators import SystemGenerator
    system_generator = SystemGenerator(forcefields=forcefields, small_molecule_forcefield=small_molecule_forcefield, 
        nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs, periodic_forcefield_kwargs=periodic_forcefield_kwargs, cache=cache)
    # Create an OpenMM System from an OpenMM Topology object
    system = system_generator.create_system(openmm_topology, molecules=[molecule])

    # Simulation parameters
    from openmm import unit
    temperature = 298.0 * unit.kelvin
    collision_rate = 5.0 / unit.picoseconds
    timestep = 4.0 * unit.femtoseconds
    n_steps = 250 # number of steps per iteration
    reassign_velocities = False # whether to reassign velocities every iteration

    #
    # Vacuum phase
    #

    # Define the region of the System to be alchemically modified.
    from openmmtools.alchemy import AlchemicalRegion, AbsoluteAlchemicalFactory, AlchemicalState
    alchemical_atoms = list(range(molecule.n_atoms))
    alchemical_region = AlchemicalRegion(alchemical_atoms=alchemical_atoms)
    factory = AbsoluteAlchemicalFactory()
    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    # Simulation parameters
    n_lambda = 6 # number of alchemical states
    n_iterations = 1000
    storage_path = 'vacuum.nc'

    # Define sampler
    from openmmtools.mcmc import LangevinDynamicsMove
    mcmc_move = LangevinDynamicsMove(timestep=timestep, collision_rate=collision_rate, reassign_velocities=reassign_velocities, n_steps=n_steps, n_restart_attempts=6)

    # Initialize compound thermodynamic states at different temperatures and alchemical states.
    import numpy as np
    protocol = {'temperature': temperature * np.ones([n_lambda]),
                'lambda_electrostatics': np.linspace(1, 0, n_lambda),
                'lambda_sterics': np.linspace(1, 0, n_lambda)}
    alchemical_state = AlchemicalState.from_system(alchemical_system)
    from openmmtools.states import create_thermodynamic_state_protocol
    compound_states = create_thermodynamic_state_protocol(alchemical_system, protocol=protocol, composable_states=[alchemical_state])

    # Initialize sampler states
    from openmmtools.states import SamplerState
    sampler_states = [ SamplerState(positions=positions) for _ in compound_states ]

    # Run the combined Hamiltonian replica exchange + parallel tempering simulation.
    from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
    sampler = ReplicaExchangeSampler(mcmc_moves=mcmc_move, number_of_iterations=n_iterations)
    reporter = MultiStateReporter(storage_path, checkpoint_interval=n_iterations)
    sampler.create(thermodynamic_states=compound_states, sampler_states=sampler_states, storage=reporter)

    # Run the sampler
    from rich.progress import track
    for iteration in track(range(n_iterations), description="Running vacuum phase..."):
        sampler.run(1)
    

cli.add_command(run)

if __name__ == '__main__':
    cli()