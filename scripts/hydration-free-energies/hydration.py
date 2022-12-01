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
@click.option('--forcefield', 
              default='openff-2.0.0',
              help='Small molecule force field to use')
@click.option('--filepath', 
              required=True,
              help='File path to store output')
@click.option('--niterations', 
              default=1000,
              help='Number of iterations to run')
def run(smiles, toolkit, method, forcefield, filepath, niterations):
    """
    Compute the hydration free energy of a specified molecule

    \b
    Parameters
    ----------
    smiles : str
        The SMILES string of the compound whose hydration free energy is to be computed.
        The SMILES string can have explicit or implicit protons.
        The compound must be neutral.
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
    forcefield : str
        Small molecule force field to use
        Valid options depend on openff force fields installed, but include
        ['gaff-1.8', 'gaff-2.1', 'gaff-2.11']
        ['openff-1.2.1', 'openff-1.3.1', 'openff-2.0.0']
            see https://github.com/openforcefield/openff-forcefields for complete list of available openff force field versions
    filepath : str
        The filepath containing the simulation NetCDF and YAML files.   
    niterations : int
        The number of iterations to run
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
    molecule.assign_partial_charges(method, toolkit_registry=toolkit_wrapper)

    # DEBUG
    print(f'Assigned partial charges from toolkit {toolkit} : {method}')
    print(molecule.partial_charges)

    # Generate positions
    molecule.generate_conformers(n_conformers=1)
 
    # Thermodynamic parameters for simulation
    from openmm import unit
    temperature = 298.0 * unit.kelvin
    pressure = 1.0 * unit.atmospheres

    # Set up vacuum and solvent simulations
    import openmm
    from openmm import app
    ffxml_forcefields = ['amber/tip3p_standard.xml'] # OpenMM ffxml solvent force field
    solvent_model = 'tip3p' # solvent model for Modeller.addSolvent()
    padding = 12.0 * unit.angstroms # padding for solvent box construction
    hydrogen_mass = 4.0 * unit.amu
    nonperiodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : hydrogen_mass }
    periodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : hydrogen_mass, 'nonbondedMethod' : app.PME }
    from openmmforcefields.generators import SystemGenerator
    barostat = openmm.MonteCarloBarostat(pressure, temperature) # reference barostat for adding barostat to generated systems
    system_generator = SystemGenerator(forcefields=ffxml_forcefields, small_molecule_forcefield=forcefield, 
        nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs, periodic_forcefield_kwargs=periodic_forcefield_kwargs)

    # Force the SystemGenerator to generate and memorize parameters for the small molecule
    system_generator.create_system(molecule.to_topology().to_openmm(), molecules=[molecule])

    # Generate OpenMM Topology and positions for each phase
    openff_topology = dict() # openff_topology[phase] is the OpenFF Topology for phase, where phase is one of ['vacuum', 'solvent']

    #
    # Vacuum phase setup
    #

    openff_topology['vacuum'] = molecule.to_topology() # OpenFF Topology

    #
    # Solvent phase setup
    #

    # Use openmm.app.Modeller to add solvent
    from openmm.app import Modeller
    modeller = Modeller(molecule.to_topology().to_openmm(), molecule.conformers[0].to_openmm())
    modeller.addSolvent(system_generator.forcefield, model=solvent_model, padding=padding)
    # Transform OpenMM Topology and positions back into OpenFF Topology
    water = Molecule.from_smiles('HOH')
    from openff.toolkit.topology import Topology
    openff_topology['solvent'] = Topology.from_openmm(modeller.topology, unique_molecules=[molecule, water])
    import openff.units.openmm
    openff_topology['solvent'].set_positions( openff.units.openmm.from_openmm(modeller.positions) )

    #
    # Run simulations
    #

    phases = ['vacuum', 'solvent']
    for phase in phases:
        # Modify the SystemGenerator to work around bug in openmmforcefields <=0.11.2 : https://github.com/openmm/openmmforcefields/issues/252
        match phase:
            case 'solvent':
                system_generator.barostat = barostat
            case 'vacuum':
                system_generator.barostat = None

        # Create the OpenMM System
        system = system_generator.create_system(openff_topology[phase].to_openmm(), molecules=[molecule])

        # Define the thermodynamic state
        from openmmtools.states import ThermodynamicState
        match phase:
            case 'solvent':
                thermodynamic_state = ThermodynamicState(system=system, temperature=temperature, pressure=pressure)
            case 'vacuum':
                thermodynamic_state = ThermodynamicState(system=system, temperature=temperature)

        # Run the alchemical free energy calculation for this phase
        run_phase(molecule, system, openff_topology[phase], thermodynamic_state, phase, filepath, niterations)

def run_phase(molecule, system, topology, thermodynamic_state, phase, filepath, niterations):
    """Run an alchemical free energy calculation for a single phase.

    Parameters
    ----------
    molecule : openff.toolkit.topology.Molecule
        The molecule to be alchemically annihilated/decoupled
    system : openmm.System
        The System to alchemically modify
    topology : openff.toolkit.topology.Topology
        OpenFF Topology associated with the System
    thermodynamic_state : openmmtools.states.ThermodynamicState
        The thermodynamic state to run the simulation at
    phase : str
        The phase to simulate: ['vacuum', 'solvent']
    filepath : str
        The filepath containing the simulation NetCDF and YAML files.  
    niterations : int
        The number of iterations to run
    """
    ALLOWED_PHASES = ['vacuum', 'solvent']
    if phase not in ALLOWED_PHASES:
        raise ValueError(f"phase must be one of {ALLOWED_PHASES}; specified '{phase}'")

    # Write out PDB file
    #with open(phase + '.pdb', 'wt') as outfile:
    #    from openmm.app import PDBFile
    #    PDBFile.writeFile(topology.to_openmm(), topology.get_positions().to_openmm(), outfile)

    # Simulation parameters
    import openmm
    from openmm import unit
    collision_rate = 1.0 / unit.picoseconds
    timestep = 4.0 * unit.femtoseconds
    n_steps = 250 # number of steps per iteration
    reassign_velocities = False # whether to reassign velocities every iteration

    # Define the region of the System to be alchemically modified.
    from openmmtools.alchemy import AlchemicalRegion, AbsoluteAlchemicalFactory, AlchemicalState
    alchemical_atoms = list(range(molecule.n_atoms))
    #alchemical_region = AlchemicalRegion(alchemical_atoms=alchemical_atoms, softcore_beta=1.0)
    #factory = AbsoluteAlchemicalFactory(alchemical_pme_treatment='direct-space')
    alchemical_region = AlchemicalRegion(alchemical_atoms=alchemical_atoms)
    factory = AbsoluteAlchemicalFactory()
    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    # Simulation parameters
    import os
    storage_path = os.path.join(filepath, phase + '.nc')
    online_analysis_interval = 25
    match phase:
        case 'vacuum':
            n_lambda = 6 # number of alchemical states
            platform = openmm.Platform.getPlatformByName('CPU')
        case 'solvent':
            n_lambda = 30 # number of alchemical states
            platform = openmm.Platform.getPlatformByName('CUDA')
            platform.setPropertyDefaultValue('Precision', 'mixed')
            platform.setPropertyDefaultValue('DeterministicForces', 'true')

    # Define sampler
    from openmmtools.mcmc import LangevinDynamicsMove
    mcmc_move = LangevinDynamicsMove(timestep=timestep, collision_rate=collision_rate, reassign_velocities=reassign_velocities, n_steps=n_steps, n_restart_attempts=6)

    # Initialize compound thermodynamic states at different temperatures and alchemical states.
    import numpy as np

    #protocol = {'temperature': thermodynamic_state.temperature * np.ones([n_lambda]),
    #            'lambda_electrostatics': np.linspace(1, 0, n_lambda),
    #            'lambda_sterics': np.linspace(1, 0, n_lambda)}
    # DEBUG
    protocol = dict()
    protocol['lambda_electrostatics'] = np.array([1.00, 0.75, 0.50, 0.37, 0.25, 0.12, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
    protocol['lambda_sterics']        = np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00])
    n_lambda = len(protocol['lambda_electrostatics'])

    protocol['temperature'] = thermodynamic_state.temperature * np.ones([n_lambda])
    if system.usesPeriodicBoundaryConditions():
        protocol['pressure'] = thermodynamic_state.pressure * np.ones([n_lambda])

    alchemical_state = AlchemicalState.from_system(alchemical_system)
    from openmmtools.states import create_thermodynamic_state_protocol
    compound_states = create_thermodynamic_state_protocol(alchemical_system, protocol=protocol, composable_states=[alchemical_state])

    # Initialize sampler states
    from openmmtools.states import SamplerState
    sampler_state = SamplerState(positions=topology.get_positions().to_openmm())
    if system.usesPeriodicBoundaryConditions():
        sampler_state.box_vectors = system.getDefaultPeriodicBoxVectors()        
    sampler_states = [ sampler_state for _ in compound_states ]

    # Run the combined Hamiltonian replica exchange + parallel tempering simulation.
    from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
    sampler = ReplicaExchangeSampler(mcmc_moves=mcmc_move, number_of_iterations=niterations)
    reporter = MultiStateReporter(storage_path, checkpoint_interval=niterations)
    sampler.create(thermodynamic_states=compound_states, sampler_states=sampler_states, storage=reporter)
    sampler.online_analysis_interval = online_analysis_interval

    # Setup context cache for multistate samplers
    from openmmtools.cache import ContextCache
    context_cache = ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.energy_context_cache = context_cache
    sampler.sampler_context_cache = context_cache

    # Minimize
    sampler.minimize()

    # Run the sampler
    from rich.progress import track
    for iteration in track(range(niterations), description=f"Running {phase} phase..."):
        sampler.run(1)

    # Clean up to free up resources
    del sampler
    del context_cache

@click.command()
@click.option('--filepath', 
              required=True,
              help='File path to store output')
def analyze(filepath):
    """
    Analyze the hydration free energy of a specified molecule

    \b
    Parameters
    ----------
    filepath : str
        The filepath containing the simulation NetCDF and YAML files.

    """
    import yaml
    import os

    # Thermodynamic parameters for simulation
    # TODO: Find a way to harmonize temperature between methods
    from openmm import unit
    temperature = 298.0 * unit.kelvin

    # TODO: Improve how statistics are handled

    free_energy = dict()    
    phases = ['vacuum', 'solvent']
    for phase in phases:
        filename = os.path.join(filepath, phase +  '_real_time_analysis.yaml')
        with open(filename, 'rt') as infile:
            estimates = yaml.safe_load(infile)
        df = estimates[-1]['mbar_analysis']['free_energy_in_kT']
        ddf = estimates[-1]['mbar_analysis']['standard_error_in_kT']
        free_energy[phase] = (df, ddf)
    
    # Compute hydration free energy in kT
    # TODO: This part should probably use uncertainties: https://pypi.org/project/uncertainties/
    import numpy as np
    df = free_energy['vacuum'][0] - free_energy['solvent'][0]
    ddf = np.sqrt(free_energy['vacuum'][1]**2 + free_energy['solvent'][1]**2)

    # Convert to kcal/mol
    from openmm import unit
    from openmmtools.constants import kB
    kT = kB * temperature
    df *= (kT / unit.kilocalories_per_mole)
    ddf *= (kT / unit.kilocalories_per_mole)

    # TODO: Format this with appropriate sig figs
    print(f'hydration free energy: {df} +- {ddf} kcal/mol')

    # TODO: Write out to a file?

    return df, ddf

cli.add_command(run)
cli.add_command(analyze)

if __name__ == '__main__':
    cli()