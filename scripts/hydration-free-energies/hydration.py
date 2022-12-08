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
              default=5000,
              help='Number of iterations to run')
@click.option('--write-pdb', 
              is_flag=True,
              help='Write PDB files of initial models')   
@click.option('--model-url', 
              default=None,
              help='If specified, provide model URL or local filepath to espaloma_charge')   
def run(smiles, toolkit, method, forcefield, filepath, niterations, write_pdb):
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
            'EspalomaCharge' : ['espaloma-am1bcc']
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
    write_pdb : bool
        If True, write PDB file of intiial models 
    model_url : str, optional, default=None
        Provide optional model_url or filepath to espaloma_charge.app.charge()       
    """
    run_all(smiles, toolkit, method, forcefield, filepath, niterations, write_pdb, model_url)

def run_all(smiles, toolkit, method, forcefield, filepath, niterations, write_pdb, model_url):
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
            'EspalomaCharge' : ['espaloma-am1bcc']
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
    write_pdb : bool
        If True, write PDB file of intiial models        
    model_url : str
        Provide optional model_url or filepath to espaloma_charge.app.charge()       
    """
    # Terminate immediately if both phases are complete
    PHASES = ['vacuum', 'solvent']
    import os
    # TODO: Check that simulations have finished by checking YAML files for expected length?
    if all([ os.path.exists(os.path.join(filepath, phase + '.nc')) for phase in PHASES ]):
        print('All phases complete already. Exiting.')
        return

    # Create an OpenFF Molecule object
    from openff.toolkit.topology import Molecule
    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    # Check to make sure molecule is neutral
    if molecule.total_charge.magnitude != 0.0:
        raise ValueError(f'Molecule {smiles} has a net charge of {molecule.total_charge}. Only neutral molecules are supported.')

    # Assign partial charges
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
    padding = 14.0 * unit.angstroms # padding for solvent box construction
    hydrogen_mass = 3.8 * unit.amu
    nonbonded_cutoff = 9.0 * unit.angstroms
    nonperiodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : hydrogen_mass }
    periodic_forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : hydrogen_mass, 'nonbondedMethod' : app.PME, 'nonbondedCutoff' : nonbonded_cutoff }
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

    PHASES = ['vacuum', 'solvent']
    for phase in PHASES:
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
        run_phase(molecule, system, openff_topology[phase], thermodynamic_state, phase, filepath, niterations, write_pdb)

def run_phase(molecule, system, topology, thermodynamic_state, phase, filepath, niterations, write_pdb):
    """Run an alchemical free energy calculation for a single phase.
    The phase will be skipped if the .nc file already exists.

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
    write_pdb : bool
        If True, write PDB file of intiial models
    """
    import logging
    logger = logging.getLogger(__name__)

    ALLOWED_PHASES = ['vacuum', 'solvent']
    if phase not in ALLOWED_PHASES:
        raise ValueError(f"phase must be one of {ALLOWED_PHASES}; specified '{phase}'")

    # Skip if phase already exists
    import os
    storage_path = os.path.join(filepath, phase + '.nc')
    # TODO: Try to resume instead of skip
    if os.path.exists(storage_path):
        print(f'{storage_path} already exists; skipping this phase.')
        return

    # Activate logging
    import logging
    import os
    if not os.path.exists(filepath):
        os.makedirs(filepath, exist_ok=True)
    log_filename = os.path.join(filepath, phase + '.log')
    logging.basicConfig(filename=log_filename, encoding='utf-8', level=logging.INFO)

    # Get initial positions
    positions = topology.get_positions().to_openmm()

    # Write out initial PDB file
    if write_pdb:
        import os
        filename = os.path.join(filepath, phase + '-initial.pdb')
        with open(filename, 'wt') as outfile:
            from openmm.app import PDBFile
            PDBFile.writeFile(topology.to_openmm(), positions, outfile)

    # Simulation parameters
    import os
    import openmm
    online_analysis_interval = 25
    checkpoint_interval = 100
    match phase:
        case 'vacuum':
            n_lambda = 8 # number of alchemical states
            platform = openmm.Platform.getPlatformByName('CPU')

            MAX_VACUUM_ITERATIONS = 200
            if niterations > MAX_VACUUM_ITERATIONS:
                print(f'capping number of vacuum iterations at {MAX_VACUUM_ITERATIONS}')
                niterations = MAX_VACUUM_ITERATIONS
        case 'solvent':
            n_lambda = 24 # number of alchemical states
            platform = openmm.Platform.getPlatformByName('CUDA')
            platform.setPropertyDefaultValue('Precision', 'mixed')
            platform.setPropertyDefaultValue('DeterministicForces', 'true')

    # Minimize
    logger.debug('Minimizing...')
    import openmm
    integrator = openmm.VerletIntegrator(0.001)
    context = openmm.Context(system, integrator, platform)
    context.setPositions(positions)
    logger.debug(f'  Initial energy: {context.getState(getEnergy=True).getPotentialEnergy()}')
    openmm.LocalEnergyMinimizer.minimize(context)
    logger.debug(f'  Final energy:   {context.getState(getEnergy=True).getPotentialEnergy()}')
    positions = context.getState(getPositions=True).getPositions(asNumpy=True)
    del context, integrator

    # Write out initial PDB file
    if write_pdb:
        import os
        filename = os.path.join(filepath, phase + '-minimized.pdb')
        with open(filename, 'wt') as outfile:
            from openmm.app import PDBFile
            PDBFile.writeFile(topology.to_openmm(), positions, outfile)

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

    # Standard alchemical factory
    alchemical_region = AlchemicalRegion(alchemical_atoms=alchemical_atoms)
    factory = AbsoluteAlchemicalFactory()
    # Fused softcore alchemical factory (experimental)
    #alchemical_region = AlchemicalRegion(alchemical_atoms=alchemical_atoms, softcore_beta=1.0)
    #factory = AbsoluteAlchemicalFactory(alchemical_pme_treatment='direct-space')

    alchemical_system = factory.create_alchemical_system(system, alchemical_region)

    # Define sampler
    from openmmtools.mcmc import LangevinDynamicsMove
    mcmc_move = LangevinDynamicsMove(timestep=timestep, collision_rate=collision_rate, reassign_velocities=reassign_velocities, n_steps=n_steps, n_restart_attempts=6)

    # Initialize compound thermodynamic states at different temperatures and alchemical states.
    import numpy as np
    protocol = dict()
    # TODO: Use thermodynamic trailblazing to automatically determine number and spacing of alchemical states

    # Two phase
    n1 = int(n_lambda/2) # number of steps for turning off electrostatics
    n2 = n_lambda - int(n_lambda/2) + 1 # number of steps for turning off Lennard-Jones
    protocol['lambda_electrostatics'] = np.concatenate( [np.linspace(1, 0, n1), np.linspace(0, 0, n2)[1:]] )
    protocol['lambda_sterics']        = np.concatenate( [np.linspace(1, 1, n1), np.linspace(1, 0, n2)[1:]] )
    assert len(protocol['lambda_electrostatics']) == n_lambda

    # One phase (requires softcore)
    #protocol['lambda_electrostatics'] = np.linspace(1, 0, n_lambda)
    #protocol['lambda_sterics']        = np.linspace(1, 0, n_lambda)
    #assert len(protocol['lambda_electrostatics']) == n_lambda

    protocol['temperature'] = thermodynamic_state.temperature * np.ones([n_lambda])
    if system.usesPeriodicBoundaryConditions():
        protocol['pressure'] = thermodynamic_state.pressure * np.ones([n_lambda])

    alchemical_state = AlchemicalState.from_system(alchemical_system)
    from openmmtools.states import create_thermodynamic_state_protocol
    compound_states = create_thermodynamic_state_protocol(alchemical_system, protocol=protocol, composable_states=[alchemical_state])

    # Initialize sampler states
    from openmmtools.states import SamplerState
    sampler_state = SamplerState(positions=positions)
    if system.usesPeriodicBoundaryConditions():
        sampler_state.box_vectors = system.getDefaultPeriodicBoxVectors()        
    sampler_states = [ sampler_state for _ in compound_states ]

    # Run the combined Hamiltonian replica exchange + parallel tempering simulation.
    from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
    sampler = ReplicaExchangeSampler(mcmc_moves=mcmc_move, number_of_iterations=niterations)
    reporter = MultiStateReporter(storage_path, checkpoint_interval=checkpoint_interval)
    sampler.create(thermodynamic_states=compound_states, sampler_states=sampler_states, storage=reporter)
    sampler.online_analysis_interval = online_analysis_interval

    # Setup context cache for multistate samplers
    from openmmtools.cache import ContextCache
    sampler.energy_context_cache = ContextCache(capacity=None, time_to_live=None, platform=platform)
    sampler.sampler_context_cache = ContextCache(capacity=None, time_to_live=None, platform=platform)
    
    # Minimize
    #sampler.minimize()

    # Run the sampler
    from rich.progress import track
    for iteration in track(range(niterations), description=f"Running {phase} phase..."):
        sampler.run(1)

    # Clean up to free up resources
    del sampler.energy_context_cache
    del sampler.sampler_context_cache
    del sampler

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
    df, ddf = analyze_single_calculation(filepath)
    # TODO: Format this with appropriate sig figs
    print(f'hydration free energy: {df} +- {ddf} kcal/mol')


def analyze_single_calculation(filepath):
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

    # TODO: Write out to a file?

    return df, ddf

@click.command()
@click.option('--index', 
              required=True,
              type=click.INT,
              help='Index (starting with 1) into freesolv compound database')
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
              default=5000,
              help='Number of iterations to run')
@click.option('--write-pdb', 
              is_flag=True,
              help='Write PDB files of initial models')   
@click.option('--model-url', 
              default=None,
              help='If specified, provide model URL or local filepath to espaloma_charge')
def freesolv(index, toolkit, method, forcefield, filepath, niterations, write_pdb, model_url):
    """Run specified molecule index from FreeSolv database
    """
    # Load FreeSolv database
    # TODO: Retrieve automatically if needed
    import json
    with open('freesolv.json', 'rt') as infile:
        freesolv = json.load(infile)

    # Extract info
    name = list(freesolv.keys())[index - 1]
    smiles = freesolv[name]['smiles']

    # Run free energy calculation
    import os
    run_all(smiles, toolkit, method, forcefield, os.path.join(filepath, name), niterations, write_pdb, model_url)

@click.command()
@click.option('--filepath', 
              required=True,
              help='Base file path containing FreeSolv calculations')
@click.option('--label', 
              required=True,
              help='Label for calculations')
@click.option('--outfile', 
              required=True,
              help='Output filename for CSV')
def analyze_freesolv(filepath, label, outfile):
    """Analyze 
    """
    # Load FreeSolv database
    # TODO: Retrieve automatically if needed
    import json
    with open('freesolv.json', 'rt') as infile:
        freesolv = json.load(infile)
    
    # DEBUG
    #freesolv = { key : freesolv[key] for key in list(freesolv.keys())[0:20] }

    records = list()
    for name, freesolv_entry in freesolv.items():
        record = dict()

        record['name'] = name
        record['SMILES'] = freesolv_entry['smiles']
        record['IUPAC name'] = freesolv_entry['iupac']
        
        record['experimental hydration free energy (kcal/mol)'] = freesolv_entry['expt'] # kcal/mol
        record['experimental hydration free energy uncertainty (kcal/mol)'] = freesolv_entry['d_expt'] # measurement uncertainty in kcal/mol
        record['experimental hydration free energy DOI'] = freesolv_entry['expt_reference']

        import os
        import numpy as np
        try:
            smiles = freesolv_entry['smiles']
            DeltaG_exp = freesolv_entry['expt'] # kcal/mol
            dDeltaG_exp = freesolv_entry['d_expt'] # measurement uncertainty in kcal/mol

            DeltaG_calc, dDeltaG_calc = analyze_single_calculation(os.path.join(filepath, name)) # in kcal/mol
            DeltaG_error = DeltaG_calc - DeltaG_exp
            dDeltaG_error = np.sqrt(dDeltaG_calc**2 + dDeltaG_exp**2)
    
            print(f'{name:20s} {smiles:85s} : exp {DeltaG_exp:8.2f} +- {dDeltaG_exp:5.2f} kcal/mol | calc {DeltaG_calc:8.2f} +- {dDeltaG_calc:5.2f} kcal/mol | calc-exp  {DeltaG_error:8.2f} +- {dDeltaG_error:5.2f} kcal/mol')

            record[f'calculated hydration free energy {label} (kcal/mol)'] = DeltaG_calc
            record[f'calculated hydration free energy uncertainty {label} (kcal/mol)'] = dDeltaG_calc

            record[f'hydration free energy {label} error (kcal/mol)'] = DeltaG_error
            record[f'hydration free energy uncertainty {label} error (kcal/mol)'] = dDeltaG_error

            records.append(record)
        except FileNotFoundError as e:
            pass

    # Write details
    import pandas as pd
    df = pd.DataFrame.from_records(records, columns=records[0].keys())
    print(df)
    df.to_csv(outfile, index=False)
    n_compounds = len(df)

    # Compute summary statistics
    import numpy as np
    import scipy
    def extract_bootstrap_sample(indices):
        "Retrieve calculated, experimental, and error data (in kcal/mol) given compound indices into dataframe"
        calculated = np.array(df.iloc[indices][f'calculated hydration free energy {label} (kcal/mol)'])
        experimental = np.array(df.iloc[indices][f'experimental hydration free energy (kcal/mol)'])
        errors = calculated - experimental
        return calculated, experimental, errors
    # Define statistics functions to use, with docstrings used as titles
    def mean_signed_error(indices):
        "Mean signed error"
        calculated, experimental, errors = extract_bootstrap_sample(indices)
        return np.mean(errors)
    def root_mean_squared_error(indices):
        "Root mean squared error"
        calculated, experimental, errors = extract_bootstrap_sample(indices)
        return np.sqrt((errors**2).mean())
    def mean_unsigned_error(indices):
        "Mean unsigned error"
        calculated, experimental, errors = extract_bootstrap_sample(indices)
        return np.abs(errors).mean()
    def kendall_tau(indices):
        "Kendall tau"
        calculated, experimental, errors = extract_bootstrap_sample(indices)
        return scipy.stats.kendalltau(calculated, experimental)[0]
    def pearson_r(indices):
        "Pearson R"
        calculated, experimental, errors = extract_bootstrap_sample(indices)
        return scipy.stats.pearsonr(calculated, experimental)[0]        
    statistics = [mean_signed_error, root_mean_squared_error, mean_unsigned_error, kendall_tau, pearson_r]
    # Compute statistics
    computed_statistics = dict()
    all_indices = np.array(range(n_compounds))
    for statistic in statistics:
        name = statistic.__doc__
        computed_statistics[name] = dict()
        computed_statistics[name]['mle'] = statistic(all_indices)
        bootstrap_result = scipy.stats.bootstrap((all_indices,), statistic, confidence_level=0.95, vectorized=False, method='percentile')
        computed_statistics[name]['95% CI low'] = bootstrap_result.confidence_interval.low
        computed_statistics[name]['95% CI high'] = bootstrap_result.confidence_interval.high

    for name, value in computed_statistics.items():
        print(f"{name:25} {value['mle']:8.4f} [{value['95% CI low']:8.4f},{value['95% CI high']:8.4f}]")

    # Generate plot
    import matplotlib as mpl
    mpl.use('Agg')
    import seaborn
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=[6,6])
    plt.errorbar(
        df[f'experimental hydration free energy (kcal/mol)'], 
        df[f'calculated hydration free energy {label} (kcal/mol)'], 
        fmt='k.',
        xerr=df[f'experimental hydration free energy uncertainty (kcal/mol)'], 
        yerr=df[f'calculated hydration free energy uncertainty {label} (kcal/mol)'],
        )
    xmin = min(df[f'experimental hydration free energy (kcal/mol)'].min(), df[f'calculated hydration free energy {label} (kcal/mol)'].min()) - 0.5
    xmax = max(df[f'experimental hydration free energy (kcal/mol)'].max(), df[f'calculated hydration free energy {label} (kcal/mol)'].max()) + 0.5
    plt.plot([xmin, xmax], [xmin, xmax], 'k-', linewidth=1)
    plt.axis([xmin, xmax, xmin, xmax])

    title = f"{label}"
    plt.title(title)

    statistics_text = f"N = {n_compounds} compounds\n"
    for name, value in computed_statistics.items():
        statistics_text += f"{name} {value['mle']:.2f} [{value['95% CI low']:.2f},{value['95% CI high']:.2f}] kcal/mol\n"
    plt.legend([statistics_text], fontsize=7)

    plt.xlabel('experimental hydration free energy (kcal/mol)')
    plt.ylabel('calculated hydration free energy (kcal/mol)')
    plt.tight_layout()
    figure_filename = label + '.pdf'
    plt.savefig(figure_filename)
    print(f'Figure written to {figure_filename}')

cli.add_command(run)
cli.add_command(analyze)
cli.add_command(freesolv)
cli.add_command(analyze_freesolv)

if __name__ == '__main__':
    cli()