# Computing explicit solvent hydration free energies with espaloma charges

This directory contains scripts for computing explicit solvent hydration free energies with esploma charges, 
and comparing these with AmberTools antechamber (sqm) charges and OpenEye Toolkit ELF10 charges
on the [FreeSolv dataset](https://github.com/MobleyLab/FreeSolv) of neutral molecules.

## Prerequisites
Install the prerequisites with `conda` or `mamba`:
```
$ mamba env create -f environment.yml -n hydration
$ mamba activate hydration
```

## Computing hydration free energies

### Running a simulation for an arbitrary molecule

To set up and run a simulation, use
```bash
$ python hydration.py run --smiles <smiles> --toolkit <toolkit> --method <method> --niterations <niterations> --forcefield <forcefield> --filepath <filepath>
```
`python hydration.py run --help` will provide full help on arguments, but briefly:
* The `--smiles <smiles>` string specifies the SMILES string of the molecule you want to compute.
* The `--toolkit <toolkit>` and `--method <method>` arguments can be used to specify the toolkit (`EspalomaCharge`, `AmberTools`, `OpenEye`, `RDKit`) and charge method (depends on toolkit).
* The `--niterations <niterations>` specifies how many iterations (1 ps/iteration) to run (default: 1000).
* The `--forcefield <forcefield>` argument can specify any force fields supported by the [`SystemGenerator`](https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator) (e.g. `gaff-2.11`).
* The `--filepath <filepath>` argument specifies the directory in which simulation files should be written

A typical calculation will take < 12 hours.

To run the hydration free energy calculation using the `EspalomaCharge` toolkit and `espaloma-am1bcc` charge method for caffeiene:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit EspalomaCharge --method espaloma-am1bcc --filepath espaloma
```
To run the hydration free energy calculation using the `AmberTools` toolkit and `am1bcc` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit AmberTools --method am1bcc --filepath ambertools
```
To run the hydration free energy calculation using the `OpenEye` toolkit and `am1bccelf10` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit OpenEye --method am1bccelf10 --filepath openeye
```

**NOTE:** Currently, simulations are not able to be resumed, and existing filepaths must be deleted before being re-run.

For debugging purposes, it is useful to reduce the number of iterations: `--niterations 100` should take only ~10 min.
Don't go below 25, which is the current interval for online analysis---the YAML file will not be written if you go below 25 iterations.

### Running FreeSolv calculations

To run the FreeSolv database, first download the [FreeSolv database](https://github.com/MobleyLab/FreeSolv):
```bash
$ wget https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.json -O freesolv.json
```
To run molecule `$LSB_JOBINDEX`, you can use the `freesolv` command:
```bash
$ python hydration.py freesolv --index $LSB_JOBINDEX --toolkit EspalomaCharge --method espaloma-am1bcc --forcefield "gaff-2.11" --filepath espaloma --niterations 1000
```
A set of scripts `submit-{espaloma,ambertools,openeye}.sh` have been provided as an example of how to run the whole FreeSolv set via the LSF batch queue system.

### Analyzing FreeSolv calculations

To analyze all simulations and generate a `.csv` file of all compounds and a plot PDF file:
```bash
$ python hydration.py analyze-freesolv --filepath espaloma --outfile espaloma.csv --label espaloma
```
