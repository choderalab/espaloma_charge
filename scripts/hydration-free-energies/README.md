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
Download the [FreeSolv database](https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt):
```bash
$ wget https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt -O freesolv.csv
```

## Computing hydration free energies

### Running the simulation

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

A typical calculation will take ~1 hour.

To run the hydration free energy calculation using the `EspalomaCharge` toolkit and `espaloma-am1bcc` charge method for caffeiene:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit EspalomaCharge --method espaloma-am1bcc --filepath caffeine-espaloma
```
To run the hydration free energy calculation using the `AmberTools` toolkit and `am1bcc` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit AmberTools --method am1bcc --filepath caffeine-ambertools
```
To run the hydration free energy calculation using the `OpenEye` toolkit and `am1bccelf10` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit OpenEye --method am1bccelf10 --filepath caffeine-openeye
```

**NOTE:** Currently, simulations are not able to be resumed, and existing filepaths must be deleted before being re-run.

### Analyzing the simulation

To analyze the simulation, use 
```bash
$ python hydration.py analyze --filepath <filepath>
```
The `analyze` method will need to be modified to write data to a file or in an easily collectable manner.

## Manifest
* `freesolv.csv` - FreeSolv dataset ([retrieved 2022-11-29](https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt) and renamed)
* `hydration.py` - compute hydration free energy for the specified compound

