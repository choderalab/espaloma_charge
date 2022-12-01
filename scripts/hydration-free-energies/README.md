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

The `python hydration.py run --smiles <smiles>` command will set up and run an alchemical hydration free energy calculation for the SMILES string `<smiles>.

The `--toolkit` and `--method` arguments can be used to specify the toolkit (`EspalomaCharge`, `AmberTools`, `OpenEye`, `RDKit`) and charge method (depends on toolkit).

To run the hydration free energy calculation using the `EspalomaCharge` toolkit and `espaloma-am1bcc` charge method for caffeiene:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit EspalomaCharge --method espaloma-am1bcc
```
To run the hydration free energy calculation using the `AmberTools` toolkit and `am1bcc` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit AmberTools --method am1bcc
```
To run the hydration free energy calculation using the `OpenEye` toolkit and `am1bccelf10` charge method:
```bash
$ python hydration.py run --smiles "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" --toolkit OpenEye --method am1bccelf10
```

The small molecule force field can be specified using the `--forcefield` flag.
All force fields supported by the [`SystemGenerator`](https://github.com/openmm/openmmforcefields#automating-force-field-management-with-systemgenerator) are supported.

**NOTE:** Existing NetCDF files will not be overwritten, but will instead be continued.

## Manifest
* `freesolv.csv` - FreeSolv dataset ([retrieved 2022-11-29](https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt) and renamed)
* `hydration.py` - compute hydration free energy for the specified compound

