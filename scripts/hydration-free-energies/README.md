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
```
```

## Manifest
* `freesolv.csv` - FreeSolv dataset ([retrieved 2022-11-29](https://raw.githubusercontent.com/MobleyLab/FreeSolv/master/database.txt) and renamed)
* `hydration.py` - compute hydration free energy for the specified compound

