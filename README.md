Espaloma Charge
=======

Standalone charge assignment from Espaloma framework. https://doi.org/10.1039/D2SC02739A

## Installation

```bash
$ pip install espaloma_charge
```

## Examples
**Option0: Assign charges to rdkit molecule.**

```python
>>> from rdkit import Chem; from espaloma_charge import charge
>>> molecule = Chem.MolFromSmiles("N#N")
>>> charge(molecule)
array([0., 0.], dtype=float32)

```

Assign charges to your favorite molecule in 
[![Google Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1e14EkNyidPI0wXBGcewh9m9LC1imSRWZ?usp=sharing)


**Option1: Use with [`openff-toolkit`](https://github.com/openforcefield/openff-toolkit)**(installation required)

```python
>>> from openff.toolkit.topology import Molecule
>>> from espaloma_charge.openff_wrapper import EspalomaChargeToolkitWrapper
>>> toolkit_registry = EspalomaChargeToolkitWrapper()
>>> molecule = Molecule.from_smiles("N#N")
>>> molecule.assign_partial_charges('espaloma-am1bcc', toolkit_registry=toolkit_registry)
>>> molecule.partial_charges
<Quantity([0. 0.], 'elementary_charge')>
```

**Option2: Use as Command Line Interface to write [`antechamber`](http://ambermd.org/antechamber/ac.html)-compatible charges.**
```bash
$ espaloma_charge -i in.mol2 -o in.crg
$ antechamber -fi mol2 -fo mol2 -i in.mol2 -o out.mol2 -c rc -cf in.crg 
```

## Reference
If you are using this little tool in your pipeline, please consider citing:

```
@Article{D2SC02739A,
author ="Wang, Yuanqing and Fass, Josh and Kaminow, Benjamin and Herr, John E. and Rufa, Dominic and Zhang, Ivy and Pulido, Iván and Henry, Mike and Bruce Macdonald, Hannah E. and Takaba, Kenichiro and Chodera, John D.",
title  ="End-to-end differentiable construction of molecular mechanics force fields",
journal  ="Chem. Sci.",
year  ="2022",
volume  ="13",
issue  ="41",
pages  ="12016-12033",
publisher  ="The Royal Society of Chemistry",
doi  ="10.1039/D2SC02739A",
url  ="http://dx.doi.org/10.1039/D2SC02739A"}

```
