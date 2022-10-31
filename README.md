# Espaloma Charge
Standalone charge assignment from Espaloma framework. https://doi.org/10.1039/D2SC02739A

## Installation

```
pip install espaloma_charge
```

## Example

```python
from rdkit import Chem; from espaloma_charge import charge
molecule = Chem.MolFromSmiles("N#N")
charge(molecule)
>>> array([0., 0.], dtype=float32)

```
