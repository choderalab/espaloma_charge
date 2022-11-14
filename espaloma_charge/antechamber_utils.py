from typing import Optional
import numpy as np
from rdkit import Chem
from .app import charge

def charge_mol2(in_path: str, out_path: Optional[str]) -> None:
    """Assign charges to mol2 and write crg file.

    Parameters
    ----------
    in_path : str
        Input path.

    out_path : str
        Output path.

    """
    assert "mol2" in in_path, "Only handles mol2."
    if len(out_path) == 0:
        out_path = in_path.replace(".mol2", ".crg")
    molecule = Chem.MolFromMol2File(in_path)
    charges = charge(molecule)
    out = " " + np.array2string(
        charges, precision=6, max_line_width=81, floatmode="fixed"
    )[1:-1]
    with open(out_path, "w") as out_file:
        out_file.write(out)
