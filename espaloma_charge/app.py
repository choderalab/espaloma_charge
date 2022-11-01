import os
import tempfile
from urllib import request
import torch
from torch.utils.model_zoo import load_url
import numpy as np
from .utils import from_rdkit_mol
MODEL_URL = """
https://github.com/choderalab/espaloma_charge/releases/download/v0.0.0/model.pt
"""

def charge(
        molecule,
        total_charge: float = 0.0,
        model_url: str = MODEL_URL,
    ) -> np.ndarray:
    """Assign machine-learned AM1-BCC partial charges to a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        Input molecule.

    total_charge : float = 0.0


    model_url : str
        Model url.

    Returns
    -------
    np.ndarray : (n_atoms, ) array of partial charges.

    """
    with tempfile.TemporaryDirectory() as tempdir:
        target_path = os.path.join(tempdir, "model.pt")
        request.urlretrieve(model_url, target_path)
        model = torch.load(target_path, map_location="cpu")
    graph = from_rdkit_mol(molecule)
    graph = model(graph)
    return graph.ndata["q"].cpu().detach().flatten().numpy()