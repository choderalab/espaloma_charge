import os
import tempfile
from urllib import request
from rdkit import Chem
import torch
from torch.utils.model_zoo import load_url
import numpy as np
from .utils import from_rdkit_mol

try:
    import dgl
    dgl.use_libxsmm(False)
except:
    pass

# TODO: Do we really want to define this at file level, rather than within some kind of class?
MODEL_URL = """
https://github.com/choderalab/espaloma_charge/releases/download/v0.0.7/model.pt
"""

def charge(
        molecule,
        total_charge: float = None,
        model_url: str = None,
    ) -> np.ndarray:
    """Assign machine-learned AM1-BCC partial charges to a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        Input molecule.

    total_charge : float = 0.0


    model_url : str, optional, default=None
        URL or filepath to retrieve the model from.
        If None, the default MODEL_URL (defined at file level) is used

    Returns
    -------
    np.ndarray : (n_atoms, ) array of partial charges.

    """
    if model_url is None:
        model_url = MODEL_URL

    if not os.path.exists(".model.pt"):
        request.urlretrieve(model_url, ".model.pt")

    model = torch.load(".model.pt")

    if total_charge is None:
        total_charge = Chem.GetFormalCharge(molecule)
    graph = from_rdkit_mol(molecule)
    graph = model(graph)
    return graph.ndata["q"].cpu().detach().flatten().numpy()
