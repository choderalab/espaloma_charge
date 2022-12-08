import os
import tempfile
from urllib import request
from rdkit import Chem
import dgl
try:
    dgl.use_libxsmm(False)
except:
    pass
import torch
from torch.utils.model_zoo import load_url
import numpy as np
from .utils import from_rdkit_mol

# TODO: Do we really want to define this at file level, rather than within some kind of class?
MODEL_URL = """
https://github.com/choderalab/espaloma_charge/releases/download/v0.0.0/model.pt
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

    # TODO: Can we memoize/cache the model so we don't have to retrieve it every time we invoke charging?
    with tempfile.TemporaryDirectory() as tempdir:
        target_path = os.path.join(tempdir, "model.pt")
        if os.path.exists(model_url):
            # model_url is a local filepath
            import shutil
            shutil.copyfile(model_url, target_path)
        else:
            # See if model_url is a URL
            request.urlretrieve(model_url, target_path)
        model = torch.load(target_path, map_location="cpu")
    if total_charge is None:
        total_charge = Chem.GetFormalCharge(molecule)
    graph = from_rdkit_mol(molecule)
    n_nodes = graph.number_of_nodes()
    graph.ndata["q_ref"] = torch.ones(n_nodes, 1) * total_charge / n_nodes
    graph = model(graph)
    return graph.ndata["q"].cpu().detach().flatten().numpy()
