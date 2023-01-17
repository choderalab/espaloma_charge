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
from typing import Sequence


# TODO: Do we really want to define this at file level, 
# rather than within some kind of class?
MODEL_URL = """
https://github.com/choderalab/espaloma_charge/releases/download/v0.0.8/model.pt
"""

MODEL_PATH = ".espaloma_charge_model.pt"


def charge(
        molecule,
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
    # switch to `charge_multiple` if a sequence is fed
    if isinstance(molecule, Sequence):
        return charge_multiple(molecule)

    # load model
    if model_url is None:
        model_url = MODEL_URL

    if not os.path.exists(MODEL_PATH):
        request.urlretrieve(model_url, MODEL_PATH)

    model = torch.load(MODEL_PATH)

    # deploy the model to get charge
    graph = from_rdkit_mol(molecule)

    if torch.cuda.is_available():
        graph = graph.to("cuda:0")
        model = model.cuda()

    graph = model(graph)
    return graph.ndata["q"].cpu().detach().flatten().numpy()


def charge_multiple(
        molecules,
        model_url: str = None,
    ) -> np.ndarray:
    """Assign machine-learned AM1-BCC partial charges to a molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.Mol
        Input molecule.

    model_url : str, optional, default=None
        URL or filepath to retrieve the model from.
        If None, the default MODEL_URL (defined at file level) is used

    Returns
    -------
    np.ndarray : (n_atoms, ) array of partial charges.

    """
    if model_url is None:
        model_url = MODEL_URL

    if not os.path.exists(MODEL_PATH):
        request.urlretrieve(model_url, MODEL_PATH)

    model = torch.load(MODEL_PATH)

    graphs = [from_rdkit_mol(molecule) for molecule in molecules]
    graph = dgl.batch(graphs)
    
    if torch.cuda.is_available():
        graph = graph.to("cuda:0")
        model = model.cuda()

    graph = model(graph)
    graph = graph.to("cpu")
    graphs = dgl.unbatch(graph)

    return [graph.ndata["q"].detach().flatten().numpy() for graph in graphs]

