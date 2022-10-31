import torch
import numpy as np
from .utils import from_rdkit_mol

def charge(
        molecule,
        total_charge: float = 0.0,
        model_url: str = "",
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
    model = torch.utils.model_zoo.load_url(model_url)
    graph = from_rdkit_mol(molecule)
    graph = model(graph)
    return graph.ndata["q"].cpu().detach().flatten()
