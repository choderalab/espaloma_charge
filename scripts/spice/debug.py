import pandas as pd
import torch
import dgl
from rdkit import Chem
from openff.toolkit.topology import Molecule
# dgl.use_libxsmm(False)

def run():
    from espaloma_charge.utils import from_rdkit_mol
    molecule = Molecule.from_smiles("CCl")
    molecule.assign_partial_charges("am1bcc")
    print(molecule.partial_charges)
    g = from_rdkit_mol(molecule.to_rdkit())
    model = torch.load("model.pt")
    model(g)
    print(g.ndata['q'])


if __name__ == "__main__":
    run()
