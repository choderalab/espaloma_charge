import torch
import dgl
from openff.toolkit.topology import Molecule

class ChargeDataset(dgl.data.DGLDataset):
    def __init__(self, graphs):
        super().__init__()
        self.graphs = graphs

    def __len__(self):
        return len(self.graphs)

    def __getitem__(self, idx):
        return self.graphs[idx]

def run(args):
    from espaloma_charge.utils import from_rdkit_mol
    molecules = Molecule.from_files(args.path)
    graphs = [from_rdkit_mol(molecule.to_rdkit()) for molecule in molecules]
    charges = [molecule.partial_charge for molecule in molecules]
    for graph, charge in zip(graphs, charges):
        graph.ndata["q_ref"] = torch.tensor(charge._value)

    config = [args.width, args.activation] * args.depth
    from espaloma_charge.models import (
        Sequential, ChargeReadout, ChargeEquilibrium
    )
    from functools import partial
    model = torch.nn.Sequential(
        Sequential(
            layer=partial(dgl.nn.SAGEConv, aggregator_type="mean"),
            config=config,
        ),
        ChargeReadout(args.width),
        ChargeEquilibrium(),
    )

    dataset = ChargeDataset(graphs)
    dataloader = dgl.dataloading.GraphDataloader(dataset, batch_size=32)

    for g in dataloader:
        print(g)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--width", type=int, default=32)
    parser.add_argument("--activation", type=str, default="relu")
    parser.add_argument("--path", type=str, default="spice.mol2")
    args = parser.parse_arguments()
    run(args)
