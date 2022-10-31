import pandas as pd
import torch
import dgl
from openff.toolkit.topology import Molecule
dgl.use_libxsmm(False)
class ChargeDataset(dgl.data.DGLDataset):
    def __init__(self, graphs):
        super().__init__(name="charge_dataset")
        self.graphs = graphs

    def __len__(self):
        return len(self.graphs)

    def __getitem__(self, idx):
        return self.graphs[idx]

def run(args):
    from espaloma_charge.utils import from_rdkit_mol
    molecules = Molecule.from_file(args.path, allow_undefined_stereo=True)
    graphs = [from_rdkit_mol(molecule.to_rdkit()) for molecule in molecules]
    charges = [molecule.partial_charges for molecule in molecules]
    for graph, charge in zip(graphs, charges):
        graph.ndata["q_ref"] = torch.tensor(charge.magnitude).unsqueeze(-1)

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
    dataset_train, dataset_valid, dataset_test = dgl.data.utils.split_dataset(dataset, random_state=2666)
    dataloader = dgl.dataloading.GraphDataLoader(dataset_train, batch_size=args.batch_size, pin_memory=True)

    dataloader_valid = dgl.dataloading.GraphDataLoader(dataset_valid, batch_size=args.batch_size)
    dataloader_test = dgl.dataloading.GraphDataLoader(dataset_test, batch_size=args.batch_size)
    

    rmse_vl_best = 9999.9
    rmse_te_best = 9999.9

    if torch.cuda.is_available(): model = model.cuda()
    optimizer = torch.optim.Adam(
        model.parameters(), 
        args.learning_rate,
        weight_decay=args.weight_decay,
    )

    for idx_epoch in range(args.n_epochs):
        for g in dataloader:
            optimizer.zero_grad()
            if torch.cuda.is_available():
                g = g.to("cuda:0")
            g = model(g)
            loss = torch.nn.MSELoss()(
                g.ndata["q_ref"],
                g.ndata["q"],
            )

            loss.backward()
            optimizer.step()


        q_vl = []
        q_vl_ref = []
        q_te = []
        q_te_ref = []
        for g in dataloader_valid:
            if torch.cuda.is_available():
                g = g.to("cuda:0")
                g = model(g)
            q_vl.append(g.ndata["q"])
            q_vl_ref.append(g.ndata["q_ref"])

        for g in dataloader_test:
            if torch.cuda.is_available():
                g = g.to("cuda:0")
            g = model(g)
            q_te.append(g.ndata["q"])
            q_te_ref.append(g.ndata["q_ref"])

        q_vl = torch.cat(q_vl)
        q_vl_ref = torch.cat(q_vl_ref)
        q_te = torch.cat(q_te)
        q_te_ref = torch.cat(q_te_ref)

        rmse_vl = (q_vl - q_vl_ref).pow(2).mean().pow(0.5).item()
        rmse_te = (q_te - q_te_ref).pow(2).mean().pow(0.5).item()

        if rmse_vl <= rmse_vl_best:
            rmse_vl_best = rmse_vl
            rmse_te_best = rmse_te
            torch.save(model.cpu(), "model.pt")
     
    import json
    import pandas as pd
    key = dict(vars(args))
    key.pop("out")
    key = json.dumps(key)

    df = pd.DataFrame.from_dict(
        {key: [rmse_vl_best, rmse_te_best]},
        orient="index",
        columns=["vl", "te"]
    )
    df.to_csv(args.out + ".csv", mode="a", header=False)




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--depth", type=int, default=4)
    parser.add_argument("--width", type=int, default=32)
    parser.add_argument("--activation", type=str, default="relu")
    parser.add_argument("--path", type=str, default="spice.mol2")
    parser.add_argument("--learning_rate", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--n_epochs", type=int, default=5000)
    parser.add_argument("--batch_size", type=int, default=32)
    parser.add_argument("--out", type=str, default="_out")
    args = parser.parse_args()
    run(args)
