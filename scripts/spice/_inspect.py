import pandas as pd
import torch
import dgl
from openff.toolkit.topology import Molecule
# dgl.use_libxsmm(False)
class ChargeDataset(dgl.data.DGLDataset):
    def __init__(self, graphs):
        super().__init__(name="charge_dataset")
        self.graphs = graphs

    def __len__(self):
        return len(self.graphs)

    def __getitem__(self, idx):
        return self.graphs[idx]

def run():

    from espaloma_charge.utils import from_rdkit_mol
    molecules = Molecule.from_file("spice.oeb", allow_undefined_stereo=True)

    import pandas as pd
    df = pd.DataFrame(columns=["SMILES"])
    for molecule in molecules:
        df = df.append(
            {
                'SMILES': molecule.to_smiles(),
            },
            ignore_index=True,
        )


    df.to_csv("inspect.csv")
    from rdkit.Chem import PandasTools
    PandasTools.AddMoleculeColumnToFrame(df, "SMILES", "MOL")
    open("inspect_all.html", "w").write(df.to_html())



if __name__ == "__main__":
    run()
