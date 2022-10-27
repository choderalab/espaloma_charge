import dgl
from . import utils

class Graph(object):
    """A unified graph object that support translation to and from
    message-passing graphs and MM factor graph.

    Methods
    -------
    save(path)
        Save graph to file.

    load(path)
        Load a graph from path.

    Note
    ----
    This object provides access to popular attributes of homograph and
    heterograph.

    This object also provides access to `ndata` and `edata` from the heterograph.

    Examples
    --------


    """

    def __init__(self, mol=None):
        homograph = self.get_homograph_from_mol(mol)

        self.mol = mol
        self.homograph = homograph

    @staticmethod
    def get_homograph_from_mol(mol):
        graph = utils.from_rdkit_mol(mol)
        return graph

    @property
    def ndata(self):
        return self.homograph.ndata

    @property
    def edata(self):
        return self.homograph.edata

    @property
    def nodes(self):
        return self.heterograph.nodes
