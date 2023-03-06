import torch

from dgllife.utils.featurizers import (
    BaseAtomFeaturizer,
    ConcatFeaturizer,
    atom_type_one_hot,
    atom_degree_one_hot,
    atom_hybridization_one_hot,
    atom_is_aromatic,
    atom_total_num_H_one_hot
)


def atom_ring_size_one_hot(atom):
    return [                    
            atom.IsInRingSize(3),
            atom.IsInRingSize(4),
            atom.IsInRingSize(5),
            atom.IsInRingSize(6),
            atom.IsInRingSize(7),
            atom.IsInRingSize(8),
    ]

class AtomFeaturizer(BaseAtomFeaturizer):
    def __init__(self, atom_data_field='h'):
        super().__init__(
            featurizer_funcs={atom_data_field: ConcatFeaturizer(
                [
                 atom_type_one_hot,
                 atom_degree_one_hot,
                 # atom_implicit_valence_one_hot,
                 # atom_formal_charge,
                 # atom_num_radical_electrons,
                 atom_hybridization_one_hot,
                 atom_is_aromatic,
                 atom_total_num_H_one_hot,
                 atom_ring_size_one_hot,
                 ]
            )})





def from_rdkit_mol(mol):
    import dgl
    from rdkit import Chem
    from dgllife.utils import mol_to_bigraph

    # g = CanonicalAtomFeaturizer("h0")(mol)
    g = mol_to_bigraph(
        mol, add_self_loop=False, node_featurizer=AtomFeaturizer("h0"),
        canonical_atom_order=False,
    )

    # # initialize graph
    # g = dgl.DGLGraph()

    # # enter nodes
    # n_atoms = mol.GetNumAtoms()
    # g.add_nodes(n_atoms)
    # g.ndata["type"] = torch.Tensor(
    #     [[atom.GetAtomicNum()] for atom in mol.GetAtoms()]
    # )

    g.ndata["q_ref"] = torch.Tensor(
        [[atom.GetFormalCharge()] for atom in mol.GetAtoms()]
    )
    # h_v = torch.zeros(g.ndata["type"].shape[0], 100, dtype=torch.float32)

    # h_v[
    #     torch.arange(g.ndata["type"].shape[0]),
    #     torch.squeeze(g.ndata["type"]).long(),
    # ] = 1.0

    # h_v_fp = torch.stack([fp_rdkit(atom) for atom in mol.GetAtoms()], axis=0)

    # if use_fp == True:
    #     h_v = torch.cat([h_v, h_v_fp], dim=-1)  # (n_atoms, 117)

    # g.ndata["h0"] = h_v

    # # enter bonds
    # bonds = list(mol.GetBonds())
    # bonds_begin_idxs = [bond.GetBeginAtomIdx() for bond in bonds]
    # bonds_end_idxs = [bond.GetEndAtomIdx() for bond in bonds]
    # bonds_types = [bond.GetBondType().real for bond in bonds]

    # # NOTE: dgl edges are directional
    # g.add_edges(bonds_begin_idxs, bonds_end_idxs)
    # g.add_edges(bonds_end_idxs, bonds_begin_idxs)

    # g.edata["type"] = torch.Tensor(bonds_types)[:, None].repeat(2, 1)

    return g
