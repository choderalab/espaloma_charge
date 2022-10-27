import torch

def fp_rdkit(atom):
    from rdkit import Chem

    HYBRIDIZATION_RDKIT = {
        Chem.rdchem.HybridizationType.SP: torch.tensor(
            [1, 0, 0, 0, 0],
            dtype=torch.get_default_dtype(),
        ),
        Chem.rdchem.HybridizationType.SP2: torch.tensor(
            [0, 1, 0, 0, 0],
            dtype=torch.get_default_dtype(),
        ),
        Chem.rdchem.HybridizationType.SP3: torch.tensor(
            [0, 0, 1, 0, 0],
            dtype=torch.get_default_dtype(),
        ),
        Chem.rdchem.HybridizationType.SP3D: torch.tensor(
            [0, 0, 0, 1, 0],
            dtype=torch.get_default_dtype(),
        ),
        Chem.rdchem.HybridizationType.SP3D2: torch.tensor(
            [0, 0, 0, 0, 1],
            dtype=torch.get_default_dtype(),
        ),
        Chem.rdchem.HybridizationType.S: torch.tensor(
            [0, 0, 0, 0, 0],
            dtype=torch.get_default_dtype(),
        ),
    }
    return torch.cat(
        [
            torch.tensor(
                [
                    atom.GetTotalDegree(),
                    atom.GetTotalValence(),
                    atom.GetExplicitValence(),
                    atom.GetFormalCharge(),
                    atom.GetIsAromatic() * 1.0,
                    atom.GetMass(),
                    atom.IsInRingSize(3) * 1.0,
                    atom.IsInRingSize(4) * 1.0,
                    atom.IsInRingSize(5) * 1.0,
                    atom.IsInRingSize(6) * 1.0,
                    atom.IsInRingSize(7) * 1.0,
                    atom.IsInRingSize(8) * 1.0,
                ],
                dtype=torch.get_default_dtype(),
            ),
            HYBRIDIZATION_RDKIT[atom.GetHybridization()],
        ],
        dim=0,
    )


def from_rdkit_mol(mol, use_fp=True):
    import dgl
    from rdkit import Chem

    # initialize graph
    g = dgl.DGLGraph()

    # enter nodes
    n_atoms = mol.GetNumAtoms()
    g.add_nodes(n_atoms)
    g.ndata["type"] = torch.Tensor(
        [[atom.GetAtomicNum()] for atom in mol.GetAtoms()]
    )
    g.ndata["formal_charge"] = torch.Tensor(
        [[atom.GetFormalCharge()] for atom in mol.GetAtoms()]
    )
    h_v = torch.zeros(g.ndata["type"].shape[0], 100, dtype=torch.float32)

    h_v[
        torch.arange(g.ndata["type"].shape[0]),
        torch.squeeze(g.ndata["type"]).long(),
    ] = 1.0

    h_v_fp = torch.stack([fp_rdkit(atom) for atom in mol.GetAtoms()], axis=0)

    if use_fp == True:
        h_v = torch.cat([h_v, h_v_fp], dim=-1)  # (n_atoms, 117)

    g.ndata["h0"] = h_v

    # enter bonds
    bonds = list(mol.GetBonds())
    bonds_begin_idxs = [bond.GetBeginAtomIdx() for bond in bonds]
    bonds_end_idxs = [bond.GetEndAtomIdx() for bond in bonds]
    bonds_types = [bond.GetBondType().real for bond in bonds]

    # NOTE: dgl edges are directional
    g.add_edges(bonds_begin_idxs, bonds_end_idxs)
    g.add_edges(bonds_end_idxs, bonds_begin_idxs)

    # g.edata["type"] = torch.Tensor(bonds_types)[:, None].repeat(2, 1)

    return g


def from_homogeneous_and_mol(g, offmol):
    r"""Build heterogeneous graph from homogeneous ones.
    Note
    ----
    For now we name single node, two-, three, and four-,
    hypernodes as `n1`, `n2`, `n3`, and `n4`. These correspond
    to atom, bond, angle, and torsion nodes in chemical graphs.
    Parameters
    ----------
    g : `espaloma.HomogeneousGraph` object
        the homogeneous graph to be translated.
    Returns
    -------
    hg : `espaloma.HeterogeneousGraph` object
        the resulting heterogeneous graph.
    """

    # initialize empty dictionary
    hg = {}

    # get adjacency matrix
    a = g.adjacency_matrix()

    # get all the indices
    idxs = relationship_indices_from_offmol(offmol)

    # make them all numpy
    idxs = {key: value.numpy() for key, value in idxs.items()}

    # also include n1
    idxs["n1"] = np.arange(g.number_of_nodes())[:, None]

    # =========================
    # neighboring relationships
    # =========================
    # NOTE:
    # here we only define the neighboring relationship
    # on atom level
    hg[("n1", "n1_neighbors_n1", "n1")] = idxs["n2"]

    # build a mapping between indices and the ordering
    idxs_to_ordering = {}

    for term in ["n1", "n2", "n3", "n4", "n4_improper"]:
        idxs_to_ordering[term] = {
            tuple(subgraph_idxs): ordering
            for (ordering, subgraph_idxs) in enumerate(list(idxs[term]))
        }

    # ===============================================
    # relationships between nodes of different levels
    # ===============================================
    # NOTE:
    # here we define all the possible
    # 'has' and 'in' relationships.
    # TODO:
    # we'll test later to see if this adds too much overhead
    #

    for small_idx in range(1, 5):
        for big_idx in range(small_idx + 1, 5):
            for pos_idx in range(big_idx - small_idx + 1):

                hg[
                    (
                        "n%s" % small_idx,
                        "n%s_as_%s_in_n%s" % (small_idx, pos_idx, big_idx),
                        "n%s" % big_idx,
                    )
                ] = np.stack(
                    [
                        np.array(
                            [
                                idxs_to_ordering["n%s" % small_idx][tuple(x)]
                                for x in idxs["n%s" % big_idx][
                                    :, pos_idx : pos_idx + small_idx
                                ]
                            ]
                        ),
                        np.arange(idxs["n%s" % big_idx].shape[0]),
                    ],
                    axis=1,
                )

                hg[
                    (
                        "n%s" % big_idx,
                        "n%s_has_%s_n%s" % (big_idx, pos_idx, small_idx),
                        "n%s" % small_idx,
                    )
                ] = np.stack(
                    [
                        np.arange(idxs["n%s" % big_idx].shape[0]),
                        np.array(
                            [
                                idxs_to_ordering["n%s" % small_idx][tuple(x)]
                                for x in idxs["n%s" % big_idx][
                                    :, pos_idx : pos_idx + small_idx
                                ]
                            ]
                        ),
                    ],
                    axis=1,
                )

    # ======================================
    # nonbonded terms
    # ======================================
    # NOTE: everything is counted twice here
    # nonbonded is where
    # $A = AA = AAA = AAAA = 0$

    # make dense
    a_ = a.to_dense().detach().numpy()

    idxs["nonbonded"] = np.stack(
        np.where(np.equal(a_ + a_ @ a_ + a_ @ a_ @ a_, 0.0)),
        axis=-1,
    )

    # onefour is the two ends of torsion
    # idxs["onefour"] = np.stack(
    #     [
    #         idxs["n4"][:, 0],
    #         idxs["n4"][:, 3],
    #     ],
    #     axis=1,
    # )

    idxs["onefour"] = np.stack(
        np.where(
            np.equal(a_ + a_ @ a_, 0.0) * np.greater(a_ @ a_ @ a_, 0.0),
        ),
        axis=-1,
    )

    # membership
    for term in ["nonbonded", "onefour"]:
        for pos_idx in [0, 1]:
            hg[(term, "%s_has_%s_n1" % (term, pos_idx), "n1")] = np.stack(
                [np.arange(idxs[term].shape[0]), idxs[term][:, pos_idx]],
                axis=-1,
            )

            hg[("n1", "n1_as_%s_in_%s" % (pos_idx, term), term)] = np.stack(
                [
                    idxs[term][:, pos_idx],
                    np.arange(idxs[term].shape[0]),
                ],
                axis=-1,
            )

    # membership of n1 in n4_improper
    for term in ["n4_improper"]:
        for pos_idx in [0, 1, 2, 3]:
            hg[(term, "%s_has_%s_n1" % (term, pos_idx), "n1")] = np.stack(
                [np.arange(idxs[term].shape[0]), idxs[term][:, pos_idx]],
                axis=-1,
            )

            hg[("n1", "n1_as_%s_in_%s" % (pos_idx, term), term)] = np.stack(
                [
                    idxs[term][:, pos_idx],
                    np.arange(idxs[term].shape[0]),
                ],
                axis=-1,
            )

    # ======================================
    # relationships between nodes and graphs
    # ======================================
    for term in [
        "n1",
        "n2",
        "n3",
        "n4",
        "n4_improper",
        "nonbonded",
        "onefour",
    ]:
        hg[(term, "%s_in_g" % term, "g",)] = np.stack(
            [np.arange(len(idxs[term])), np.zeros(len(idxs[term]))],
            axis=1,
        )

        hg[("g", "g_has_%s" % term, term)] = np.stack(
            [
                np.zeros(len(idxs[term])),
                np.arange(len(idxs[term])),
            ],
            axis=1,
        )

    import dgl

    hg = dgl.heterograph(
        {key: value.astype(np.int32).tolist() for key, value in hg.items()}
    )

    hg.nodes["n1"].data["h0"] = g.ndata["h0"]
    hg.nodes["g"].data["sum_q"] = g.ndata["sum_q"][0].reshape(1, 1)
    # include indices in the nodes themselves
    for term in [
        "n1",
        "n2",
        "n3",
        "n4",
        "n4_improper",
        "onefour",
        "nonbonded",
    ]:
        hg.nodes[term].data["idxs"] = torch.tensor(idxs[term])

    return hg
