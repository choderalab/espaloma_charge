"""Neural network components of espaloma charge."""

import torch

class _Sequential(torch.nn.Module):
    """Sequentially staggered neural networks."""

    def __init__(
        self,
        layer,
        config,
        in_features,
        model_kwargs={},
    ):
        super(_Sequential, self).__init__()

        self.exes = []

        # init dim
        dim = in_features

        # parse the config
        for idx, exe in enumerate(config):

            try:
                exe = float(exe)

                if exe >= 1:
                    exe = int(exe)
            except BaseException:
                pass

            # int -> feedfoward
            if isinstance(exe, int):
                setattr(self, "d" + str(idx), layer(dim, exe, **model_kwargs))

                dim = exe
                self.exes.append("d" + str(idx))

            # str -> activation
            elif isinstance(exe, str):
                if exe == "bn":
                    setattr(self, "a" + str(idx), torch.nn.BatchNorm1d(dim))

                else:
                    activation = getattr(torch.nn.functional, exe)
                    setattr(self, "a" + str(idx), activation)

                self.exes.append("a" + str(idx))

            # float -> dropout
            elif isinstance(exe, float):
                dropout = torch.nn.Dropout(exe)
                setattr(self, "o" + str(idx), dropout)

                self.exes.append("o" + str(idx))

    def forward(self, g, x):
        for exe in self.exes:
            if exe.startswith("d"):
                if g is not None:
                    x = getattr(self, exe)(g, x)
                else:
                    x = getattr(self, exe)(x)
            else:
                x = getattr(self, exe)(x)

        return x


class Sequential(torch.nn.Module):
    """Sequential neural network with input layers.

    Parameters
    ----------
    layer : torch.nn.Module
        DGL graph convolution layers.

    config : List
        A sequence of numbers (for units) and strings (for activation functions)
        denoting the configuration of the sequential model.

    feature_units : int(default=117)
        The number of input channels.

    Methods
    -------
    forward(g, x)
        Forward pass.
    """

    def __init__(
        self,
        layer,
        config,
        feature_units=117,
        input_units=128,
        model_kwargs={},
    ):
        super(Sequential, self).__init__()

        # initial featurization
        self.f_in = torch.nn.Sequential(
            torch.nn.Linear(feature_units, input_units), torch.nn.Tanh()
        )

        self._sequential = _Sequential(
            layer, config, in_features=input_units, model_kwargs=model_kwargs
        )

    def _forward(self, g, x):
        """Forward pass with graph and features."""
        for exe in self.exes:
            if exe.startswith("d"):
                x = getattr(self, exe)(g, x)
            else:
                x = getattr(self, exe)(x)

        return x

    def forward(self, g, x=None):
        """Forward pass.

        Parameters
        ----------
        g : `dgl.DGLHeteroGraph`,
            input graph

        Returns
        -------
        g : `dgl.DGLHeteroGraph`
            output graph
        """
        import dgl

        if x is None:
            # get node attributes
            x = g.ndata["h0"]
            x = self.f_in(x)

        # message passing on homo graph
        x = self._sequential(g, x)

        # put attribute back in the graph
        g.ndata["h"] = x

        return g

def get_charges(node):
    """ Solve the function to get the absolute charges of atoms in a
    molecule from parameters.
    Parameters
    ----------
    e : tf.Tensor, dtype = tf.float32,
        electronegativity.
    s : tf.Tensor, dtype = tf.float32,
        hardness.
    Q : tf.Tensor, dtype = tf.float32, shape=(),
        total charge of a molecule.
    We use Lagrange multipliers to analytically give the solution.
    $$
    U({\bf q})
    &= \sum_{i=1}^N \left[ e_i q_i +  \frac{1}{2}  s_i q_i^2\right]
        - \lambda \, \left( \sum_{j=1}^N q_j - Q \right) \\
    &= \sum_{i=1}^N \left[
        (e_i - \lambda) q_i +  \frac{1}{2}  s_i q_i^2 \right
        ] + Q
    $$
    This gives us:
    $$
    q_i^*
    &= - e_i s_i^{-1}
    + \lambda s_i^{-1} \\
    &= - e_i s_i^{-1}
    + s_i^{-1} \frac{
        Q +
         \sum\limits_{i=1}^N e_i \, s_i^{-1}
        }{\sum\limits_{j=1}^N s_j^{-1}}
    $$
    """
    e = node.data["e"]
    s = node.data["s"]
    sum_e_s_inv = node.data["sum_e_s_inv"]
    sum_s_inv = node.data["sum_s_inv"]
    sum_q = node.data["sum_q"]

    return {
        "q": -e * s**-1
        + (s**-1) * torch.div(sum_q + sum_e_s_inv, sum_s_inv)
    }

class ChargeReadout(torch.nn.Module):
    def __init__(self, in_features):
        super().__init__()
        self.fc_params = torch.nn.Linear(in_features, 2)

    def forward(self, g):
        h = self.fc_params(g.ndata["h"])
        e, s = h.split(1, -1)
        g.ndata["e"], g.ndata["s"] = e, s
        return g

class ChargeEquilibrium(torch.nn.Module):
    """Charge equilibrium within batches of molecules."""

    def __init__(self):
        super(ChargeEquilibrium, self).__init__()

    def forward(self, g, total_charge=0.0):
        """apply charge equilibrium to all molecules in batch"""
        # calculate $s ^ {-1}$ and $ es ^ {-1}$
        import dgl

        g.apply_nodes(
            lambda node: {"s_inv": node.data["s"] ** -1},
        )

        g.apply_nodes(
            lambda node: {"e_s_inv": node.data["e"] * node.data["s"] ** -1},
        )

        if "q_ref" in g.ndata:
            total_charge = dgl.sum_nodes(g, "q_ref")
        else:
            total_charge = torch.ones(1, 1) * total_charge
        
        g.ndata["sum_q"] = dgl.broadcast_nodes(g, total_charge)

        sum_s_inv = dgl.sum_nodes(g, "s_inv")
        sum_e_s_inv = dgl.sum_nodes(g, "e_s_inv")
        g.ndata["sum_s_inv"] = dgl.broadcast_nodes(g, sum_s_inv)
        g.ndata["sum_e_s_inv"] = dgl.broadcast_nodes(g, sum_e_s_inv)

        # g.update_all(
        #     dgl.function.copy_src(src="sum_q", out="m_sum_q"),
        #     dgl.function.sum(msg="m_sum_q", out="sum_q"),
        #     etype="g_has_n1",
        # )
        #
        # # get the sum of $s^{-1}$ and $m_s^{-1}$
        # g.update_all(
        #     dgl.function.copy_src(src="s_inv", out="m_s_inv"),
        #     dgl.function.sum(msg="m_s_inv", out="sum_s_inv"),
        #     etype="n1_in_g",
        # )
        #
        # g.update_all(
        #     dgl.function.copy_src(src="e_s_inv", out="m_e_s_inv"),
        #     dgl.function.sum(msg="m_e_s_inv", out="sum_e_s_inv"),
        #     etype="n1_in_g",
        # )
        #
        # g.update_all(
        #     dgl.function.copy_src(src="sum_s_inv", out="m_sum_s_inv"),
        #     dgl.function.sum(msg="m_sum_s_inv", out="sum_s_inv"),
        #     etype="g_has_n1",
        # )
        #
        # g.update_all(
        #     dgl.function.copy_src(src="sum_e_s_inv", out="m_sum_e_s_inv"),
        #     dgl.function.sum(msg="m_sum_e_s_inv", out="sum_e_s_inv"),
        #     etype="g_has_n1",
        # )

        g.apply_nodes(get_charges)

        return g
