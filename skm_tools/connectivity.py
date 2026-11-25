'''Are a set of nodes (e.g. differentially expressed genes) more connected than random nodes?

The test compares the subnetwork induced by the nodes with the subnetworks induced by random node sets of the same
size and degree distribution, drawn from a universe (e.g. the genes measured in the experiment). Statistics: the
number of edges among the nodes, and the size of the largest connected component (LCC).
'''

import logging

import numpy as np

logger = logging.getLogger(__name__)


def _components_largest(n, u, v):
    '''Size of the largest connected component of a graph with nodes 0..n-1 and edges (u, v) (numpy arrays).'''
    if n == 0:
        return 0
    lab = np.arange(n)
    while True:
        m = np.minimum(lab[u], lab[v])
        new = lab.copy()
        np.minimum.at(new, u, m)
        np.minimum.at(new, v, m)
        new = new[new]                       # pointer jumping: follow the label to its own label
        if np.array_equal(new, lab):
            break
        lab = new
    return int(np.bincount(lab).max())


class _Pairs:
    '''Undirected node pairs of a graph restricted to a universe, as index arrays.'''

    def __init__(self, g, universe):
        self.nodes = list(universe)
        self.index = {x: i for i, x in enumerate(self.nodes)}
        pairs = {(min(a, b), max(a, b)) for a, b in
                 ((self.index[x], self.index[y]) for x, y in g.edges() if x in self.index and y in self.index)
                 if a != b}
        arr = np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)
        self.u, self.v = arr[:, 0], arr[:, 1]
        self.degree = np.bincount(np.concatenate([self.u, self.v]), minlength=len(self.nodes))

    def stats(self, idx):
        '''Number of edges and LCC size among the node indices idx.'''
        mask = np.zeros(len(self.nodes), bool)
        mask[idx] = True
        keep = mask[self.u] & mask[self.v]
        local = np.full(len(self.nodes), -1)
        local[idx] = np.arange(len(idx))
        return int(keep.sum()), _components_largest(len(idx), local[self.u[keep]], local[self.v[keep]])


def induced_connectivity(g, nodes, universe=None, n_random=500, n_bins=10, seed=None):
    '''Are `nodes` more connected among themselves in `g` than random nodes?

    `g` is treated as undirected: edge directions, parallel edges (e.g. a CKN binding and a regulation edge
    between the same genes) and self-loops count once / not at all, so the statistics are about node pairs.
    Random node sets have the same size as `nodes` and the same degree distribution: the universe's nodes are
    split into `n_bins` quantile bins of their degree (within the universe), and each random set draws as many
    nodes from each bin as `nodes` has there.

    Parameters
    ----------
    g : networkx.Graph
        The network, e.g. CKN filtered to the edge ranks of interest.
    nodes : iterable
        The nodes to test, e.g. the DEGs. Nodes that are not in `g` or not in the universe are left out.
    universe : iterable, optional
        The nodes the random sets are drawn from, e.g. the genes measured in the experiment (default: all nodes of
        `g`). Only nodes of `g` count.
    n_random : int
        Number of random node sets (default 500). The smallest possible p-value is 1 / (n_random + 1).
    n_bins : int
        Number of degree bins (default 10).
    seed : int or numpy.random.Generator, optional
        For reproducible random sets.

    Returns
    -------
    dict
        ``n_nodes`` (nodes tested), ``n_universe``; for ``edges`` and ``lcc``: the observed value, the mean over
        the random sets (``edges_random``, ``lcc_random``), the ratio (``edges_fold``, ``lcc_fold``) and the
        empirical p-value, (number of random sets at least as high + 1) / (n_random + 1) (``edges_p``,
        ``lcc_p``). With fewer than 2 nodes only ``n_nodes`` and ``n_universe`` are returned.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.complete_graph(5)
    >>> g.add_edges_from((i, i + 1) for i in range(5, 40))
    >>> r = induced_connectivity(g, [0, 1, 2, 3, 4], n_random=99, seed=1)
    >>> r["edges"], r["lcc"]
    (10, 5)
    '''
    universe = set(g) if universe is None else set(universe) & set(g)
    pairs = _Pairs(g, sorted(universe, key=str))
    sel = np.array(sorted({pairs.index[x] for x in nodes if x in pairs.index}), dtype=np.int64)
    res = dict(n_nodes=len(sel), n_universe=len(pairs.nodes))
    if len(sel) < 2:
        logger.info("Fewer than 2 of the nodes are in the universe: nothing to test.")
        return res
    rng = np.random.default_rng(seed)
    deg = pairs.degree
    edges_q = np.unique(np.quantile(deg, np.linspace(0, 1, n_bins + 1))[1:-1])
    bins = np.digitize(deg, edges_q)
    by_bin = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    need = {b: int(k) for b, k in zip(*np.unique(bins[sel], return_counts=True))}
    obs = pairs.stats(sel)
    null = np.array([pairs.stats(np.concatenate([rng.choice(by_bin[b], k, replace=False) for b, k in need.items()]))
                     for _ in range(n_random)])
    for i, name in enumerate(("edges", "lcc")):
        mean = null[:, i].mean()
        res.update({name: obs[i], f"{name}_random": float(mean),
                    f"{name}_fold": float(obs[i] / mean) if mean > 0 else float("inf") if obs[i] else float("nan"),
                    f"{name}_p": float(((null[:, i] >= obs[i]).sum() + 1) / (n_random + 1))})
    logger.info("%d nodes: %d edges (random %.1f), LCC %d (random %.1f)", len(sel), obs[0], res["edges_random"],
                obs[1], res["lcc_random"])
    return res
