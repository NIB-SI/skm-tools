'''Export networks as DiNAR node and edge tables.

`DiNAR <https://github.com/NIB-SI/DiNAR>`_ (Differential Network Analysis in R) shows
experimental data on a knowledge network, one cluster of nodes at a time. Its custom
networks are tab-separated node and edge tables with fixed column names (as DiNAR's
``Networks/PSS/Ath-based_*.txt``); missing values are ``-``.
'''

import math

import networkx as nx
import pandas as pd


_MISSING = "-"

DINAR_NODE_COLUMNS = [
    "geneID", "shortDescription", "shortName", "MapManBin", "clusterID", "clusterName",
    "x", "y", "clusterSimplifiedNodeDegree", "expressed",
]
DINAR_EDGE_COLUMNS = [
    "geneID1", "geneID2", "reactionType", "clusterID_geneID1", "clusterID_geneID2",
    "clusterSimplifiedNodeDegree_geneID1", "clusterSimplifiedNodeDegree_geneID2", "exists",
]


def _text(x):
    if x is None or x == "":
        return _MISSING
    return str(x)


def _first(x):
    '''The first value of a list attribute (e.g. a gene network gene's first cluster).'''
    if isinstance(x, (list, tuple)):
        return next((v for v in x if v is not None and v != ""), None)
    return x


def _cluster_labels(g, clusters):
    '''node -> cluster label (None: not in a cluster)'''
    if clusters is None:
        return {n: 1 for n in g}
    if isinstance(clusters, str):
        return {n: _first(d.get(clusters)) for n, d in g.nodes(data=True)}
    return {n: clusters.get(n) for n in g}


def _cluster_ids(labels):
    '''Cluster label -> DiNAR clusterID (1, 2, ...; largest cluster first, ties by label).'''
    sizes = {}
    for label in labels.values():
        if label is not None:
            sizes[label] = sizes.get(label, 0) + 1
    order = sorted(sizes, key=lambda label: (-sizes[label], str(label)))
    return {label: i for i, label in enumerate(order, start=1)}


def _cluster_layout(g, cluster_of, seed):
    '''Spring layout of each cluster on its own, clusters (and unclustered nodes) on a grid.'''
    groups = {}
    for n in g:
        groups.setdefault(cluster_of[n], []).append(n)
    order = sorted(groups, key=lambda c: (c == 0, c))  # unclustered (0) last

    columns = math.ceil(math.sqrt(len(order)))
    simple = nx.Graph(g)
    pos = {}
    for i, c in enumerate(order):
        cx, cy = 3 * (i % columns), -3 * (i // columns)
        layout = nx.spring_layout(simple.subgraph(groups[c]), seed=seed)  # within [-1, 1]
        pos.update({n: (cx + x, cy + y) for n, (x, y) in layout.items()})
    return pos


def _positions(g, positions, cluster_of, seed):
    if positions is None:
        if all("pos" in d for _, d in g.nodes(data=True)):
            positions = "pos"
        else:
            return _cluster_layout(g, cluster_of, seed)
    if isinstance(positions, str):
        positions = {n: d.get(positions) for n, d in g.nodes(data=True)}
    missing = [n for n in g if positions.get(n) is None]
    if missing:
        raise ValueError(f"No position for {len(missing)} nodes, e.g. {missing[:3]}.")
    return {n: (float(positions[n][0]), float(positions[n][1])) for n in g}


def to_dinar(g, clusters=None, positions=None, seed=0):
    '''Node and edge tables of `g` in DiNAR's format.

    Parameters
    ----------
    g : networkx.Graph
        PSS or CKN network (any networkx graph type), e.g. a PSS gene network.
    clusters : str or dict, optional
        The clusters DiNAR shows one at a time: a node attribute (e.g. ``"pathway"``; for
        list attributes, the first value is used), or a ``{node: cluster}`` mapping (e.g.
        from community detection). Nodes without a cluster get ``clusterID`` 0, which DiNAR
        doesn't show. Default: all nodes in one cluster.
    positions : str or dict, optional
        Node coordinates: a ``{node: (x, y)}`` mapping (as returned by the networkx layout
        functions), or a node attribute holding them. Default: the ``pos`` node attribute
        (the networkx convention, ``nx.set_node_attributes(g, pos, "pos")``) if every node
        has one; otherwise a spring layout of each cluster, with the clusters on a grid.
    seed : int, optional
        Random seed for the default layout (default 0, so the layout is reproducible).

    Returns
    -------
    nodes : pandas.DataFrame
        One row per node: ``geneID`` (the node id), ``shortDescription``
        (``description``), ``shortName`` (``display_label``, or the node id),
        ``MapManBin`` (the ``mapman`` bins joined with ``" | "``), ``clusterID`` (1, 2, ...,
        largest cluster first), ``clusterName`` (the cluster label), ``x``, ``y``,
        ``clusterSimplifiedNodeDegree`` (number of neighbours in the same cluster) and
        ``expressed`` (1).
    edges : pandas.DataFrame
        ``geneID1``, ``geneID2``, ``reactionType`` (the edge's ``interaction``, e.g.
        ``positive-influence``), the clusters and cluster degrees of both nodes, and
        ``exists`` (1). Parallel edges with the same ``interaction`` (e.g. one per PSS
        reaction) are one row.

    Missing values are ``-``.

    Raises
    ------
    ValueError
        If `positions` has no coordinates for some nodes.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.MultiDiGraph()
    >>> g.add_node("AT2G38470", display_label="WRKY33", mapman=["26.11.3.2.1_External stimuli"])
    >>> g.add_node("camalexin")
    >>> _ = g.add_edge("AT2G38470", "camalexin", interaction="positive-influence")
    >>> nodes, edges = to_dinar(g)
    >>> nodes.loc[0, ["geneID", "shortName", "MapManBin"]].tolist()
    ['AT2G38470', 'WRKY33', '26.11.3.2.1_External stimuli']
    '''
    labels = _cluster_labels(g, clusters)
    ids = _cluster_ids(labels)
    cluster_of = {n: ids.get(label, 0) if label is not None else 0 for n, label in labels.items()}
    pos = _positions(g, positions, cluster_of, seed)

    # degree within the cluster, in the simplified network (undirected, no parallel edges
    # or self-loops)
    simple = nx.Graph(g)
    simple.remove_edges_from(nx.selfloop_edges(simple))
    degree = {
        n: sum(1 for m in simple[n] if cluster_of[m] == cluster_of[n]) for n in simple
    }

    nodes = pd.DataFrame(
        [
            [
                str(n),
                _text(d.get("description")),
                _text(d.get("display_label") or n),
                _text(" | ".join(m for m in (d.get("mapman") or []) if m)),
                cluster_of[n],
                _text(labels[n]),
                pos[n][0],
                pos[n][1],
                degree[n],
                1,
            ]
            for n, d in g.nodes(data=True)
        ],
        columns=DINAR_NODE_COLUMNS,
    )

    edges = pd.DataFrame(
        [
            [str(u), str(v), _text(d.get("interaction")), cluster_of[u], cluster_of[v],
             degree[u], degree[v], 1]
            for u, v, d in g.edges(data=True)
        ],
        columns=DINAR_EDGE_COLUMNS,
    ).drop_duplicates(ignore_index=True)

    return nodes, edges


def write_dinar(g, nodes_path, edges_path, clusters=None, positions=None, seed=0):
    '''Write `g` as DiNAR node and edge tables (tab-separated), see :func:`to_dinar`.

    Parameters
    ----------
    g : networkx.Graph
        Network to write.
    nodes_path, edges_path : str or pathlib.Path
        Files to write (overwritten if they exist).
    clusters, positions, seed
        As for :func:`to_dinar`.
    '''
    nodes, edges = to_dinar(g, clusters=clusters, positions=positions, seed=seed)
    nodes.to_csv(nodes_path, sep="\t", index=False)
    edges.to_csv(edges_path, sep="\t", index=False)
