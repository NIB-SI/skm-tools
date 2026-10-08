'''Export networks as DiNAR node and edge tables.

`DiNAR <https://github.com/NIB-SI/DiNAR>`_ (Differential Network Analysis in R) shows
experimental data on a knowledge network, one cluster of nodes at a time. Its custom
networks are tab-separated node and edge tables with fixed column names (as DiNAR's
``Networks/PSS/Ath-based_*.txt``); missing values are ``-``.

DiNAR shows one cluster at a time (cluster 0 not at all), and can't show clusters with more
than 2000 edges.
'''

import csv
import math
import re

import networkx as nx
import pandas as pd


_MISSING = "-"
MAX_CLUSTER_EDGES = 2000

DINAR_NODE_COLUMNS = [
    "geneID", "shortDescription", "shortName", "MapManBin", "clusterID",
    "x", "y", "clusterSimplifiedNodeDegree", "expressed",
]
DINAR_EDGE_COLUMNS = [
    "geneID1", "geneID2", "reactionType", "clusterID_geneID1", "clusterID_geneID2",
    "clusterSimplifiedNodeDegree_geneID1", "clusterSimplifiedNodeDegree_geneID2", "exists",
]


def _missing(x):
    return x is None or x == "" or (isinstance(x, float) and math.isnan(x))


def _text(x):
    '''A table value: "-" if missing; tabs and line breaks (which would break DiNAR's
    tables) replaced by spaces.'''
    if _missing(x):
        return _MISSING
    return re.sub(r"[\t\r\n]+", " ", str(x))


def _first(x):
    '''The first value of a list attribute (e.g. a gene network gene's first cluster; sets
    sorted), None if missing.'''
    if isinstance(x, (set, frozenset)):
        x = sorted(x, key=str)
    if isinstance(x, (list, tuple)):
        return next((v for v in x if not _missing(v)), None)
    return None if _missing(x) else x


def _reaction_type(d):
    '''DiNAR draws edges by the start of reactionType: "act" (an arrow), "inh" (a bar),
    anything else dashed; as DiNAR's PSS tables, act_/inh_/unk_ followed by the PSS
    reaction type (e.g. "act_protein activation"), or "influence" (CKN).'''
    sign = {"positive-influence": "act", "negative-influence": "inh"}.get(d.get("interaction"), "unk")
    reaction_type = d.get("reaction_type")
    if isinstance(reaction_type, (list, tuple, set)):
        reaction_type = "/".join(sorted(map(str, reaction_type)))
    return f"{sign}_{_text(reaction_type) if not _missing(reaction_type) else 'influence'}"


def _cluster_labels(g, clusters):
    '''node -> cluster label (None: not in a cluster)'''
    if clusters is None:
        return {n: 1 for n in g}
    if isinstance(clusters, str):
        return {n: _first(d.get(clusters)) for n, d in g.nodes(data=True)}
    return {n: _first(clusters.get(n)) for n in g}


def _cluster_ids(labels):
    '''Cluster label -> DiNAR clusterID (1, 2, ...; largest cluster first, ties by label;
    nodes without a label (None) together, last).'''
    sizes = {}
    for label in labels.values():
        if label is not None:
            sizes[label] = sizes.get(label, 0) + 1
    order = sorted(sizes, key=lambda label: (-sizes[label], str(label)))
    ids = {label: i for i, label in enumerate(order, start=1)}
    ids[None] = len(order) + 1
    return ids


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
        from community detection). Nodes without a cluster (None or NaN) are put together
        in one more cluster. Default: all nodes in one cluster (fine for small networks,
        e.g. PSS). Nodes without edges get ``clusterID`` 0, which DiNAR doesn't show.
    positions : str or dict, optional
        Node coordinates: a ``{node: (x, y)}`` mapping (as returned by the networkx layout
        functions), or a node attribute holding them. Default: the ``pos`` node attribute
        (the networkx convention, ``nx.set_node_attributes(g, pos, "pos")``) if every node
        has one; otherwise a spring layout of each cluster, with the clusters on a grid.
    seed : int, optional
        Random seed for the default layout (default 0, so the layout is the same each
        time, with the same networkx and numpy versions).

    Returns
    -------
    nodes : pandas.DataFrame
        One row per node: ``geneID`` (the node id), ``shortDescription``
        (``description``), ``shortName`` (``display_label``, or the node id),
        ``MapManBin`` (the ``mapman`` bins joined with ``" | "``), ``clusterID`` (1, 2, ...,
        largest cluster first; 0 for nodes without edges), ``x``, ``y``,
        ``clusterSimplifiedNodeDegree`` (number of neighbours in the same cluster) and
        ``expressed`` (1).
    edges : pandas.DataFrame
        ``geneID1``, ``geneID2``, ``reactionType``, the clusters and cluster degrees of both
        nodes, and ``exists`` (1). ``reactionType`` is ``act_``, ``inh_`` or ``unk_`` (from
        the ``interaction``: positive, negative or unknown influence; DiNAR draws edges by
        it), followed by the PSS ``reaction_type`` (e.g. ``act_protein activation``), or
        ``influence`` for networks without one (CKN). Parallel edges with the same
        ``reactionType`` are one row.

    Missing values are ``-``.

    Raises
    ------
    ValueError
        If `positions` has no coordinates for some nodes, two node ids are the same as
        text (e.g. ``1`` and ``"1"``), or a cluster has more than 2000 edges (DiNAR can't
        show it: use smaller clusters).

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
    if len({str(n) for n in g}) < g.number_of_nodes():
        raise ValueError("Node ids that are the same as text (e.g. 1 and '1').")

    # the simplified network: undirected, no parallel edges or self-loops
    simple = nx.Graph(g)
    simple.remove_edges_from(nx.selfloop_edges(simple))

    labels = _cluster_labels(g, clusters)
    # nodes without edges: cluster 0 (not shown); the others numbered by cluster size
    ids = _cluster_ids({n: label for n, label in labels.items() if simple.degree(n)})
    cluster_of = {n: ids[label] if simple.degree(n) else 0 for n, label in labels.items()}
    pos = _positions(g, positions, cluster_of, seed)

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
            [str(u), str(v), _reaction_type(d), cluster_of[u], cluster_of[v],
             degree[u], degree[v], 1]
            for u, v, d in g.edges(data=True)
        ],
        columns=DINAR_EDGE_COLUMNS,
    ).drop_duplicates(ignore_index=True)

    within = edges[edges["clusterID_geneID1"] == edges["clusterID_geneID2"]]
    too_large = within["clusterID_geneID1"].value_counts()
    too_large = too_large[too_large > MAX_CLUSTER_EDGES]
    if len(too_large):
        raise ValueError(f"DiNAR can't show clusters with more than {MAX_CLUSTER_EDGES} edges; "
                         f"cluster(s) {dict(too_large)} (clusterID: edges). Use smaller clusters "
                         "(`clusters`).")

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
    # DiNAR reads the tables without quoting (values have no tabs or line breaks, see to_dinar)
    nodes.to_csv(nodes_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)
    edges.to_csv(edges_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)
