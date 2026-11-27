'''Augment PSS with CKN: the CKN interactions around the PSS genes, and what they add.

PSS is curated and small; CKN is large and mostly from high-throughput data. Augmenting PSS
with the CKN edges that touch its genes shows which PSS genes CKN connects that PSS doesn't
(directly, or through one CKN node), as candidates for curation or hypotheses.

- :func:`build_augmented_pss`: PSS plus the CKN edges touching it, each node and edge
  tagged with its ``origin`` (``"PSS"`` or ``"CKN"``);
- :func:`find_new_direct_edges`, :func:`find_new_intermediate_pairs`: the PSS node pairs
  newly connected by CKN;
- :func:`extract_pair_subnetwork`, :func:`extract_node_neighbourhood`: small subnetworks to
  inspect a pair or a node, e.g. in Cytoscape with :func:`show_augmented`.

The PSS network and CKN must use the same node ids: use the PSS gene network of Arabidopsis
(:func:`skm_tools.pss.pss_gene_network_to_networkx`, ``species="ath"``), whose nodes are
Arabidopsis gene ids, as in CKN.

Contributed by Vid Modic.
'''

import logging
import random
from itertools import combinations

import networkx as nx
import pandas as pd

from skm_tools.pss import _check_not_reaction_graph

logger = logging.getLogger(__name__)

PSS, CKN = "PSS", "CKN"

# colours of the origins in Cytoscape (show_augmented)
ORIGIN_COLORS = {PSS: "#2171B5", CKN: "#E6550D"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_edge_either_way(g, u, v):
    '''Whether `g` has an edge between u and v, in either direction.'''
    return g.has_edge(u, v) or g.has_edge(v, u)


def _edge_dicts(g, u, v):
    '''The attribute dicts of the edges u -> v (several in a multigraph).'''
    if not g.has_edge(u, v):
        return []
    return list(g[u][v].values()) if g.is_multigraph() else [g[u][v]]


def _has_ckn_edge(g, u, v):
    '''Whether `g` has a CKN edge u -> v.'''
    return any(d.get("origin") == CKN for d in _edge_dicts(g, u, v))


def _ckn_edge_ranks(g, u, v):
    '''The ranks of the CKN edges between u and v, in either direction.'''
    return [d["rank"] for a, b in ((u, v), (v, u)) for d in _edge_dicts(g, a, b)
            if d.get("origin") == CKN and d.get("rank") is not None]


def _undirected_neighbours(g, n):
    '''Successors and predecessors of `n` (its neighbours, if `g` is undirected).'''
    if g.is_directed():
        return set(g.successors(n)) | set(g.predecessors(n))
    return set(g.neighbors(n))


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------

def build_augmented_pss(pss, ckn, between_added=True):
    '''PSS plus the CKN nodes and edges that touch it.

    A CKN edge u -> v is added if

    - u or v is a PSS node, unless both are and PSS has an edge between them (in either
      direction): PSS edges are never replaced or duplicated by CKN edges;
    - with `between_added`, also if u and v are both CKN nodes added by the rule above (so
      the CKN connections among the new nodes are kept too).

    Every node and edge gets an ``origin`` attribute, ``"PSS"`` or ``"CKN"``, and a ``step``:
    ``"input"`` (from PSS), ``"step1"`` (CKN, touching a PSS node) or ``"step1+step2"``
    (CKN, between two CKN nodes added in step 1). Other attributes are copied as they are;
    PSS nodes keep their PSS attributes.

    Parameters
    ----------
    pss : networkx.DiGraph or networkx.MultiDiGraph
        PSS network, with the node ids of CKN: the Arabidopsis gene network (or a part of it).
        Not changed.
    ckn : networkx.DiGraph or networkx.MultiDiGraph
        CKN (e.g. :func:`skm_tools.ckn.ckn_to_networkx`, filtered as needed). Not changed.
    between_added : bool
        Also add the CKN edges between added CKN nodes (default True).

    Returns
    -------
    networkx.MultiDiGraph
        The augmented network, with the graph attributes of `pss`.

    Raises
    ------
    ValueError
        If `pss` is the PSS reaction graph.
    '''
    _check_not_reaction_graph(pss, "build_augmented_pss")
    augmented = nx.MultiDiGraph()
    augmented.graph.update(pss.graph)

    for n, d in pss.nodes(data=True):
        augmented.add_node(n, **{**d, "origin": PSS, "step": "input"})
    if pss.is_multigraph():
        for u, v, k, d in pss.edges(keys=True, data=True):
            augmented.add_edge(u, v, key=k, **{**d, "origin": PSS, "step": "input"})
    else:
        for u, v, d in pss.edges(data=True):
            augmented.add_edge(u, v, **{**d, "origin": PSS, "step": "input"})

    pss_nodes = set(pss)
    if not any(n in ckn for n in pss_nodes):
        logger.warning("No PSS node is in CKN: the networks need the same node ids, e.g. the "
                       "Arabidopsis gene network (pss_gene_network_to_networkx(species='ath')).")

    added_nodes = set()
    n_touching = 0
    for u, v, d in ckn.edges(data=True):
        u_in, v_in = u in pss_nodes, v in pss_nodes
        if not (u_in or v_in) or (u_in and v_in and _has_edge_either_way(pss, u, v)):
            continue
        for n in (u, v):
            if n not in augmented:
                augmented.add_node(n, **{**ckn.nodes[n], "origin": CKN, "step": "step1"})
                added_nodes.add(n)
        augmented.add_edge(u, v, **{**d, "origin": CKN, "step": "step1"})
        n_touching += 1

    n_between = 0
    if between_added:
        for u, v, d in ckn.edges(data=True):
            if u in added_nodes and v in added_nodes:
                augmented.add_edge(u, v, **{**d, "origin": CKN, "step": "step1+step2"})
                n_between += 1

    logger.info("Augmented PSS (%d nodes, %d edges) with %d CKN nodes and %d CKN edges (%d touching "
                "PSS, %d between added CKN nodes); %d PSS nodes got no CKN edge.",
                pss.number_of_nodes(), pss.number_of_edges(), len(added_nodes),
                n_touching + n_between, n_touching, n_between,
                len(check_unaugmented_nodes(pss, augmented)))
    return augmented


def check_unaugmented_nodes(pss, augmented):
    '''The PSS nodes without a CKN edge in the augmented network.

    A PSS node gets none if it isn't in CKN, or if all its CKN edges are to PSS nodes it
    already has a PSS edge with. Worth checking for small, hand-made subnetworks, e.g. for
    mistyped node ids.

    Parameters
    ----------
    pss : networkx.Graph
        The PSS network given to :func:`build_augmented_pss`.
    augmented : networkx.MultiDiGraph
        Its result.

    Returns
    -------
    set
    '''
    touched = {n for u, v, d in augmented.edges(data=True) if d.get("origin") == CKN for n in (u, v)}
    return set(pss) - touched


def get_augmentation_summary(pss, augmented):
    '''Node and edge counts, before and after augmenting.

    Returns
    -------
    dict
        ``pss_nodes``, ``pss_edges``, ``augmented_nodes``, ``augmented_edges``,
        ``ckn_nodes_added``, ``ckn_edges_added``.
    '''
    return {
        "pss_nodes": pss.number_of_nodes(),
        "pss_edges": pss.number_of_edges(),
        "augmented_nodes": augmented.number_of_nodes(),
        "augmented_edges": augmented.number_of_edges(),
        "ckn_nodes_added": sum(1 for _, d in augmented.nodes(data=True) if d.get("origin") == CKN),
        "ckn_edges_added": sum(1 for *_, d in augmented.edges(data=True) if d.get("origin") == CKN),
    }


def augmentation_diff(augmented):
    '''What CKN added, as a table: a row per CKN node and per CKN edge.

    Returns
    -------
    pandas.DataFrame
        Columns ``kind`` (``"node"`` or ``"edge"``), ``source`` (the node, for nodes),
        ``target`` (None for nodes), and the attributes.
    '''
    rows = [{"kind": "node", "source": n, "target": None, **d}
            for n, d in augmented.nodes(data=True) if d.get("origin") == CKN]
    rows += [{"kind": "edge", "source": u, "target": v, **d}
             for u, v, d in augmented.edges(data=True) if d.get("origin") == CKN]
    return pd.DataFrame(rows, columns=None if rows else ["kind", "source", "target"])


# ---------------------------------------------------------------------------
# What the augmentation adds
# ---------------------------------------------------------------------------

def find_new_direct_edges(augmented, pss, max_rank=None):
    '''PSS node pairs connected by a CKN edge, without a PSS edge between them (self-loops
    are not pairs).

    Parameters
    ----------
    augmented : networkx.MultiDiGraph
        From :func:`build_augmented_pss`.
    pss : networkx.Graph
        The PSS network it was built from.
    max_rank : int, optional
        Keep only the pairs whose best (lowest) CKN rank is at most this.

    Returns
    -------
    pandas.DataFrame
        A row per pair: ``u``, ``v``, ``rank`` (the best CKN rank between them, either
        direction), ``direction`` (``"u->v"``, ``"v->u"`` or ``"both"``); by rank.
    '''
    pss_nodes = set(pss)
    seen = set()
    rows = []
    for u, v, d in augmented.edges(data=True):
        if (d.get("origin") != CKN or u == v or u not in pss_nodes or v not in pss_nodes
                or _has_edge_either_way(pss, u, v)):
            continue
        pair = tuple(sorted((u, v), key=str))
        if pair in seen:
            continue
        seen.add(pair)

        ranks = _ckn_edge_ranks(augmented, u, v)
        best = min(ranks) if ranks else None
        if max_rank is not None and (best is None or best > max_rank):
            continue
        forward, backward = _has_ckn_edge(augmented, u, v), _has_ckn_edge(augmented, v, u)
        rows.append({"u": u, "v": v, "rank": best,
                     "direction": "both" if forward and backward else "u->v" if forward else "v->u"})

    df = pd.DataFrame(rows, columns=["u", "v", "rank", "direction"])
    return df.sort_values(["rank", "u", "v"]).reset_index(drop=True)


def find_new_intermediate_pairs(augmented, pss, max_rank=None):
    '''PSS node pairs connected through one CKN node, without a PSS edge between them.

    Parameters
    ----------
    augmented : networkx.MultiDiGraph
        From :func:`build_augmented_pss`.
    pss : networkx.Graph
        The PSS network it was built from.
    max_rank : int, optional
        Keep only the paths whose worse step has at most this (best) rank.

    Returns
    -------
    pandas.DataFrame
        A row per pair and intermediate: ``u``, ``v``, ``intermediate``, ``rank_ux``,
        ``rank_xv`` (the best CKN rank of each step, either direction) and ``path_rank``
        (the higher of the two); by path rank. Directions are not considered.
    '''
    pss_nodes = set(pss)
    rows = []
    for x, d in augmented.nodes(data=True):
        if d.get("origin") != CKN:
            continue
        partners = sorted((n for n in _undirected_neighbours(augmented, x) if n in pss_nodes), key=str)
        for a, b in combinations(partners, 2):
            if _has_edge_either_way(pss, a, b):
                continue
            ranks_ax, ranks_xb = _ckn_edge_ranks(augmented, a, x), _ckn_edge_ranks(augmented, x, b)
            if not ranks_ax or not ranks_xb:
                continue
            rank_ax, rank_xb = min(ranks_ax), min(ranks_xb)
            path_rank = max(rank_ax, rank_xb)
            if max_rank is not None and path_rank > max_rank:
                continue
            rows.append({"u": a, "v": b, "intermediate": x, "rank_ux": rank_ax, "rank_xv": rank_xb,
                         "path_rank": path_rank})

    df = pd.DataFrame(rows, columns=["u", "v", "intermediate", "rank_ux", "rank_xv", "path_rank"])
    return df.sort_values(["path_rank", "u", "v", "intermediate"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Inspecting pairs and nodes
# ---------------------------------------------------------------------------

def extract_pair_subnetwork(augmented, pss, a, b, intermediate=None, max_neighbours=None):
    '''A small subnetwork around a pair (a, b) of the augmented network, e.g. for Cytoscape.

    Nodes: `a` and `b` (node attribute ``role`` ``"focus"``), `intermediate`
    (``"intermediate"``), and the neighbours of `a` and `b` in PSS (``"context"``; at most
    `max_neighbours` each, the first by name). Edges: all edges of `augmented` between them,
    with ``edge_status`` ``"new_from_ckn"`` (CKN edges between the focus and intermediate
    nodes), ``"pss_original"`` or ``"ckn_context"``.

    Parameters
    ----------
    augmented : networkx.MultiDiGraph
    pss : networkx.Graph
    a, b : node
    intermediate : node, optional
        E.g. from :func:`find_new_intermediate_pairs`.
    max_neighbours : int, optional

    Returns
    -------
    networkx.MultiDiGraph
        A copy.
    '''
    focal = {a, b} | ({intermediate} if intermediate is not None else set())

    def context(node):
        if node not in pss:
            return set()
        neighbours = _undirected_neighbours(pss, node) - focal
        if max_neighbours is not None:
            neighbours = set(sorted(neighbours, key=str)[:max_neighbours])
        return neighbours

    nodes = (focal | context(a) | context(b)) & set(augmented)
    sub = augmented.subgraph(nodes).copy()
    for n, d in sub.nodes(data=True):
        d["role"] = "focus" if n in (a, b) else "intermediate" if n == intermediate else "context"
    for u, v, d in sub.edges(data=True):
        origin = d.get("origin", PSS)
        if origin == CKN and u in focal and v in focal:
            d["edge_status"] = "new_from_ckn"
        else:
            d["edge_status"] = "pss_original" if origin == PSS else "ckn_context"
    return sub


def pick_random_node(augmented, origin=None, seed=None):
    '''A random node of the augmented network, e.g. to inspect.

    Parameters
    ----------
    augmented : networkx.MultiDiGraph
    origin : {"PSS", "CKN"}, optional
        Pick only from these nodes (default: any).
    seed : int, optional
        For a reproducible pick.

    Raises
    ------
    ValueError
        If there is no such node.
    '''
    candidates = sorted((n for n, d in augmented.nodes(data=True)
                         if origin is None or d.get("origin") == origin), key=str)
    if not candidates:
        raise ValueError(f"No nodes of origin {origin!r}.")
    return random.Random(seed).choice(candidates)


def extract_node_neighbourhood(augmented, node, radius=1, max_neighbours=None, seed=0):
    '''The neighbourhood of a node in the augmented network (edges in either direction).

    Radius 1 shows whether CKN changed what the node is connected to. Radius 2 can show
    crosstalk through a CKN node, but CKN hubs connect hundreds of nodes: use it with
    `max_neighbours` (e.g. 15 to 25).

    Parameters
    ----------
    augmented : networkx.MultiDiGraph
    node : node
    radius : int
        Steps from `node` (default 1).
    max_neighbours : int, optional
        At most this many new neighbours of each node per step (a random sample).
    seed : int
        For the sample (default 0).

    Returns
    -------
    networkx.MultiDiGraph
        A copy, with node attribute ``role``: ``"focus"`` (`node`) or ``"context"``.

    Raises
    ------
    ValueError
        If `node` is not in `augmented`.
    '''
    if node not in augmented:
        raise ValueError(f"{node!r} is not in the augmented network.")
    undirected = augmented.to_undirected(as_view=True)
    if max_neighbours is None:
        nodes = set(nx.ego_graph(undirected, node, radius=radius))
    else:
        rng = random.Random(seed)
        nodes, frontier = {node}, {node}
        for _ in range(radius):
            new = set()
            for n in sorted(frontier, key=str):
                neighbours = sorted(set(undirected.neighbors(n)) - nodes - new, key=str)
                if len(neighbours) > max_neighbours:
                    neighbours = rng.sample(neighbours, max_neighbours)
                new.update(neighbours)
            nodes |= new
            frontier = new

    sub = augmented.subgraph(nodes).copy()
    for n, d in sub.nodes(data=True):
        d["role"] = "focus" if n == node else "context"
    return sub


def get_node_neighbourhood_augmentation(pss, augmented, node, radius=1, max_neighbours=None, seed=0):
    '''How much CKN added to the neighbourhood of a node.

    Parameters are as for :func:`extract_node_neighbourhood`.

    Returns
    -------
    dict
        ``node``, ``radius``, ``pss_neighbours`` (within `radius` in PSS),
        ``augmented_neighbours`` (in the augmented network), ``ckn_added_neighbours`` (of
        them, CKN nodes) and ``pct_new`` (their percentage).
    '''
    sub = extract_node_neighbourhood(augmented, node, radius, max_neighbours, seed)
    augmented_neighbours = sub.number_of_nodes() - 1
    pss_neighbours = (nx.ego_graph(pss.to_undirected(as_view=True), node, radius=radius).number_of_nodes() - 1
                      if node in pss else 0)
    ckn_added = sum(1 for _, d in sub.nodes(data=True) if d.get("origin") == CKN)
    return {
        "node": node,
        "radius": radius,
        "pss_neighbours": pss_neighbours,
        "augmented_neighbours": augmented_neighbours,
        "ckn_added_neighbours": ckn_added,
        "pct_new": 100 * ckn_added / augmented_neighbours if augmented_neighbours else 0.0,
    }


# ---------------------------------------------------------------------------
# Cytoscape
# ---------------------------------------------------------------------------

def show_augmented(g, title, collection=None, style="SKM-augmented"):
    '''Load (a part of) the augmented network into Cytoscape, coloured by origin.

    The SKM style, with edges and node borders coloured by ``origin`` (PSS blue, CKN
    orange, see :data:`ORIGIN_COLORS`); node borders wider for PSS nodes, or, for the
    subnetworks of :func:`extract_pair_subnetwork` and :func:`extract_node_neighbourhood`, by
    ``role`` (focus, intermediate, context); edges of ``edge_status`` ``"new_from_ckn"``
    wider. Use it on subnetworks: Cytoscape is slow with tens of thousands of nodes.

    Needs a running Cytoscape and the ``cytoscape`` extra.

    Parameters
    ----------
    g : networkx.Graph
        From :func:`build_augmented_pss`, or a subnetwork of it.
    title : str
        Network name.
    collection : str, optional
        Collection to add it to (default: a new one named `title`).
    style : str
        Name of the style made (a copy of the SKM style; replaced if it exists).

    Returns
    -------
    int
        SUID of the network.
    '''
    from skm_tools import cytoscape_utils as cu

    p4c = cu.p4c
    suid = cu.load_network(g, title, collection=collection)
    base = cu.apply_builtin_style("skm", network=suid)
    cu.copy_style(base, style, networks=[suid], overwrite=True)

    origins, colors = list(ORIGIN_COLORS), list(ORIGIN_COLORS.values())
    p4c.set_edge_color_mapping("origin", origins, colors, mapping_type="d", style_name=style, network=suid)
    p4c.set_node_border_color_mapping("origin", origins, colors, mapping_type="d", style_name=style,
                                      network=suid)
    if any(r is not None for _, r in g.nodes(data="role")):
        p4c.set_node_border_width_mapping("role", ["focus", "intermediate", "context"], [10, 6, 2],
                                          mapping_type="d", style_name=style, network=suid)
    else:
        p4c.set_node_border_width_mapping("origin", origins, [5, 2], mapping_type="d", style_name=style,
                                          network=suid)
    if any(s is not None for *_, s in g.edges(data="edge_status")):
        p4c.set_edge_line_width_mapping("edge_status", ["pss_original", "ckn_context", "new_from_ckn"],
                                        [2, 2, 6], mapping_type="d", style_name=style, network=suid)
    p4c.layout_network("force-directed", network=suid)
    p4c.fit_content(network=suid)
    return suid
