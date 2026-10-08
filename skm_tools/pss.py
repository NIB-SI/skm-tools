'''PSS (Plant Stress Signalling): load the network exports, filter, simplify and rewire.

The PSS exports (reaction graph, interaction network, per-species gene networks) are made by
skm-pss-export (https://github.com/NIB-SI/skm-pss-export).

The loaders record which network a graph is in ``g.graph["pss_network"]``
(``"reaction_graph"``, ``"interaction_network"`` or ``"gene_network"``, with the species of a
gene network in ``g.graph["species"]``); the functions here use it to check that they are
applied to a network they work on.
'''

import logging
from collections import Counter
from pathlib import Path

import networkx as nx
import pandas as pd

from .skm_download_urls import (
    PSS_GENE_NETWORK_EDGE_FILE,
    PSS_GENE_NETWORK_EDGE_URL,
    PSS_GENE_NETWORK_NODE_FILE,
    PSS_GENE_NETWORK_NODE_URL,
    PSS_INTERACTION_NETWORK_EDGE_FILE,
    PSS_INTERACTION_NETWORK_EDGE_URL,
    PSS_INTERACTION_NETWORK_NODE_FILE,
    PSS_INTERACTION_NETWORK_NODE_URL,
    PSS_REACTION_GRAPH_EDGE_FILE,
    PSS_REACTION_GRAPH_EDGE_URL,
    PSS_REACTION_GRAPH_NODE_FILE,
    PSS_REACTION_GRAPH_NODE_URL,
)
from .utils import (
    as_list,
    download_if_missing,
    merge_values,
    read_skm_table,
    remove_isolate_nodes,
    resolve_nodes,
)

logger = logging.getLogger(__name__)

REACTION_GRAPH = "reaction_graph"
INTERACTION_NETWORK = "interaction_network"
GENE_NETWORK = "gene_network"


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

# PSS exports (pss-export): see skm_tools.utils.read_skm_table (complex names, which
# contain "|", are never in a list)
_PSS_LIST_COLUMNS = {
    "synonyms",
    "all_pathways",
    "mapman",
    "external_links",
    "components",
    "component_cluster_ids",
    "genes",
}
# Gene network only: a gene in several functional clusters has one entry per cluster
# (in the same order) in each of these. Loaded as lists for every node, so the type
# doesn't depend on the row. (display_label stays a string: it's what to show.)
_PSS_GENE_CLUSTER_COLUMNS = {
    "short_name",
    "pathway",
    "functional_cluster_id",
}
_PSS_INT_COLUMNS = {
    "rank",
}
_PSS_BOOL_COLUMNS = {
    "directed",
    "location_putative",
    "source_location_putative",
    "target_location_putative",
}


def _read_pss_table(path, list_columns=()):
    '''Read a PSS export table: lists, booleans and integers converted, None for empty.'''
    df = pd.read_csv(path, sep="\t", nrows=0)
    list_columns = (_PSS_LIST_COLUMNS | set(list_columns)
                    | {c for c in df.columns if c.endswith("_homologues")})
    return read_skm_table(path, list_columns, _PSS_BOOL_COLUMNS, _PSS_INT_COLUMNS)


def _pss_export_to_networkx(edge_path, node_path, edge_url, node_url, edge_key, network,
                            node_list_columns=()):
    '''Build a MultiDiGraph from a PSS export: nodes (with attributes) from the node file,
    edges from the edge file, keyed by the `edge_key` column. Missing files are downloaded
    from the URLs.'''
    edge_df = _read_pss_table(download_if_missing(edge_path, edge_url))
    node_df = _read_pss_table(download_if_missing(node_path, node_url), node_list_columns)

    g = nx.MultiDiGraph(pss_network=network)
    g.add_nodes_from(
        (data.pop("id"), data) for data in node_df.to_dict("records")
    )
    g.add_edges_from(
        (data.pop("source"), data.pop("target"), data[edge_key], data)
        for data in edge_df.to_dict("records")
    )
    return g


def _paths(edge_path, node_path, data_dir, edge_file, node_file):
    data_dir = Path(data_dir)
    return (Path(edge_path) if edge_path is not None else data_dir / edge_file,
            Path(node_path) if node_path is not None else data_dir / node_file)


_PATH_PARAMETERS = '''
    edge_path : str or pathlib.Path, optional
        Path to the edge file (default ``<data_dir>/{edge_file}``); if the file does not
        exist, it is downloaded from skm.nib.si.

    node_path : str or pathlib.Path, optional
        Path to the node file (default ``<data_dir>/{node_file}``); downloaded as for
        `edge_path`.

    data_dir : str or pathlib.Path
        Folder for the default file names (default: the current folder).
'''


def pss_reaction_graph_to_networkx(edge_path=None, node_path=None, data_dir="."):
    edge_path, node_path = _paths(edge_path, node_path, data_dir,
                                  PSS_REACTION_GRAPH_EDGE_FILE, PSS_REACTION_GRAPH_NODE_FILE)
    return _pss_export_to_networkx(edge_path, node_path, PSS_REACTION_GRAPH_EDGE_URL,
                                   PSS_REACTION_GRAPH_NODE_URL, edge_key="role",
                                   network=REACTION_GRAPH)


pss_reaction_graph_to_networkx.__doc__ = ''' Load the PSS reaction graph export to a networkx
    directed multigraph, including node attributes.

    Entities and reactions are both nodes (reactions have `node_type` == "reaction"), with one edge
    per reaction participant: participant -> reaction for inputs and modifiers,
    reaction -> participant for products. This is the lossless form of PSS.
    Edges are keyed by the participant's `role`, since an entity can take part in the same
    reaction twice (e.g. as template and stimulator).

    Parameters
    ----------
''' + _PATH_PARAMETERS.format(edge_file=PSS_REACTION_GRAPH_EDGE_FILE,
                              node_file=PSS_REACTION_GRAPH_NODE_FILE)


def pss_interaction_network_to_networkx(edge_path=None, node_path=None, data_dir="."):
    edge_path, node_path = _paths(edge_path, node_path, data_dir,
                                  PSS_INTERACTION_NETWORK_EDGE_FILE, PSS_INTERACTION_NETWORK_NODE_FILE)
    return _pss_export_to_networkx(edge_path, node_path, PSS_INTERACTION_NETWORK_EDGE_URL,
                                   PSS_INTERACTION_NETWORK_NODE_URL, edge_key="reaction_id",
                                   network=INTERACTION_NETWORK)


pss_interaction_network_to_networkx.__doc__ = ''' Load the PSS interaction network export to a
    networkx directed multigraph, including node attributes.

    Entities are nodes, and edges are entity -> entity influences through reactions
    (`interaction`: positive-influence, negative-influence or unknown-influence).
    Edges are keyed by `reaction_id`, as several reactions can link the same node pair.
    Mutual influences (`directed` == False, e.g. binding partners) are already listed in
    both directions.

    Parameters
    ----------
''' + _PATH_PARAMETERS.format(edge_file=PSS_INTERACTION_NETWORK_EDGE_FILE,
                              node_file=PSS_INTERACTION_NETWORK_NODE_FILE)


def pss_gene_network_to_networkx(edge_path=None, node_path=None, species="ath", data_dir="."):
    edge_path, node_path = _paths(edge_path, node_path, data_dir,
                                  PSS_GENE_NETWORK_EDGE_FILE.format(species),
                                  PSS_GENE_NETWORK_NODE_FILE.format(species))
    g = _pss_export_to_networkx(edge_path, node_path,
                                PSS_GENE_NETWORK_EDGE_URL.format(species),
                                PSS_GENE_NETWORK_NODE_URL.format(species),
                                edge_key="reaction_id", network=GENE_NETWORK,
                                node_list_columns=_PSS_GENE_CLUSTER_COLUMNS)
    g.graph["species"] = species
    return g


pss_gene_network_to_networkx.__doc__ = ''' Load a PSS gene network export (one species) to a
    networkx directed multigraph, including node attributes.

    As the interaction network, but with functional clusters expanded into their genes
    of one species. Genes have the `node_type` of their cluster (``PlantCoding`` or
    ``PlantNonCoding``), and are the nodes with a `species`.
    Edges are keyed by `reaction_id`.

    A gene can be in several functional clusters, so `short_name`, `pathway` and
    `functional_cluster_id` are lists for every node (one entry per cluster, in the same order;
    a single entry for nodes that aren't genes); `display_label` stays a single string for
    display. A cluster's name (as in the interaction network, and in the edges'
    `source_entity` / `target_entity`) is ``short_name[functional_cluster_id]``. Complexes'
    `components` are interaction-network names; match genes to them through
    `component_cluster_ids` and `functional_cluster_id`.

    Parameters
    ----------
''' + _PATH_PARAMETERS.format(edge_file=PSS_GENE_NETWORK_EDGE_FILE.format("<species>"),
                              node_file=PSS_GENE_NETWORK_NODE_FILE.format("<species>")) + '''
    species : str
        Species code of the gene network (default ``"ath"``, see the SKM translations): for
        the default file names and the download. Stored in ``g.graph["species"]``.
'''


# ---------------------------------------------------------------------------
# Filtering, simplifying and rewiring
# ---------------------------------------------------------------------------

def _network_kind(g):
    '''Which PSS network `g` is: as recorded by the loaders in ``g.graph["pss_network"]``, or
    else from its nodes (reaction nodes: the reaction graph; genes with a species: a gene
    network).'''
    kind = g.graph.get("pss_network")
    if kind is not None:
        return kind
    kind = INTERACTION_NETWORK
    for _, d in g.nodes(data=True):
        if d.get("node_type") == "reaction":
            return REACTION_GRAPH
        if d.get("species"):
            kind = GENE_NETWORK
    return kind


def _check_not_reaction_graph(g, func):
    '''The filtering and simplifying functions work on entity -> entity influences, so not on
    the reaction graph (reactions as nodes, no ``interaction`` on the edges).'''
    if _network_kind(g) == REACTION_GRAPH:
        raise ValueError(f"{func} works on the PSS interaction network or a gene network, "
                         "not the reaction graph (it has reaction nodes).")


def _reaction_ids(edge_data):
    '''The reaction ids of an edge, as a list (merged edges, e.g. from simplify_pss, have
    several).'''
    r = edge_data.get("reaction_id")
    if r is None:
        return []
    return [r] if isinstance(r, str) else list(r)


def _node_reactions(g, n, reaction_graph):
    '''The reactions node `n` takes part in.'''
    if reaction_graph:
        return {m for m in nx.all_neighbors(g, n) if g.nodes[m].get("node_type") == "reaction"}
    edges = list(g.in_edges(n, data=True)) + list(g.out_edges(n, data=True))
    return {r for _, _, d in edges for r in _reaction_ids(d)}


def _merge_edge_data(edges, differing=None, log_level=logging.DEBUG):
    '''Merge the attribute dicts of edges into one: all their reaction ids (a sorted list), and
    for every other attribute the value of :func:`skm_tools.utils.merge_values`, applied in
    reaction id order. Counts the attributes whose values differed in `differing`, and logs
    them at `log_level`.'''
    edges = sorted(edges, key=lambda d: _reaction_ids(d))
    reaction_ids = sorted({r for d in edges for r in _reaction_ids(d)})
    data = {}
    for k in sorted(set().union(*(d.keys() for d in edges)) - {"reaction_id"}):
        value, differ = merge_values(k, [d.get(k) for d in edges])
        data[k] = value
        if differ:
            if differing is not None:
                differing[k] += 1
            logger.log(log_level, "%s --> %s: values %s, keeping %s", reaction_ids, k,
                         [d.get(k) for d in edges], value)
    data["reaction_id"] = reaction_ids
    return data


def remove_reactions(g, reaction_ids, remove_isolates=True):
    '''Remove PSS reactions, in place.

    Works on all three PSS networks: in the reaction graph, the reaction nodes are removed
    (with their participant edges); in the interaction network and the gene networks, the
    edges with these ``reaction_id`` values. Merged edges (e.g. from :func:`simplify_pss`)
    that also come from other reactions are kept, with only those reaction ids.

    Parameters
    ----------
    g : networkx.MultiDiGraph or networkx.DiGraph
        PSS network. Changed in place.
    reaction_ids : str or iterable of str
        Reactions to remove (e.g. ``["rx00001"]``); reactions not in `g` are ignored.
    remove_isolates : bool
        Also remove nodes without edges afterwards (default True); this includes nodes that
        had no edges before.

    Returns
    -------
    dict
        Removed node -> reason (``"reaction removed"`` for reaction nodes of the reaction
        graph, ``"isolate"``).

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.MultiDiGraph()
    >>> _ = g.add_edge("A", "B", key="rx1", reaction_id="rx1")
    >>> _ = g.add_edge("B", "C", key="rx2", reaction_id="rx2")
    >>> remove_reactions(g, "rx1")
    {'A': 'isolate'}
    '''
    reaction_ids = set(as_list(reaction_ids))
    reasons = {}

    if _network_kind(g) == REACTION_GRAPH:
        reactions = [r for r in reaction_ids
                     if r in g and g.nodes[r].get("node_type") == "reaction"]
        g.remove_nodes_from(reactions)
        reasons.update({r: "reaction removed" for r in reactions})
    elif reaction_ids:
        edges = g.edges(keys=True, data=True) if g.is_multigraph() else g.edges(data=True)
        to_remove = []
        for *e, d in edges:
            ids = _reaction_ids(d)
            if not reaction_ids.intersection(ids):
                continue
            kept = [r for r in ids if r not in reaction_ids]
            if kept:
                d["reaction_id"] = kept
            else:
                to_remove.append(tuple(e))
        g.remove_edges_from(to_remove)

    if remove_isolates:
        reasons.update(remove_isolate_nodes(g))

    return reasons


def remove_deadend_complexes(g):
    '''Remove complexes without outgoing edges, in place.

    A complex that influences nothing (is not an input or modifier of any reaction) is a
    dead end in directed analyses. Repeated until no dead-end complexes are left, since
    removing one can leave another without outgoing edges.

    In the interaction network and the gene networks, only the complexes are removed: the
    binding partners' mutual edges from the same reactions are kept. In the reaction graph,
    the reactions producing the complexes are removed (see :func:`remove_reactions`), and
    then the nodes left without edges.

    Parameters
    ----------
    g : networkx.DiGraph or networkx.MultiDiGraph
        PSS network. Changed in place.

    Returns
    -------
    dict
        Removed node -> reason (``"dead-end complex"``; in the reaction graph also
        ``"reaction removed"`` and ``"isolate"``).
    '''
    reaction_graph = _network_kind(g) == REACTION_GRAPH
    reasons = {}

    while True:
        deadends = [n for n, data in g.nodes(data=True)
                    if data.get("node_type") == "Complex" and g.out_degree(n) == 0]
        if not deadends:
            break
        if reaction_graph:
            reactions = set().union(*(_node_reactions(g, c, True) for c in deadends))
            reasons.update(remove_reactions(g, reactions, remove_isolates=True))
        g.remove_nodes_from(deadends)
        reasons.update({n: "dead-end complex" for n in deadends})

    logger.info("Removed %d dead-end complexes.",
                sum(r == "dead-end complex" for r in reasons.values()))
    return reasons


def filter_pss_nodes(g, node_types=None, species=None, remove_isolates=True):
    '''Remove PSS nodes, in place.

    With `species`, the functional clusters without genes in the species, and the complexes
    they are components of, are removed with their reactions (see :func:`remove_reactions`),
    as in the species' gene network: applied to the interaction network, the result has the
    same reactions and edges as the gene network, with clusters instead of genes. With
    `node_types`, complexes with a removed component are removed too (using the
    ``components`` and ``component_cluster_ids`` node attributes; components not in `g` are
    ignored).

    Parameters
    ----------
    g : networkx.Graph
        PSS network. Changed in place.
    node_types : str or list of str, optional
        Keep only nodes of these ``node_type`` values (e.g. ``"PlantCoding"``, ``"Complex"``;
        in a gene network, genes have their cluster's class, e.g. ``"PlantCoding"``). Not for
        the reaction graph, where it would leave reactions with missing participants.
    species : str or list of str, optional
        Interaction network or reaction graph: remove the functional clusters without genes
        in any of these species (``<species>_homologues`` attributes, e.g. ``["stu"]``), the
        complexes with such a cluster among their components (``component_cluster_ids``),
        and their reactions. Abstract clusters (``PlantAbstract``) have no genes and are
        kept, as metabolites. Merged edges (e.g. from :func:`simplify_pss`) are kept if
        they also come from other reactions. A gene network is already for one species.
    remove_isolates : bool
        Also remove nodes without edges afterwards (default True); this includes nodes that
        had no edges before.

    Returns
    -------
    dict
        Removed node -> reason (``"species missing"``, ``"wrong node type"``,
        ``"complex component removed"``, ``"reaction removed"`` or ``"isolate"``).

    Raises
    ------
    ValueError
        If `node_types` is given for the reaction graph, or `species` for a gene network.
    '''
    node_types = as_list(node_types)
    species = as_list(species)
    kind = _network_kind(g)
    if node_types and kind == REACTION_GRAPH:
        raise ValueError("filter_pss_nodes(node_types=...) works on the PSS interaction network "
                         "or a gene network, not the reaction graph (it has reaction nodes).")
    if species and kind == GENE_NETWORK:
        raise ValueError("species filtering is for the interaction network; a gene network "
                         "is already for one species (load the gene network of the species "
                         "instead).")
    og_size = g.number_of_nodes()

    to_remove = set()
    reasons = {}

    if species:
        # gene clusters (functional clusters with genes: PlantCoding, PlantNonCoding; not
        # PlantAbstract, which has no genes and doesn't decide, as metabolites) without
        # genes in any of the species
        homologue_properties = [f"{sp}_homologues" for sp in species]
        missing_clusters = {
            data["functional_cluster_id"]: n for n, data in g.nodes(data=True)
            if data.get("functional_cluster_id") and data.get("node_type") != "PlantAbstract"
            and not any(data.get(h) for h in homologue_properties)
        }
        # complexes with such a cluster among their components (components not in `g` can't
        # be checked, and are kept)
        missing_complexes = [
            n for n, data in g.nodes(data=True)
            if set(data.get("component_cluster_ids") or []).intersection(missing_clusters)
        ]
        reasons.update({n: "species missing" for n in missing_clusters.values()})
        reasons.update({n: "complex component removed" for n in missing_complexes})

        # their reactions are not in the species: as in the species' gene network, remove
        # the whole reactions (not just these nodes' edges), then the nodes
        no_species = list(missing_clusters.values()) + missing_complexes
        reaction_graph = kind == REACTION_GRAPH
        removed_reactions = set().union(*(_node_reactions(g, n, reaction_graph) for n in no_species))
        reasons.update(remove_reactions(g, removed_reactions, remove_isolates=False))
        g.remove_nodes_from(no_species)

    if node_types:
        wrong_type = [n for n, data in g.nodes(data=True) if data.get("node_type") not in node_types]
        to_remove.update(wrong_type)
        reasons.update({n: "wrong node type" for n in wrong_type if n not in reasons})

    # now remove complexes with a component that is in the network, but would no longer be.
    # Components are matched by node id (`components`) and by functional cluster id
    # (`component_cluster_ids`): in a gene network, clusters are expanded into genes, so a
    # cluster is gone only once all its genes are.
    def present(nodes):
        ids = set()
        for n in nodes:
            ids.add(n)
            fc = g.nodes[n].get("functional_cluster_id")
            ids.update([fc] if isinstance(fc, str) else [x for x in (fc or []) if x])
        return ids

    before = present(g.nodes())
    gone = before - present(n for n in g.nodes() if n not in to_remove)
    complex_component_missing = [
        n for n, data in g.nodes(data=True)
        if gone.intersection((data.get("components") or []) + (data.get("component_cluster_ids") or []))
    ]
    to_remove.update(complex_component_missing)
    reasons.update({n: "complex component removed" for n in complex_component_missing if n not in reasons})

    g.remove_nodes_from(to_remove)

    # all nodes left without edges (also those without edges before filtering)
    if remove_isolates:
        reasons.update({n: r for n, r in remove_isolate_nodes(g).items() if n not in reasons})

    logger.info("Removed %d nodes from the network.", og_size - g.number_of_nodes())
    return reasons


def simplify_pss(g, split_on_attrs=None, verbose=False):
    '''Merge parallel edges (one per reaction) into one edge per node pair.

    Returns a new graph; `g` is unchanged. For the interaction network or a gene network,
    not the reaction graph (raises ValueError).

    Merged edges get all ``reaction_id`` values, as a sorted list (edges that aren't merged get
    a list too, so the attribute has one type). For every other attribute one value is kept,
    by this rule (:func:`skm_tools.utils.merge_values`), ignoring missing values (None):

    - all values equal: that value;
    - ``interaction``: ``"unknown-influence"`` if they differ (e.g. a positive and a negative
      influence);
    - ``directed``: True if any edge is directed;
    - ``rank``: the lowest;
    - anything else: the value of the edge with the first reaction id.

    Parameters
    ----------
    g : networkx.MultiDiGraph
        PSS interaction network or gene network.
    split_on_attrs : str or list of str, optional
        Edge attributes (e.g. ``["interaction"]``) that must not be merged away: parallel
        edges are only merged with others that have the same values for all of them.
    verbose : bool
        Log every merged edge whose attributes differed, with the value kept (INFO; otherwise
        DEBUG). A summary, with the number of merged edges per differing attribute, is always
        logged (INFO).

    Returns
    -------
    networkx.DiGraph or networkx.MultiDiGraph
        A DiGraph, or a MultiDiGraph with `split_on_attrs` (as edges that differ on them
        stay separate). Graph attributes (e.g. ``pss_network``) are copied.

    Raises
    ------
    ValueError
        If `g` is the reaction graph, or not a multigraph (e.g. already simplified).
    '''
    _check_not_reaction_graph(g, "simplify_pss")
    if not g.is_multigraph():
        raise ValueError("simplify_pss merges the parallel edges of a multigraph; this graph "
                         "has none (already simplified?).")
    split_on_attrs = as_list(split_on_attrs) or []

    new_g = nx.MultiDiGraph() if split_on_attrs else nx.DiGraph()
    new_g.graph.update(g.graph)
    new_g.add_nodes_from(g.nodes(data=True))
    differing = Counter()  # attribute -> number of merged edges with differing values

    log_level = logging.INFO if verbose else logging.DEBUG

    for source, target in dict.fromkeys(g.edges()):
        groups = {}
        for d in g[source][target].values():
            groups.setdefault(tuple(d.get(a) for a in split_on_attrs), []).append(d)
        for group_edges in groups.values():
            if len(group_edges) == 1:
                data = {**group_edges[0], "reaction_id": _reaction_ids(group_edges[0])}
            else:
                data = _merge_edge_data(group_edges, differing, log_level)
            new_g.add_edge(source, target, **data)

    if differing:
        logger.info("Merged edges with differing values: %s. See the docstring for the value "
                    "kept, and verbose=True for details.",
                    ", ".join(f"{k} ({n})" for k, n in differing.most_common()))

    return new_g


# signs of a chain A -> B -> C
_COMPOSED_INTERACTION = {
    ("positive-influence", "positive-influence"): "positive-influence",
    ("positive-influence", "negative-influence"): "negative-influence",
    ("negative-influence", "positive-influence"): "negative-influence",
    ("negative-influence", "negative-influence"): "positive-influence",
}


def _propagates(e):
    '''Whether rewiring may go through edge `e`: not through the mutual edges between binding
    partners (partner <-> partner), only the complex-forming ones (partner -> complex).'''
    return not (e.get("reaction_type") == "binding/oligomerisation"
                and e.get("target_role") != "product")


def _rewired_edge(up, down, node):
    '''The edge replacing up (A -> node) and down (node -> C).'''
    data = {}
    for k in set(up) | set(down):
        if k.startswith("source_"):
            data[k] = up.get(k)
        elif k.startswith("target_"):
            data[k] = down.get(k)
        else:
            data[k] = merge_values(k, [up.get(k), down.get(k)])[0]
    data["reaction_id"] = sorted(set(_reaction_ids(up) + _reaction_ids(down)))
    data["interaction"] = _COMPOSED_INTERACTION.get(
        (up.get("interaction"), down.get("interaction")), "unknown-influence")
    data["directed"] = True
    data["note"] = f"rewired through {node}"
    return data


def remove_and_rewire(g, nodes, dry_run=False):
    '''Remove nodes, connecting each of their upstream nodes to each of their downstream nodes.

    Changes `g` in place. The nodes are removed one at a time, in the given order, each
    against the graph as rewired so far, so chains of removed nodes (A -> B -> C -> D,
    removing B and C) are bridged (A -> D).

    Rewiring doesn't go through the mutual edges between binding partners
    (partner <-> partner), in either direction, only through the complex-forming ones
    (partner -> complex): A -> B and B <-> D (binding) give no A -> D, but A -> B and
    B -> B|D give A -> B|D.

    A new edge A -> C, replacing A -> B -> C:

    - ``reaction_id``: the reaction ids of both edges (a sorted list);
    - ``interaction``: the sign of the chain (positive and negative: negative; two negatives:
      positive; anything with an unknown influence: unknown);
    - ``directed``: True;
    - ``source_*`` attributes from A -> B, ``target_*`` attributes from B -> C, and other
      attributes as in :func:`simplify_pss` (:func:`skm_tools.utils.merge_values`);
    - ``note``: ``"rewired through B"``.

    If A -> C already exists, the new edge is merged into it, as in :func:`simplify_pss`.

    Parameters
    ----------
    g : networkx.DiGraph
        A simple directed graph, e.g. from :func:`simplify_pss`. Multigraphs are not supported.
    nodes : node or iterable of nodes
        Nodes to remove (nodes not in `g` are ignored, with a warning).
    dry_run : bool
        Don't change `g`; return the edges that would be added instead.

    Returns
    -------
    dict or list
        Removed node -> reason (``"rewired"``, or ``"nothing to rewire"`` if no edge replaced
        it); with `dry_run`, the edges ``(u, v, data)`` that would be added or changed.

    Raises
    ------
    NotImplementedError
        If `g` is a multigraph or undirected.
    ValueError
        If `g` is the reaction graph.
    '''
    _check_not_reaction_graph(g, "remove_and_rewire")
    if g.is_multigraph() or not g.is_directed():
        raise NotImplementedError("Currently only implemented for DiGraph, "
                                  "see simplify_pss.")
    nodes, _ = resolve_nodes(g, nodes)

    h = g.copy() if dry_run else g
    reasons = {}
    changed = set()  # (u, v) edges added or changed

    for node in dict.fromkeys(nodes):
        upstream = [u for u in h.predecessors(node) if u != node and _propagates(h[u][node])]
        downstream = [v for v in h.successors(node) if v != node and _propagates(h[node][v])]
        new_edges = [(u, v, _rewired_edge(h[u][node], h[node][v], node))
                     for u in upstream for v in downstream if u != v]

        h.remove_node(node)
        for u, v, data in new_edges:
            if h.has_edge(u, v):
                data = _merge_edge_data([h[u][v], data])
                h[u][v].clear()
            h.add_edge(u, v, **data)
            changed.add((u, v))
        reasons[node] = "rewired" if new_edges else "nothing to rewire"

    changed = sorted((u, v) for u, v in changed if h.has_edge(u, v))
    logger.info("Removed %d nodes, added or changed %d edges.", len(reasons), len(changed))

    if dry_run:
        return [(u, v, h[u][v]) for u, v in changed]
    return reasons


def remove_duplicated_binding_edges(g):
    '''Keep one direction of each mutual binding edge, in place.

    Binding partners are linked in both directions (A -> B and B -> A, same reaction); this
    removes one of the two (which one is arbitrary), e.g. for a less cluttered drawing.
    Don't use it before directed path or neighbourhood analysis.

    Parameters
    ----------
    g : networkx.DiGraph or networkx.MultiDiGraph
        PSS interaction network or gene network (not the reaction graph; each parallel edge
        of a multigraph is considered separately). Changed in place.

    Returns
    -------
    list
        The removed edges, ``(u, v)`` or, for a multigraph, ``(u, v, key)``.
    '''
    _check_not_reaction_graph(g, "remove_duplicated_binding_edges")
    is_multi = g.is_multigraph()

    def edge_items(u, v):
        '''(key, data) pairs for edges u->v; key is None for a non-multigraph'''
        if v not in g[u]:
            return []
        if is_multi:
            return list(g[u][v].items())
        return [(None, g[u][v])]

    edges_to_remove = set()
    for node in g.nodes():
        for upstream_node in g.predecessors(node):
            for u_key, e in edge_items(upstream_node, node):
                if e.get("reaction_type") != "binding/oligomerisation":
                    continue
                # the same reaction in the other direction
                for d_key, e2 in edge_items(node, upstream_node):
                    if e2.get("reaction_id") != e.get("reaction_id"):
                        continue
                    forward = (upstream_node, node, u_key) if is_multi else (upstream_node, node)
                    reverse = (node, upstream_node, d_key) if is_multi else (node, upstream_node)
                    # keep the "first" one: don't remove both
                    if reverse not in edges_to_remove:
                        edges_to_remove.add(forward)

    edges_to_remove = sorted(edges_to_remove, key=str)
    g.remove_edges_from(edges_to_remove)
    logger.info("Removed %d edges from the network.", len(edges_to_remove))
    return edges_to_remove
