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
    "identifiers",
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
    "organ_putative",
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


# edge attributes with one value per merged edge (skm_tools.utils.merge_values); the others (but
# reaction_id) are per reaction, as lists aligned with reaction_id
_SUMMARY_ATTRS = ("interaction", "directed", "rank")


def _reaction_key(r):
    '''Sort key of reaction ids (None, for edges without one, first).'''
    return (r is not None, str(r))


def _per_reaction(edge_data):
    '''The reactions of an edge: reaction id -> its per-reaction attributes (all but
    ``reaction_id`` and the summary attributes). In a merged edge (a list ``reaction_id``) the
    attributes are lists aligned with ``reaction_id``; otherwise each reaction gets the values
    as they are. An edge without reaction ids is one reaction, None.'''
    ids = _reaction_ids(edge_data) or [None]
    aligned = isinstance(edge_data.get("reaction_id"), (list, tuple))
    keys = [k for k in edge_data if k != "reaction_id" and k not in _SUMMARY_ATTRS]
    reactions = {}
    for i, r in enumerate(ids):
        reactions.setdefault(r, {
            k: (v[i] if aligned and isinstance(v, (list, tuple)) and len(v) == len(ids) else v)
            for k, v in ((k, edge_data[k]) for k in keys)})
    return reactions


def _edge_data(reactions, summary):
    '''Edge attributes from per-reaction attributes (reaction id -> attributes) and the summary
    attributes: ``reaction_id`` a sorted list, the per-reaction attributes lists aligned with it.'''
    ids = sorted(reactions, key=_reaction_key)
    keys = sorted(set().union(*reactions.values())) if reactions else []
    data = dict(summary)
    data.update({k: [reactions[r].get(k) for r in ids] for k in keys})
    data["reaction_id"] = ids
    return data


def _edge_subset(edge_data, reaction_ids):
    '''The edge restricted to `reaction_ids` (a merged edge, or the edge itself if it has no
    others), or None if it has none of them. The summary attributes are kept as they are.'''
    reactions = _per_reaction(edge_data)
    keep = {r: a for r, a in reactions.items() if r in reaction_ids}
    if not keep:
        return None
    if len(keep) == len(reactions):
        return edge_data
    return _edge_data(keep, {k: edge_data[k] for k in _SUMMARY_ATTRS if k in edge_data})


def _merge_edge_data(edges, differing=None, log_level=logging.DEBUG, single=()):
    '''Merge the attribute dicts of edges into one: all their reaction ids (a sorted list), with
    the per-reaction attributes as lists aligned with them (a reaction in several of the edges:
    from the first), and the summary attributes (``interaction``, ``directed``, ``rank``) by
    :func:`skm_tools.utils.merge_values`. Counts the summary attributes whose values differed in
    `differing`, and logs them at `log_level`. The attributes in `single` (equal in all the
    edges) keep their one value.'''
    reactions = {}
    for d in edges:
        for r, attrs in _per_reaction(d).items():
            reactions.setdefault(r, attrs)
    summary = {}
    for k in _SUMMARY_ATTRS:
        values = [d.get(k) for d in edges if k in d]
        if not values:
            continue
        summary[k], differ = merge_values(k, values)
        if differ:
            if differing is not None:
                differing[k] += 1
            logger.log(log_level, "%s --> %s: values %s, keeping %s",
                       sorted(reactions, key=_reaction_key), k, values, summary[k])
    data = _edge_data(reactions, summary)
    data.update({k: edges[0][k] for k in single if k in edges[0]})
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
            kept = _edge_subset(d, set(_reaction_ids(d)) - reaction_ids)
            if kept is not None:
                d.clear()
                d.update(kept)
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
        If `node_types` is given for the reaction graph, `species` for a gene network, or a
        species without a ``<species>_homologues`` attribute in `g` (e.g. a misspelt code,
        or an export for some species only, which has ``genes`` instead).
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
    if species:
        attributes = set().union(*(data.keys() for _, data in g.nodes(data=True)))
        unknown = [sp for sp in species if f"{sp}_homologues" not in attributes]
        if unknown:
            hint = (" (the network has a 'genes' attribute instead: it was exported for one "
                    "species or several, not with all species)" if "genes" in attributes else "")
            raise ValueError(f"filter_pss_nodes(species=...): no <species>_homologues attribute "
                             f"for {', '.join(unknown)} in the network{hint}.")
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

    Each edge gets the ``reaction_id`` values of the edges merged into it, as a sorted list,
    and their other attributes per reaction: lists aligned with ``reaction_id`` (e.g.
    ``reaction_type``, ``source_role``, ``evidence_sentence``, ``influence_sbo``). Edges that
    aren't merged get lists too, so each attribute has one type. Three attributes are
    summarised in one value instead (:func:`skm_tools.utils.merge_values`), ignoring missing
    values (None):

    - all values equal: that value;
    - ``interaction``: ``"unknown-influence"`` if they differ (e.g. a positive and a negative
      influence; ``influence_sbo`` has the influence of each reaction);
    - ``directed``: True if any edge is directed;
    - ``rank``: the lowest.

    For example, the edge of reactions rx1 (catalysis) and rx2 (binding/oligomerisation):
    ``reaction_id`` ``["rx1", "rx2"]``, ``reaction_type``
    ``["catalysis", "binding/oligomerisation"]``. To compare an attribute with a value, compare
    its elements: ``"catalysis" in d["reaction_type"]``.

    Parameters
    ----------
    g : networkx.MultiDiGraph
        PSS interaction network or gene network.
    split_on_attrs : str or list of str, optional
        Edge attributes (e.g. ``["interaction"]``) that must not be merged away: parallel
        edges are only merged with others that have the same values for all of them. These
        keep their one value (not a list).
    verbose : bool
        Log every merged edge whose summary attributes differed, with the value kept (INFO;
        otherwise DEBUG). A summary, with the number of merged edges per differing attribute,
        is always logged (INFO).

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
            new_g.add_edge(source, target, **_merge_edge_data(group_edges, differing, log_level,
                                                              single=split_on_attrs))

    if differing:
        logger.info("Merged edges with differing summary values: %s. See the docstring for the "
                    "value kept, and verbose=True for details.",
                    ", ".join(f"{k} ({n})" for k, n in differing.most_common()))

    return new_g


# signs of a chain A -> B -> C
_COMPOSED_INTERACTION = {
    ("positive-influence", "positive-influence"): "positive-influence",
    ("positive-influence", "negative-influence"): "negative-influence",
    ("negative-influence", "positive-influence"): "negative-influence",
    ("negative-influence", "negative-influence"): "positive-influence",
}


def _propagating(e):
    '''Edge `e` with the reactions rewiring may go through (None if none): not the mutual edges
    between binding partners (partner <-> partner), only the complex-forming ones
    (partner -> complex).'''
    return _edge_subset(e, {r for r, a in _per_reaction(e).items()
                            if not (a.get("reaction_type") == "binding/oligomerisation"
                                    and a.get("target_role") != "product")})


def _catalysis_reactions(e, source_role, target_role):
    '''The reactions of `e` in which it is a catalysis edge from `source_role` to `target_role`
    (none for the reactions of a rewired edge, in which its source or its target takes no part).'''
    return {r for r, a in _per_reaction(e).items()
            if a.get("reaction_type") == "catalysis" and a.get("source_role") == source_role
            and a.get("target_role") == target_role}


def _consumes(e):
    '''The reactions of `e` in which it is a catalysis consumption edge, catalyst -| substrate (the
    enzyme uses up its substrate): not rewired through as it is; see :func:`_feeds`.'''
    return _catalysis_reactions(e, "catalyst", "substrate")


def _catalysis(e):
    '''The reactions of `e` in which it is a catalysis edge catalyst -> product.'''
    return _catalysis_reactions(e, "catalyst", "product")


def _conversion(e):
    '''The reactions of `e` in which it is a catalysis conversion edge, substrate -> product.'''
    return _catalysis_reactions(e, "substrate", "product")


def _feeds(e):
    '''A catalysis consumption edge (catalyst -| substrate) read the other way, for rewiring through
    the substrate: substrate -> catalyst, positive (the substrate feeds the catalyst's reaction).'''
    data = {}
    for k, value in e.items():
        if k.startswith("source_"):
            data["target_" + k[len("source_"):]] = value
        elif k.startswith("target_"):
            data["source_" + k[len("target_"):]] = value
        else:
            data[k] = value
    data["interaction"] = "positive-influence"
    return data


def _rewired_edge(up, down, node):
    '''The edge replacing up (A -> node) and down (node -> C).'''
    up_reactions, down_reactions = _per_reaction(up), _per_reaction(down)
    reactions = {}
    for r in set(up_reactions) | set(down_reactions):
        a_up, a_down = up_reactions.get(r, {}), down_reactions.get(r, {})
        attrs = {}
        for k in set(a_up) | set(a_down):
            if k.startswith("source_"):   # A's, in the reactions of A -> node
                attrs[k] = a_up.get(k)
            elif k.startswith("target_"):  # C's, in the reactions of node -> C
                attrs[k] = a_down.get(k)
            else:
                attrs[k] = (a_up if r in up_reactions else a_down).get(k)
        attrs["note"] = f"rewired through {node}"
        reactions[r] = attrs
    summary = {
        "interaction": _COMPOSED_INTERACTION.get(
            (up.get("interaction"), down.get("interaction")), "unknown-influence"),
        "directed": True,
    }
    if "rank" in up or "rank" in down:
        summary["rank"] = merge_values("rank", [up.get("rank"), down.get("rank")])[0]
    return _edge_data(reactions, summary)


def remove_and_rewire(g, nodes, dry_run=False):
    '''Remove nodes, connecting each of their upstream nodes to each of their downstream nodes.

    Changes `g` in place. The nodes are removed one at a time, each against the graph as
    rewired so far, so chains of removed nodes (A -> B -> C -> D, removing B and C) are
    bridged (A -> D). They are removed upstream first (in a topological order of the network,
    nodes on a cycle by name), whatever the order they are given in, so the result doesn't
    depend on it.

    Rewiring doesn't go through the mutual edges between binding partners
    (partner <-> partner), in either direction, only through the complex-forming ones
    (partner -> complex): A -> B and B <-> D (binding) give no A -> D, but A -> B and
    B -> B|D give A -> B|D.

    Catalysis: the catalyst -| substrate edges (the enzyme uses up its substrate) are not
    rewired through as they are. Removing a substrate connects the nodes upstream of it to its
    catalysts instead, as substrate -> catalyst (positive); and not to the product of a
    catalysed reaction, which it reaches through the catalyst (the substrate -> product edge of
    a reaction with a catalyst is not rewired through, from the substrate's side; nor from the
    product's, when the substrate is removed too, so the result doesn't depend on the order).
    So removing the metabolites of a pathway connects its enzymes, in order: E1 catalysing
    A -> B and E2 catalysing B -> C give E1 -> E2 (through B; positive), without shortcuts
    past E2. A reaction without a catalyst
    keeps its substrate -> product route (E1 -> B, B -> C uncatalysed, C -> D by E3: E1 -> E3).
    The edges are recognised by their ``reaction_type``, roles and reaction ids, per reaction:
    each reaction of a merged edge (e.g. from :func:`simplify_pss`) is rewired as its own.

    A new edge A -> C, replacing A -> B -> C:

    - ``reaction_id``: the reaction ids of both edges (a sorted list), and the other attributes
      per reaction, as in :func:`simplify_pss`: ``source_*`` attributes from A -> B (None for
      the reactions of B -> C only, in which A takes no part), ``target_*`` attributes from
      B -> C (None for those of A -> B only), and the others from the edge of the reaction;
    - ``note``: ``"rewired through B"`` (for each reaction);
    - ``interaction``: the sign of the chain (positive and negative: negative; two negatives:
      positive; anything with an unknown influence: unknown);
    - ``directed``: True; ``rank``: the lower.

    If A -> C already exists, the new edge is merged into it, as in :func:`simplify_pss` (a
    reaction of both keeps the attributes of the existing edge).

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
    to_remove = set(nodes)
    # upstream first: the topological order of the condensation (a cycle's nodes by name)
    c = nx.condensation(g)
    rank = {n: (i, str(n)) for i, scc in enumerate(nx.topological_sort(c)) for n in c.nodes[scc]["members"]}
    nodes = sorted(dict.fromkeys(nodes), key=rank.get)

    h = g.copy() if dry_run else g
    reasons = {}
    changed = set()  # (u, v) edges added or changed

    for node in dict.fromkeys(nodes):
        # (edge into the node, edge out of it) by upstream / downstream node; a catalyst using
        # up the node (catalyst -| node) is downstream of it (node -> catalyst)
        upstream, downstream = {}, {}
        catalysed = set()       # reactions with a catalyst using up the node
        made = {r for u in h.predecessors(node) for r in _catalysis(h[u][node])}  # by a catalyst
        for u in h.predecessors(node):
            e = _propagating(h[u][node]) if u != node else None
            if e is None:
                continue
            drop = _consumes(e)
            if drop:
                downstream[u] = _feeds(_edge_subset(e, drop))
                catalysed.update(drop)
            # from a substrate also removed, of a catalysed reaction: it reaches the node through
            # the catalyst
            if u in to_remove:
                drop = drop | (_conversion(e) & made)
            e = _edge_subset(e, set(_per_reaction(e)) - drop)
            if e is not None:
                upstream[u] = e
        for v in h.successors(node):
            e = _propagating(h[node][v]) if v != node else None
            # the product of a catalysed reaction: reached through its catalyst
            if e is not None:
                e = _edge_subset(e, set(_per_reaction(e)) - (_conversion(e) & catalysed))
            if e is not None:
                downstream.setdefault(v, e)
        new_edges = [(u, v, _rewired_edge(up, down, node))
                     for u, up in upstream.items() for v, down in downstream.items() if u != v]

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
    removes one of the two (which one is arbitrary), e.g. for a less cluttered drawing. A merged
    edge (e.g. from :func:`simplify_pss`) is removed only if all its reactions are bindings
    that the opposite edge has too. Don't use it before directed path or neighbourhood
    analysis.

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
                reactions = _per_reaction(e)
                if any(a.get("reaction_type") != "binding/oligomerisation"
                       for a in reactions.values()):
                    continue
                # the same reactions in the other direction
                for d_key, e2 in edge_items(node, upstream_node):
                    if not set(reactions) <= set(_per_reaction(e2)):
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


# ---------------------------------------------------------------------------
# Metabolism and pathways
# ---------------------------------------------------------------------------

def metabolite_enzymes(g):
    '''The enzymes of the reactions that produce and that consume each metabolite.

    From the catalysis reactions: a catalyst's edge to a reaction's product (``target_role`` ``"product"``) and
    to its substrate (``"substrate"``), for each reaction of an edge (also of merged edges, e.g. from
    :func:`simplify_pss`).

    Parameters
    ----------
    g : networkx.MultiDiGraph or networkx.DiGraph
        PSS interaction network or gene network.

    Returns
    -------
    pandas.DataFrame
        One row per enzyme, metabolite and reaction: ``enzyme``, ``metabolite``, ``role`` (``"producer"`` or
        ``"consumer"``), ``reaction_id``, ``metabolite_type`` (``Metabolite`` or ``MetaboliteFamily``).
    '''
    _check_not_reaction_graph(g, "metabolite_enzymes")
    role = {"product": "producer", "substrate": "consumer"}
    rows = []
    for u, v, d in g.edges(data=True):
        if g.nodes[v].get("node_type") not in ("Metabolite", "MetaboliteFamily"):
            continue
        for rx, a in _per_reaction(d).items():
            if (a.get("reaction_type") == "catalysis" and a.get("source_role") == "catalyst"
                    and a.get("target_role") in role):
                rows.append((u, v, role[a["target_role"]], rx, g.nodes[v]["node_type"]))
    return (pd.DataFrame(rows, columns=["enzyme", "metabolite", "role", "reaction_id", "metabolite_type"])
            .drop_duplicates().reset_index(drop=True))


def pathway_gene_sets(g, attribute="all_pathways", genes_only=True):
    '''PSS pathways as node sets, e.g. for an enrichment analysis.

    Parameters
    ----------
    g : networkx.Graph
        PSS gene network (or interaction network).
    attribute : str
        Node attribute with the pathways (default ``"all_pathways"``; ``"pathway"`` for the main ones).
    genes_only : bool
        Only genes (nodes with a ``species``; default), not metabolites or complexes.

    Returns
    -------
    dict
        Pathway -> sorted list of nodes.
    '''
    sets = {}
    for n, d in g.nodes(data=True):
        if genes_only and not d.get("species"):
            continue
        for pw in as_list(d.get(attribute)) or []:
            sets.setdefault(pw, set()).add(n)
    return {pw: sorted(nodes, key=str) for pw, nodes in sorted(sets.items())}

