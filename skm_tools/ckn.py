'''CKN (Comprehensive Knowledge Network): load and filter.'''

from collections import defaultdict
from pathlib import Path
from urllib.request import urlretrieve

import networkx as nx
import pandas as pd

from .skm_download_urls import CKN_EDGE_URL, CKN_NODE_URL
from .utils import remove_isolate_nodes


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def ckn_to_networkx(
        edge_path=None,
        node_path=None,
        add_reciprocal_edges=True,
        directed=False,
        create_using=nx.DiGraph
    ):
    ''' Load CKN to a networkx directed graph, including node attributes.

    Downloads the CKN files from skm.nib.si if they don't exist yet.

    Parameters
    ----------
    edge_path : str or pathlib.Path
        Path to the edge list file (tab-separated, optionally gzipped);
        if the file does not exist, it is downloaded from skm.nib.si (gzipped, with a
        ``.tsv.gz`` suffix).
    node_path : str or pathlib.Path
        Path to the node annotation file; downloaded as for `edge_path`.
    add_reciprocal_edges : bool
        Add the reverse of undirected edges (``directed`` False), so directed path
        searches can use them in both directions (default True). Not meant to be used
        together with `directed`. For example A -> B (undirected) becomes A -> B and B -> A.
    directed : bool
        Remove undirected edges, and the nodes left without edges (default False).
    create_using : networkx graph class
        Graph class to create (default ``networkx.DiGraph``; ``networkx.MultiDiGraph``
        keeps parallel edges).

    Returns
    -------
    networkx.DiGraph
        CKN (or the `create_using` type), in the CKN v2.0.1 format: node attributes
        ``node_type`` (the PSS class, e.g. ``PlantCoding``, ``Metabolite``), ``locus_type``
        (genes and RNAs: the TAIR locus type, e.g. ``protein_coding``, ``mirna``),
        ``species``, ``TAIR``, ``display_label``, ``short_name``, ``synonyms``,
        ``description``, ``mapman``, ``note`` and ``tissue``; edge attributes ``interaction``
        (``positive-influence``, ``negative-influence`` or ``unknown-influence``), ``directed``,
        ``rank``, ``effect``, ``type``, ``species``, ``isTFregulation`` and ``interactionSources``.
        ``synonyms``, ``mapman``, ``tissue`` and ``interactionSources`` are lists; empty values
        are None.

    Notes
    -----
    Files in the older CKN v2 format (``node_ID``, ``GMM``, ``full_name``, ``isDirected``, ...)
    are converted to the v2.0.1 attributes when loaded, and so are older node types
    (``protein_coding``, ``metabolite``, ``biotic``, ...): to the PSS class, with the locus
    type of genes and RNAs in ``locus_type``.
    '''
    edge_path = Path(edge_path)
    node_path = Path(node_path)

    if not edge_path.exists() and edge_path.suffix != '.gz':
        edge_path = edge_path.with_suffix(".tsv.gz")
    if not edge_path.exists():
        print(f"Attempting to download the edge list to {edge_path}.", end=" ")
        urlretrieve(CKN_EDGE_URL, edge_path)
        print("Success.")

    if not node_path.exists() and node_path.suffix != '.gz':
        node_path = node_path.with_suffix(".tsv.gz")
    if not node_path.exists():
        print(f"Attempting to download the node annotations to {node_path}.", end=" ")
        urlretrieve(CKN_NODE_URL, node_path)
        print("Success.")

    edge_df = _read_ckn_table(edge_path)
    node_df = _read_ckn_table(node_path)
    if "node_ID" in node_df.columns:
        node_df, edge_df = _ckn_v2_to_v2_0_1(node_df, edge_df)
    if "locus_type" not in node_df.columns:
        node_df = _ckn_pss_node_types(node_df)

    for df in (node_df, edge_df):
        for c in _CKN_LIST_COLUMNS.intersection(df.columns):
            df[c] = df[c].map(lambda x: [y.strip() for y in x.split("|")] if x is not None else None)
    for c in ("rank", "isTFregulation"):
        edge_df[c] = edge_df[c].map(lambda x: int(x) if x is not None else None)
    edge_df["directed"] = edge_df["directed"].map({"True": True, "False": False})

    g = create_using()
    g.add_edges_from(
        (data.pop("source"), data.pop("target"), data) for data in edge_df.to_dict("records")
    )
    node_df = node_df.set_index("id")
    nx.set_node_attributes(g, node_df[node_df.index.isin(g.nodes)].to_dict('index'))

    if add_reciprocal_edges:
        edges_to_add = []
        for u, v, data in g.edges(data=True):
            if (not data["directed"]) and (not g.has_edge(v, u)):
                edges_to_add.append((v, u, data))
        _ = g.add_edges_from(edges_to_add)

    if directed:
        to_remove = [(u, v) for u, v, d in g.edges(data=True) if not d["directed"]]
        g.remove_edges_from(to_remove)

        # remove isolates resulting from filtering
        isolates = list(nx.isolates(g))
        g.remove_nodes_from(isolates)

    return g


# CKN files: tab-separated, header, lists joined with "|", empty = no value
_CKN_LIST_COLUMNS = {"synonyms", "mapman", "tissue", "interactionSources"}

# CKN v2 effect -> v2.0.1 interaction (anything else is an unknown influence)
_CKN_EFFECT_TO_INTERACTION = {"act": "positive-influence", "inh": "negative-influence"}

# older node_type (TAIR locus types and CKN's own) -> PSS class, as in CKN v2.0.1
# (skm-ckn scripts/ckn_v2.0.1.py); the locus types of genes and RNAs are kept as locus_type
_CKN_LOCUS_TYPES = {
    "protein_coding": "PlantCoding", "transposable_element_gene": "PlantCoding",
    "pseudogene": "PlantPseudogene",
    "mirna": "PlantNonCoding", "antisense_long_noncoding_rna": "PlantNonCoding",
    "pre_trna": "PlantNonCoding", "other_rna": "PlantNonCoding",
    "small_nuclear_rna": "PlantNonCoding", "small_nucleolar_rna": "PlantNonCoding",
}
_CKN_OTHER_TYPES = {
    "metabolite": "Metabolite", "complex": "Complex", "process": "Process", "abiotic": "ForeignAbiotic",
}
# biotic nodes, by id, with the classes of these entities in PSS
_CKN_BIOTIC = {
    **dict.fromkeys(["bacteria_flg22", "virus_6K2", "virus_CI", "virus_CP", "virus_HC-Pro",
                     "virus_NIa-Pro", "virus_NIb", "virus_P1", "virus_P3", "virus_VPg"], "ForeignCoding"),
    **dict.fromkeys(["bacteria", "virus_PVY"], "ForeignEntity"),
    **dict.fromkeys(["virus_dsRNA", "virus_vsiRNA", "virus_me-vsiRNA"], "ForeignNonCoding"),
}


def _read_ckn_table(path):
    '''Read a CKN file (optionally gzipped) to a DataFrame of str, with None for empty.'''
    # "N/A" is the empty species of metabolites in the CKN v2 node file
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["", "N/A"])
    return df.astype(object).where(df.notna(), None)


def _ckn_v2_to_v2_0_1(node_df, edge_df):
    '''Convert CKN v2 tables (AtCKN-v2-2023.06) to the v2.0.1 columns.'''
    node_df = node_df.rename(columns={"node_ID": "id", "full_name": "description", "GMM": "mapman"})
    # v2 joined tissues with ","
    node_df["tissue"] = node_df["tissue"].map(lambda x: x.replace(",", "|") if x is not None else None)
    node_df["display_label"] = [s if s is not None else i for s, i in zip(node_df["short_name"], node_df["id"])]

    edge_df = edge_df.rename(columns={"isDirected": "directed"})
    edge_df["directed"] = edge_df["directed"].map({"1": "True", "0": "False"})
    edge_df["interaction"] = edge_df["effect"].map(lambda x: _CKN_EFFECT_TO_INTERACTION.get(x, "unknown-influence"))

    return node_df, edge_df


def _ckn_pss_node_types(node_df):
    '''Convert older node types (CKN v2, and v2.0.1 before 2026-10-06) to the PSS classes,
    with the TAIR locus type of genes and RNAs in a new ``locus_type`` column.'''
    def convert(node_id, node_type):
        if node_type in _CKN_LOCUS_TYPES:
            return _CKN_LOCUS_TYPES[node_type], node_type
        if node_type in _CKN_OTHER_TYPES:
            return _CKN_OTHER_TYPES[node_type], None
        if node_type == "biotic" and node_id in _CKN_BIOTIC:
            return _CKN_BIOTIC[node_id], None
        raise ValueError(f"No PSS class for CKN node {node_id} (node_type {node_type}).")

    converted = [convert(i, t) for i, t in zip(node_df["id"], node_df["node_type"])]
    node_df = node_df.copy()
    node_df["node_type"] = [c for c, _ in converted]
    node_df.insert(node_df.columns.get_loc("node_type") + 1, "locus_type", [l for _, l in converted])
    return node_df


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

ckn_ranks = [0, 1, 2, 3, 4]


def rank_counts(g):
    '''Count (and print) the edges of each CKN rank.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`.

    Returns
    -------
    collections.defaultdict
        Rank -> number of edges. Edges without a ``rank`` attribute are not counted.
    '''
    counts = defaultdict(int)
    for _, _, data in g.edges(data=True):
        if 'rank' in data:
            counts[data['rank']] += 1
    for i in range(len(ckn_ranks)):
        print(f"rank {i}:\t {counts[i]:,}")

    return counts

def filter_ckn_edges(g,
                     keep_edge_ranks=None,
                     keep_edge_types=None,
                     filter_function=None,
                     remove_isolates=True):
    '''Remove CKN edges, in place.

    The filters are combined: an edge is kept only if it passes all of them.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`. Changed in place.
    keep_edge_ranks : int or list of int, optional
        Keep only edges of these ranks (0: best supported, to 4). Edges without a
        ``rank`` attribute are kept.
    keep_edge_types : str or list of str, optional
        Keep only edges of these ``type`` values (e.g. ``"binding"``).
    filter_function : callable, optional
        Called with each edge's attribute dict; return True to keep the edge.
    remove_isolates : bool
        Also remove nodes left without edges (default True).
    '''
    # / easier in reverse, but less intuitive...
    og_size = g.number_of_edges()

    to_remove = set()

    if isinstance(keep_edge_ranks, int):
        keep_edge_ranks = [keep_edge_ranks]

    if isinstance(keep_edge_types, str):
        keep_edge_types = [keep_edge_types]

    if isinstance(keep_edge_ranks, list):
        to_remove = [(u, v) for u, v, d in g.edges(data=True, )
                     if 'rank' in d and not (d["rank"] in keep_edge_ranks)]
        g.remove_edges_from(to_remove)

    if isinstance(keep_edge_types, list):
        to_remove = [(u, v) for u, v, d in g.edges(data=True, )
                     if not (d["type"] in keep_edge_types)]
        g.remove_edges_from(to_remove)

    if filter_function is not None:
        to_remove = ([(u, v) for u, v, d in g.edges(data=True)
                      if not filter_function(d)])
        g.remove_edges_from(list(to_remove))

    # remove isolates due to filtering
    if remove_isolates:
        isolate_reasons = remove_isolate_nodes(g)

    now_size = g.number_of_edges()
    print(f"Removed {og_size - now_size} edges from network.")


def filter_ckn_nodes(g,
                     node_types=None,
                     species=None,
                     tissues=None,
                     remove_isolates=True):
    '''Remove CKN nodes, in place.

    Complexes with a removed component are removed too.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`. Changed in place.
    node_types : list of str, optional
        Keep only nodes of these ``node_type`` values (the PSS classes, e.g. ``"PlantCoding"``,
        ``"Metabolite"``).
    species : list of str, optional
        Keep only nodes of these species (e.g. ``["ath"]``); nodes without a species
        (e.g. metabolites) are kept.
    tissues : list of str, optional
        Keep only nodes annotated with at least one of these tissues.
    remove_isolates : bool
        Also remove nodes left without edges (default True).

    Returns
    -------
    dict
        Removed node -> reason (``"wrong species"``, ``"wrong node type"``,
        ``"wrong tissue type"``, ``"complex component removed"`` or ``"isolate"``).
    '''
    og_size = g.number_of_nodes()

    to_remove = set()
    reasons = {}

    if species:
        no_species = [
            n for n, data in g.nodes(data=True)
            # missing requested species
            if not (data['species'] in species)
            # but not a "nan" species (e.g. metabolites)
            and data['species'] is not None
        ]
        to_remove.update(no_species)
        reasons = {
            **reasons,
            **{
                n: "wrong species"
                for n in no_species if not n in reasons
            }
        }

    if node_types:
        # nodes not in keep_types
        wrong_type = [
            n for n, data in g.nodes(data=True)
            if not (data['node_type'] in node_types)
        ]
        to_remove.update(wrong_type)
        reasons = {
            **reasons,
            **{
                n: "wrong node type"
                for n in wrong_type if not n in reasons
            }
        }

    if tissues:
        wrong_tissue = [
            n for n, data in g.nodes(data=True) if (not data['tissue']) or (
                not (len([aa for aa in data['tissue'] if aa in tissues]) > 0))
        ]
        to_remove.update(wrong_tissue)
        reasons = {
            **reasons,
            **{
                n: "wrong tissue type"
                for n in wrong_tissue if not n in reasons
            }
        }

    # now remove complexes that contain any nodes to be removed entities
    as_components = []
    for n in to_remove:
        x = g.nodes(data=True)[n]
        if isinstance(x["short_name"], str):
            as_components.append(x['short_name'])
        else:
            as_components.append(n)

    def complex_to_remove(x):
        for n in x.split("|"):
            if n in as_components:
                return True

    complex_component_missing = [n for n in g.nodes() if complex_to_remove(n)]
    to_remove.update(complex_component_missing)
    reasons = {
        **reasons,
        **{
            n: "complex component removed"
            for n in complex_component_missing if not n in reasons
        }
    }

    # remove the nodes
    g.remove_nodes_from(to_remove)

    # remove isolates due to filtering
    if remove_isolates:
        isolate_reasons = remove_isolate_nodes(g)
        reasons = {**reasons, **{n:r for n, r in isolate_reasons.items() if not n in reasons}}

    now_size = g.number_of_nodes()
    print(f"Removed {og_size - now_size} nodes from network.")

    return reasons


def to_graph_tool(g):
    """Convert a networkx graph to a graph-tool graph (structure only, no attributes).

    Requires graph-tool (https://graph-tool.skewed.de), which is not a dependency of
    skm-tools. Based on https://bbengfort.github.io/2016/06/graph-tool-from-networkx/.

    Parameters
    ----------
    g : networkx.Graph

    Returns
    -------
    gtG : graph_tool.Graph
        With a vertex property ``id`` holding the networkx node id as a string.
    vertices : dict
        networkx node -> graph-tool vertex.
    """

    import graph_tool as gt

    # Phase 0: Create a directed or undirected graph-tool Graph
    gtG = gt.Graph(directed=g.is_directed())

    # Also add the node id: in NetworkX a node can be any hashable type, but
    # in graph-tool node are defined as indices. So we capture any strings
    # in a special PropertyMap called 'id' -- modify as needed!
    gtG.vertex_properties['id'] = gtG.new_vertex_property('string')

    # Phase 2: Actually add all the nodes and vertices with their properties
    # Add the nodes
    vertices = {}  # vertex mapping for tracking edges later
    for node, data in g.nodes(data=True):

        # Create the vertex and annotate for our edges later
        v = gtG.add_vertex()
        vertices[node] = v

        # Set the vertex properties, not forgetting the id property
        gtG.vp["id"][v] = str(node)

    # Add the edges
    for src, dst, data in g.edges(data=True):

        # Look up the vertex structs from our vertices mapping and add edge.
        e = gtG.add_edge(vertices[src], vertices[dst])

    # Done, finally!
    return gtG, vertices
