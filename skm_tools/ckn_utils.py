'''CKN-specific filtering and annotation queries.'''

import re
from collections import defaultdict
import networkx as nx
import pandas as pd

from .utils import lists_intersect, is_listlike, remove_isolate_nodes

ckn_ranks = [0, 1, 2, 3, 4]


def rank_counts(g):
    '''Count (and print) the edges of each CKN rank.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`skm_tools.load_networks.ckn_to_networkx`.

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
        CKN, e.g. from :func:`skm_tools.load_networks.ckn_to_networkx`. Changed in place.
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
        CKN, e.g. from :func:`skm_tools.load_networks.ckn_to_networkx`. Changed in place.
    node_types : list of str, optional
        Keep only nodes of these ``node_type`` values (e.g. ``"protein_coding"``, ``"metabolite"``).
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
            and pd.notna(data['species'])
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


def get_all_annotations(g, key):
    '''All values of a list-valued node attribute.

    Parameters
    ----------
    g : networkx.Graph
    key : str
        Node attribute holding a list (or None), e.g. ``"GMM"`` or ``"tissue"``.

    Returns
    -------
    set
    '''
    annots = {
        x
        for n, d in g.nodes(data=True) if d[key] is not None for x in d[key]
    }
    return annots


def get_nodes_by_annotation(g, gmm=None, children=True):
    '''Nodes with any of the given GMM (MapMan) annotations.

    Parameters
    ----------
    g : networkx.Graph
        CKN, with list-valued ``GMM`` node attributes.
    gmm : list of str
        GMM bins, e.g. ``["27.3"]`` or ``["27.3_RNA.regulation of transcription"]``.
    children : bool
        Also match the sub-bins of each bin (default True), e.g. ``27.3.1``, ``27.3.2``.

    Returns
    -------
    list
        Matching nodes (empty if `gmm` is not a list).
    '''

    if is_listlike(gmm):

        if children:
            # automatically extract children annotations
            all_gmms = get_all_annotations(g, "GMM")
            gmm = [x.split("_")[0] for x in gmm]
            gmm = sorted([
                x for x in all_gmms
                if any([re.match(fr"^{re.escape(a)}[\.|_]", x) for a in gmm])
            ])
            print(
                f"skm-tools: Also using children annotations. Complete list is now:",
                end="\n\t")
            print('\n\t'.join(gmm))

        return [
            n for n, d in g.nodes(data=True) if lists_intersect(d["GMM"], gmm)
        ]

    return []


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
