'''E.g. filter ckn on edge ranks'''

from collections import defaultdict
import networkx as nx

ckn_ranks = [0, 1, 2, 3, 4]

def rank_counts(g):
    counts = defaultdict(int)
    for _, _, data in g.edges(data=True):
        counts[data['rank']] += 1
    for i in range(len(ckn_ranks)):
        print(f"rank {i}:\t {counts[i]:,}")

    return counts

def filter_ckn_edges(g, keep_edge_ranks, remove_isolates=True):

    '''
    Inplace function

    '''
    # / easier in reverse, but less intuitive...
    og_size = g.number_of_edges()

    if isinstance(keep_edge_ranks, int):
        keep_edge_ranks = [keep_edge_ranks]

    to_remove_ranks = []
    for r in ckn_ranks:
        if not (r in keep_edge_ranks):
            to_remove_ranks.append(r)

    to_remove = []
    for r in to_remove_ranks:
        to_remove += [(u,v) for u, v, d in g.edges(data=True,) if d["rank"]==r]

    g.remove_edges_from(to_remove)

    # remove isolates resulting from filtering
    if remove_isolates:
        isolates = list(nx.isolates(g))
        g.remove_nodes_from(isolates)

    now_size = g.number_of_edges()
    print(f"Removed {og_size - now_size} edges from network.")

def filter_ckn_nodes(
        g,
        node_types=None,
        species=None,
        tissues=None,
        remove_isolates=True
    ):
    '''
    Inplace function
    '''
    og_size = g.number_of_nodes()

    to_remove = set()
    reasons = {}

    if species:
        no_species = [
            n for n, data in g.nodes(data=True) if not (
                data['species'] in species
            )
        ]
        to_remove.update(no_species)
        reasons = {**reasons, **{n:"wrong species" for n in no_species if not n in reasons}}

    if node_types:
        # nodes not in keep_types
        wrong_type = [n for n, data in g.nodes(data=True) if not (data['node_type'] in node_types)]
        to_remove.update(wrong_type)
        reasons = {**reasons, **{n:"wrong node type" for n in wrong_type if not n in reasons}}

    if tissues:
        wrong_tissue = [n for n, data in g.nodes(data=True) if ( not data['tissue'] ) or ( not (len([aa for aa in data['tissue'] if aa in tissues])>0) )]
        to_remove.update(wrong_tissue)
        reasons = {**reasons, **{n:"wrong tissue type" for n in wrong_tissue if not n in reasons}}


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
    reasons = {**reasons, **{n:"complex component removed" for n in complex_component_missing if not n in reasons}}

    # remove the nodes
    g.remove_nodes_from(to_remove)

    # remove isolates due to filtering
    if remove_isolates:
        isolates = list(nx.isolates(g))
        g.remove_nodes_from(isolates)
        reasons = {**reasons, **{n:"isolate" for n in isolates if not n in reasons}}

    now_size = g.number_of_nodes()
    print(f"Removed {og_size - now_size} nodes from network.")

    return reasons
