'''PSS-specific filtering, simplification and rewiring.

For the PSS interaction network and gene networks loaded with
:mod:`skm_tools.load_networks`.
'''

import networkx as nx

from .utils import remove_isolate_nodes, unique_item


def remove_deadend_complexes(g):
    '''Remove complexes without outgoing edges, in place.

    A complex that influences nothing is a dead end in directed analyses. Repeated
    (up to five times) since removing one complex can leave another without outgoing edges.

    Parameters
    ----------
    g : networkx.DiGraph or networkx.MultiDiGraph
        PSS network. Changed in place.

    Returns
    -------
    list
        The removed complexes.
    '''
    removed_complexes = []

    # do five times
    for i in range(5):
        complexes = [n for n, data in g.nodes(data=True) if data["node_type"] == "Complex"]

        # g is a directed graph, so we can just use "neighbors" to find complexes without out/downstream edges
        to_remove = []
        for c in complexes:
            if len(list(g.neighbors(c))) == 0:
                to_remove.append(c)

        if len(to_remove) == 0:
            break

        g.remove_nodes_from(to_remove)
        removed_complexes += to_remove

    print(f"Number of complexes removed: {len(removed_complexes)}")

    return removed_complexes



def filter_pss_nodes(g, node_types=None, species=None, remove_isolates=True):
    '''Remove PSS nodes, in place.

    Complexes with a removed component are removed too (using the ``components``
    node attribute; components not in `g` are ignored).

    Parameters
    ----------
    g : networkx.Graph
        PSS interaction network or gene network. Changed in place.
    node_types : list of str, optional
        Keep only nodes of these ``node_type`` values (e.g. ``"PlantCoding"``, ``"Complex"``).
    species : list of str, optional
        Remove ``PlantCoding`` and ``PlantNonCoding`` nodes without homologues in any of
        these species (``<species>_homologues`` attributes, e.g. ``["stu"]``).
    remove_isolates : bool
        Also remove nodes left without edges (default True).

    Returns
    -------
    dict
        Removed node -> reason (``"species missing"``, ``"wrong node type"``,
        ``"complex component removed"`` or ``"isolate"``).
    '''
    og_size = g.number_of_nodes()

    to_remove = set()
    reasons = {}

    if species:
        # nodes that are PlantCoding or PlantNonCoding,
        # and do not have required homologues
        homologue_properties = [f"{sp}_homologues" for sp in species]
        no_species = [
            n for n, data in g.nodes(data=True) if not
            (
                 len([h for h in homologue_properties if data.get(h)])>0
                 or
                 not (data['node_type'] in ['PlantCoding', 'PlantNonCoding'])
            )
        ]
        to_remove.update(no_species)
        reasons = {**reasons, **{n:"species missing" for n in no_species if not n in reasons}}

    if node_types:
        # nodes not in keep_types
        wrong_type = [n for n, data in g.nodes(data=True) if not (data['node_type'] in node_types)]
        to_remove.update(wrong_type)
        reasons = {**reasons, **{n:"wrong node type" for n in wrong_type if not n in reasons}}

    # now remove complexes with a component that is in the network, but would no longer be.
    # Components are node ids (interaction network) or entity ids (gene network, where a
    # functional cluster is gone only once all its genes are).
    def entities(nodes):
        ids = set()
        for n in nodes:
            ids.add(n)
            ids.update(e for e in (g.nodes[n].get("entity") or []) if e)
        return ids

    before = entities(g.nodes())
    gone = before - entities(n for n in g.nodes() if n not in to_remove)
    complex_component_missing = [
        n for n, data in g.nodes(data=True)
        if gone.intersection(data.get("components") or [])
    ]
    to_remove.update(complex_component_missing)
    reasons = {**reasons, **{n:"complex component removed" for n in complex_component_missing if not n in reasons}}

    # remove the nodes
    g.remove_nodes_from(to_remove)

    # remove isolates due to filtering
    if remove_isolates:
        isolate_reasons = remove_isolate_nodes(g)
        reasons = {**reasons, **{n:r for n, r in isolate_reasons.items() if not n in reasons}}

    now_size = g.number_of_nodes()
    print(f"Removed {og_size - now_size} nodes from network.")

    return reasons


def simplify_pss(g, split_on_attrs=None):
    '''Merge parallel edges (one per reaction) into one edge per node pair.

    Returns a new graph; `g` is unchanged. For the interaction network or a gene network,
    not the reaction graph.

    Parameters
    ----------
    g : networkx.MultiDiGraph
        PSS interaction network or gene network.
    split_on_attrs : list of str, optional
        Edge attributes (e.g. ``["interaction"]``) that must not be merged away: parallel
        edges are only merged with others that have the same values for all of them.

    Returns
    -------
    networkx.DiGraph or networkx.MultiDiGraph
        A DiGraph, or a MultiDiGraph with `split_on_attrs` (as edges that differ on them
        stay separate). Merged edges get all ``reaction_id`` values, comma-joined; for every
        other attribute a single value is kept, with a printed warning when the merged
        edges had different values (e.g. a positive and a negative influence).

    Notes
    -----
    TODO - hierarchy for keeping attributes?
    '''

    split_on_attrs = split_on_attrs or []

    new_g = nx.MultiDiGraph() if split_on_attrs else nx.DiGraph()
    new_g.add_nodes_from(g.nodes(data=True))

    for source in g.nodes():
        edges_to_add = []
        for target in g[source]:
            edges = g[source][target]
            if len(edges) == 1:
                edges_to_add.append((source, target, next(iter(edges.values()))))
            else:
                if split_on_attrs:
                    groups = {}
                    for d in edges.values():
                        key = tuple(d[attr] for attr in split_on_attrs)
                        groups.setdefault(key, []).append(d)
                else:
                    groups = {None: list(edges.values())}

                for group_edges in groups.values():
                    if len(group_edges) == 1:
                        edges_to_add.append((source, target, group_edges[0]))
                        continue
                    data = {}
                    reaction_ids = ",".join([d['reaction_id'] for d in group_edges])
                    data['reaction_id'] = reaction_ids
                    # merge every other attribute present on any of the edges being combined,
                    # rather than a fixed whitelist, so nothing is silently dropped
                    other_attrs = set().union(*(d.keys() for d in group_edges)) - {'reaction_id'}
                    for k in sorted(other_attrs):
                        v, m = unique_item([d.get(k) for d in group_edges])
                        if not (m is None):
                            print(f"{reaction_ids} --> {k}: {m}\n\tKeeping: {v}.")
                        data[k] = v
                    edges_to_add.append((source, target, data))
        new_g.add_edges_from(edges_to_add)


    return new_g


def remove_and_rewire(g, nodes, dry_run=False):
    '''Remove nodes, connecting each of their upstream nodes to each downstream node.

    Changes `g` in place. Mutual binding edges (partner <-> partner) are not propagated,
    only the complex-forming ones (partner -> complex). Prints a summary of removed nodes
    for which no replacement edges were created.

    Parameters
    ----------
    g : networkx.DiGraph
        A simple directed graph, e.g. from :func:`simplify_pss`. Multigraphs are not supported.
    nodes : iterable
        Nodes to remove (nodes not in `g` are ignored).
    dry_run : bool
        Only print the summary; don't change `g`.

    Raises
    ------
    NotImplementedError
        If `g` is a multigraph or undirected.

    Notes
    -----
    New edges get all ``reaction_id`` values of the two edges they replace, comma-joined, and a
    single value of ``reaction_type``, ``reaction_effect`` and ``interaction`` (without a warning
    if they differ), plus a ``note``. TODO - hierarchy for keeping attributes?
    '''

    def generate_dict():
        return {"reaction_id": [], "reaction_type": [], "reaction_effect": [], "interaction": []}

    nodes_to_remove_and_rewire = set(nodes).intersection(g.nodes)

    all_new_edges = []

    if g.is_multigraph() or not g.is_directed():
        raise NotImplementedError("Currently only implemented for DiGraph, "
                                  "see pss_utils.simplify_pss.")

    # can remove without issue / rewiring
    in_pendants = [node for node in nodes_to_remove_and_rewire if g.in_degree(node) == 0]
    out_pendants = [node for node in nodes_to_remove_and_rewire if g.out_degree(node) == 0]

    for node in (remaining_nodes_to_remove := nodes_to_remove_and_rewire - set(in_pendants + out_pendants)):
        new_edges = []
        upstream = {}
        downstream = {}
        failed_reasons = {}

        for upstream_node in g.predecessors(node):

            if upstream_node in nodes_to_remove_and_rewire:
                failed_reasons[upstream_node] = "To be removed"
                upstream[upstream_node] = None
                continue

            e = g[upstream_node][node]

            # if binding interaction, only keep if this is the complex forming interaction
            # (partner -> complex, target_role "product"), not the mutual partner <-> partner one
            if e["reaction_type"] == "binding/oligomerisation":
                if e["target_role"] != "product":
                    failed_reasons[upstream_node] = "Not propagating binding"
                    continue

            upstream[upstream_node] = generate_dict()
            upstream[upstream_node]["node_type"] = g.nodes()[upstream_node]['node_type']
            upstream[upstream_node]["reaction_id"].append(e["reaction_id"])
            upstream[upstream_node]["reaction_type"].append(e["reaction_type"])
            upstream[upstream_node]["reaction_effect"].append(e["reaction_effect"])
            upstream[upstream_node]["interaction"].append(e["interaction"])

        for downstream_node in g.successors(node):

            if downstream_node in nodes_to_remove_and_rewire:
                failed_reasons[downstream_node] = "To be removed"
                downstream[downstream_node] = None
                continue

            e = g[node][downstream_node]




            downstream[downstream_node] = generate_dict()
            downstream[downstream_node]["node_type"] = g.nodes()[downstream_node]['node_type']
            downstream[downstream_node]["reaction_id"].append(e["reaction_id"])
            downstream[downstream_node]["reaction_type"].append(e["reaction_type"])
            downstream[downstream_node]["reaction_effect"].append(e["reaction_effect"])
            downstream[downstream_node]["interaction"].append(e["interaction"])

            if downstream_node in nodes_to_remove_and_rewire:
                failed_reasons[downstream_node] = "To be removed"

        for source in set(upstream):
            if source in nodes_to_remove_and_rewire:
                continue

            for target in set(downstream):
                if target in nodes_to_remove_and_rewire:
                    continue

                if source == target:
                    failed_reasons[target] = "Same source/target"
                    failed_reasons[source] = "Same source/target"
                    continue

                data = {
                    "reaction_id": ",".join(downstream[target]["reaction_id"] + upstream[source]["reaction_id"]),
                    "reaction_type": unique_item(downstream[target]["reaction_type"] + upstream[source]["reaction_type"])[0],
                    "reaction_effect": unique_item(downstream[target]["reaction_effect"] + upstream[source]["reaction_effect"])[0],
                    "interaction": unique_item(downstream[target]["interaction"] + upstream[source]["interaction"])[0],
                    "note": f"rewired {node} from {source} to {target}",
                }
                new_edges.append((source, target, data))

        all_new_edges += new_edges

        # print out any issues:
        if len(new_edges) == 0:
            # test if theres a problem
            if upstream == downstream:
                # don't need to connect self - no problem
                continue

            print(node)
            for n in downstream:
                if n in failed_reasons:
                    print(f"  Downstream: {n} -- {failed_reasons[n]}")
                else:
                    print(f"  Downstream {n}")
            print("     --->")
            for n in upstream:
                if n in failed_reasons:
                    print(f"  Upstream: {n} -- {failed_reasons[n]}")
                else:
                    print(f"  Upstream {n}")
            print()

    if dry_run:
        print("End of dry run.")
        return

    og_size_n = g.number_of_nodes()
    og_size_e = g.number_of_edges()

    g.remove_nodes_from(in_pendants)
    g.remove_nodes_from(out_pendants)
    step1_size_n = g.number_of_nodes()
    step1_size_e = g.number_of_edges()
    print(f"Removed {og_size_n - step1_size_n} pendant nodes, with {og_size_e - step1_size_e} edges.")

    g.remove_nodes_from(remaining_nodes_to_remove)
    step2_size_n = g.number_of_nodes()
    step2_size_e = g.number_of_edges()
    print(f"Removed further {step1_size_n - step2_size_n} nodes, with {step1_size_e - step2_size_e} edges.")

    g.add_edges_from(all_new_edges)
    step3_size_e = g.number_of_edges()
    print(f"Added {step3_size_e - step2_size_e} edges in rewiring")


def remove_duplicated_binding_edges(g):
    '''Keep one direction of each mutual binding edge, in place.

    Binding partners are linked in both directions (A -> B and B -> A, same reaction); this
    removes one of the two (which one is arbitrary), e.g. for a less cluttered drawing.
    Don't use it before directed path or neighbourhood analysis.

    Parameters
    ----------
    g : networkx.DiGraph or networkx.MultiDiGraph
        PSS network (each parallel edge of a multigraph is considered separately). Changed in place.
    '''

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

        # does not matter if we check upstream or downstream first
        for upstream_node in g.predecessors(node):
            for u_key, e in edge_items(upstream_node, node):
                if e["reaction_type"] != "binding/oligomerisation":
                    continue

                # check if it is also downstream
                for d_key, e2 in edge_items(node, upstream_node):
                    if e2["reaction_id"] != e["reaction_id"]:
                        continue

                    forward = (upstream_node, node, u_key) if is_multi else (upstream_node, node)
                    reverse = (node, upstream_node, d_key) if is_multi else (node, upstream_node)
                    if reverse not in edges_to_remove:
                        edges_to_remove.add(forward) # keep this "first" one, but make sure we're not removing both!


    og_size_e = g.number_of_edges()

    g.remove_edges_from(edges_to_remove)

    now_size = g.number_of_edges()
    print(f"Removed {og_size_e - now_size} edges from network.")
