'''E.g. filter PSS by species '''

import networkx as nx

from .utils import remove_isolate_nodes, unique_item


def remove_deadend_complexes(g):

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
    '''
    Inplace function
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
                 len([h for h in homologue_properties if data[h]])>0
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
        isolate_reasons = remove_isolate_nodes(g)
        reasons = {**reasons, **{n:r for n, r in isolate_reasons.items() if not n in reasons}}

    now_size = g.number_of_nodes()
    print(f"Removed {og_size - now_size} nodes from network.")

    return reasons


def simplify_pss(g):
    '''
    From MultiDiGraph to DiGraph by merging parallel edges
    Returns a new graph (not inplace function!)

    All reaction_ids are kept, for other edge attributes, a single value is kept, and a warning if
    multiple non-unique were observed.

    TODO - hierarchy for keeping attributes?
    '''

    new_g = nx.DiGraph()
    new_g.add_nodes_from(g.nodes(data=True))

    for source in g.nodes():
        edges_to_add = []
        for target in g[source]:
            edges = g[source][target]
            if len(edges) == 1:
                edges_to_add.append((source, target, edges[0]))
            else:
                data = {}
                reaction_ids = ",".join([d['reaction_id'] for d in edges.values()])
                data['reaction_id'] = reaction_ids
                for k in ['interaction_type', 'directed', 'reaction_type', 'reaction_effect', 'source_edge_type', 'source_location', 'source_form', 'target_edge_type', 'target_location', 'target_form']:
                    v, m = unique_item([d[k] for d in edges.values()])
                    if not (m is None):
                        print(f"{reaction_ids} --> {k}: {m}\n\tKeeping: {v}.")
                    data[k] = v
                edges_to_add.append((source, target, data))
        new_g.add_edges_from(edges_to_add)


    return new_g


def remove_and_rewire(g, nodes, dry_run=False):
    '''
    Removes nodes, replacing with edges from all upstream to all downstream nodes.

    Prints a summary of nodes that are removed, but no edges created to replace them.

    Use `dry_run=True` to get the summary, but not actually remove any of the nodes.

    For new edge attributes, all reaction_ids are kept, for other edge attributes, a single value is kept,
    NO WARNING if multiple non-unique were observed.

    TODO - hierarchy for keeping attributes?
    '''

    def generate_dict():
        return {"reaction_id": [], "reaction_type": [], "reaction_effect": [], "interaction_type": []}

    nodes_to_remove_and_rewire = set(nodes).intersection(g.nodes)

    all_new_edges = []

    if not isinstance(g, nx.DiGraph):
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
            if e["reaction_type"] == "binding/oligomerisation":

                # only keep as an upstream if source_edge_type = "SUBSTRATE" and target_edge_type == "PRODUCT"
                if not( e["source_edge_type"] == "SUBSTRATE" and e["target_edge_type"] == "PRODUCT" ):
                    failed_reasons[target] = "Not propagating binding"
                    failed_reasons[source] = "Not propagating binding"
                    continue

            upstream[upstream_node] = generate_dict()
            upstream[upstream_node]["node_type"] = g.nodes()[upstream_node]['node_type']
            upstream[upstream_node]["reaction_id"].append(e["reaction_id"])
            upstream[upstream_node]["reaction_type"].append(e["reaction_type"])
            upstream[upstream_node]["reaction_effect"].append(e["reaction_effect"])
            upstream[upstream_node]["interaction_type"].append(e["interaction_type"])

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
            downstream[downstream_node]["interaction_type"].append(e["interaction_type"])

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
                    "interaction_type": unique_item(downstream[target]["interaction_type"] + upstream[source]["interaction_type"])[0],
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
    '''Removes one if there are two binding edges between the same nodes from a "binding/oligomerisation" reaction
    Which edge is removed is arbitrary.

    Do not run this before doing directed paths/neighbour analysis.
    '''

    edges_to_remove = set()
    for node in g.nodes():

        # does not matter if we check upstream or downstream first
        for upstream_node in g.predecessors(node):
            e = g[upstream_node][node]
            if e["reaction_type"] == "binding/oligomerisation":
                # check if it is also downstream
                if upstream_node in g[node] and g[node][upstream_node]["reaction_id"] == e["reaction_id"]:
                    if not (node, upstream_node) in edges_to_remove:
                        edges_to_remove.update([(upstream_node, node)]) # keep this "first" one, but make sure we're not removing both!


    og_size_e = g.number_of_edges()

    g.remove_edges_from(edges_to_remove)

    now_size = g.number_of_edges()
    print(f"Removed {og_size_e - now_size} edges from network.")
