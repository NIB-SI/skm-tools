'''E.g. filter PSS by species '''

import networkx as nx

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
            n for n, data in g.nodes(data=True) if not (
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
        isolates = list(nx.isolates(g))
        g.remove_nodes_from(isolates)
        reasons = {**reasons, **{n:"isolate" for n in isolates if not n in reasons}}

    now_size = g.number_of_nodes()
    print(f"Removed {og_size - now_size} nodes from network.")

    return reasons



