'''Minimum cuts between sets of nodes.'''

import networkx as nx


def get_cutset(sources, targets, g):
    '''Minimum edge cut separating `sources` from `targets` (max-flow / min-cut).

    The sources are joined to a super-source and the targets to a super-sink with
    effectively unlimited capacity, and the minimum cut between them is computed with
    Edmonds-Karp. Prints the maximum flow.

    Parameters
    ----------
    sources : iterable
        Source nodes (nodes not in `g` are ignored).
    targets : iterable
        Target nodes (nodes not in `g`, or also in `sources`, are ignored).
    g : networkx.DiGraph
        Graph with a ``capacity`` edge attribute. Edges without it have infinite capacity
        in networkx, so set it on every edge (e.g. ``nx.set_edge_attributes(g, 1, "capacity")``
        to count edges). `g` is not changed.

    Returns
    -------
    list of tuple
        The cut edges (u, v), sorted.
    '''
    source_sink_graph = g.copy()
    
    source_sink_graph.add_node("source")
    for node in sources:
        if source_sink_graph.has_node(node):
            c = 99999
    #         c = 0
    #         for n in [x for x in source_sink_graph.successors(node)]:
    #             c += source_sink_graph.get_edge_data(node, n)['capacity']        
            source_sink_graph.add_edge('source', node, capacity=c)

    source_sink_graph.add_node('sink')
    for node in targets:
        if source_sink_graph.has_node(node) and not (node in sources):
            c = 99999
    #         c = 0
    #         for n in [x for x in source_sink_graph.predecessors(node)]:
    #             c += source_sink_graph.get_edge_data(n, node)['capacity']                
            source_sink_graph.add_edge(node, 'sink', capacity=c)
    
    r = nx.flow.edmonds_karp(source_sink_graph, 'source', 'sink')
    print(f"max_flow = {r.graph['flow_value']}")    
    
    cut_value, partition = nx.minimum_cut(source_sink_graph, 'source', 'sink', flow_func=nx.flow.edmonds_karp)
    reachable, non_reachable = partition
    
    cutset = set()
    for u, nbrs in ((n, source_sink_graph[n]) for n in reachable):
        cutset.update((u, v) for v in nbrs if v in non_reachable)

    return sorted(cutset)