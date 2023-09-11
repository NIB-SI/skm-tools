import networkx as nx


def get_cutset(sources, targets, g):
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

    print(cut_value == sum(source_sink_graph.edges[u, v]["capacity"] for (u, v) in cutset)   )
    
    return sorted(cutset)