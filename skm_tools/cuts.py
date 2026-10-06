'''Minimum cuts between sets of nodes.'''

import networkx as nx


def get_cutset(sources, targets, g):
    '''Minimum edge cut separating `sources` from `targets` (max-flow / min-cut).

    The sources are joined to a super-source and the targets to a super-sink with
    unlimited capacity, and the minimum cut between them is computed with Edmonds-Karp.
    Prints the maximum flow (the total capacity of the cut edges).

    Parameters
    ----------
    sources : iterable
        Source nodes (nodes not in `g` are ignored).
    targets : iterable
        Target nodes (nodes not in `g`, or also in `sources`, are ignored).
    g : networkx.DiGraph
        Graph with a ``capacity`` edge attribute. Edges without it have infinite capacity
        in networkx, so set it on every edge (e.g. ``nx.set_edge_attributes(g, 1, "capacity")``
        to count edges). `g` is not changed. Multigraphs are not supported (for PSS, use
        :func:`skm_tools.pss.simplify_pss` first).

    Returns
    -------
    list of tuple
        The cut edges (u, v), sorted. Empty if no target can be reached from the sources.

    Raises
    ------
    networkx.NetworkXUnbounded
        If a path from a source to a target has only edges without a ``capacity``.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D"), ("D", "E")])
    >>> nx.set_edge_attributes(g, 1, "capacity")
    >>> get_cutset(["A"], ["E"], g)
    max_flow = 1
    [('D', 'E')]
    '''
    sources = {n for n in sources if n in g}
    targets = {n for n in targets if n in g} - sources

    # super-source and super-sink: new objects, so they can't clash with node names;
    # edges without a capacity have unlimited capacity
    source, sink = object(), object()
    h = g.copy()
    h.add_edges_from((source, n) for n in sources)
    h.add_edges_from((n, sink) for n in targets)
    h.add_nodes_from([source, sink])

    flow, (reachable, non_reachable) = nx.minimum_cut(h, source, sink,
                                                      flow_func=nx.flow.edmonds_karp)
    print(f"max_flow = {flow}")

    return sorted((u, v) for u in reachable for v in h[u] if v in non_reachable)
