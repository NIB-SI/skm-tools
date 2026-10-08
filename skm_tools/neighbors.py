'''Neighbourhoods of nodes of interest.'''

from .utils import resolve_nodes


_DIRECTIONS = ("both", "out", "in")


def neighborhood_nodes(g, nodes, depth=1, direction="both"):
    '''Nodes within `depth` steps of `nodes`.

    Parameters
    ----------
    g : networkx.Graph
        Graph to search.
    nodes : node or iterable of nodes
        Start node(s). Nodes not in `g` are ignored, with a warning.
    depth : int or None
        Number of steps to expand: a positive integer (default 1: first neighbours), or
        None for no limit (everything reachable).
    direction : {"both", "out", "in"}
        For directed graphs: follow edges in both directions (default, as an undirected
        graph), only downstream (successors), or only upstream (predecessors).
        Ignored for undirected graphs.

    Returns
    -------
    dict
        Node -> its distance (number of steps) from the closest start node, for the start
        nodes (distance 0) and every node reached.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("B", "C"), ("D", "A")])
    >>> neighborhood_nodes(g, "A", depth=2, direction="out")
    {'A': 0, 'B': 1, 'C': 2}
    '''
    if direction not in _DIRECTIONS:
        raise ValueError(f"direction must be one of {_DIRECTIONS}, not {direction!r}.")
    if depth is not None and (isinstance(depth, bool) or not isinstance(depth, int) or depth < 1):
        raise ValueError(f"depth must be a positive integer or None (no limit), not {depth!r}.")

    def step(n):
        if not g.is_directed():
            return g.neighbors(n)
        if direction == "out":
            return g.successors(n)
        if direction == "in":
            return g.predecessors(n)
        return list(g.successors(n)) + list(g.predecessors(n))

    distances = dict.fromkeys(resolve_nodes(g, nodes)[0], 0)
    frontier = list(distances)
    d = 0
    while frontier and (depth is None or d < depth):
        d += 1
        next_frontier = []
        for n in frontier:
            for m in step(n):
                if m not in distances:
                    distances[m] = d
                    next_frontier.append(m)
        frontier = next_frontier

    return distances


def get_neighborhood(g, nodes, depth=1, direction="both", induced=True, distance_attr="distance"):
    '''Subgraph of `g` around `nodes` (a copy).

    Parameters
    ----------
    g : networkx.Graph
        Graph to search.
    nodes : node or iterable of nodes
        Start node(s). Nodes not in `g` are ignored, with a warning.
    depth : int or None
        Number of steps to expand (default 1: first neighbours; None: no limit).
    direction : {"both", "out", "in"}
        For directed graphs: follow edges in both directions (default), only downstream,
        or only upstream. See :func:`neighborhood_nodes`.
    induced : bool
        If True (default), node-induced: every edge of `g` between the nodes found.
        If False, only the edges followed while expanding (each node reached through
        the edges from the previous step), so e.g. two second neighbours that happen to
        interact aren't linked.
    distance_attr : str or None
        Node attribute for each node's number of steps from the closest start node
        (default ``"distance"``; None: don't add it).

    Returns
    -------
    networkx.Graph
        Same type as `g`, with node and edge attributes, and the distances (see
        `distance_attr`).

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("B", "C"), ("C", "D")])
    >>> sorted(get_neighborhood(g, "B").nodes())
    ['A', 'B', 'C']
    '''
    distances = neighborhood_nodes(g, nodes, depth=depth, direction=direction)

    if induced:
        h = g.subgraph(distances).copy()
    else:
        # an edge is "followed" if it goes one step outwards in the allowed direction
        def followed(u, v):
            if not (u in distances and v in distances):
                return False
            if not g.is_directed() or direction == "both":
                return abs(distances[u] - distances[v]) == 1
            if direction == "out":
                return distances[v] == distances[u] + 1
            return distances[u] == distances[v] + 1

        # the edges of the nodes found (directed: their out-edges, which include every
        # edge between them)
        if g.is_multigraph():
            edges = [(u, v, k) for u, v, k in g.edges(list(distances), keys=True) if followed(u, v)]
        else:
            edges = [(u, v) for u, v in g.edges(list(distances)) if followed(u, v)]
        h = g.edge_subgraph(edges).copy()
        # keep the start nodes even if they have no edges
        h.add_nodes_from((n, g.nodes[n]) for n, d in distances.items() if d == 0)

    if distance_attr is not None:
        for n in h:
            h.nodes[n][distance_attr] = distances[n]

    return h
