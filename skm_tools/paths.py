'''Shortest paths between nodes of interest.'''

import networkx as nx

from .utils import to_node_list


def get_paths(g, sources, targets, directed=True, all_shortest=True, shortest_overall=False):
    '''Shortest paths from each source to each target.

    Parameters
    ----------
    g : networkx.Graph
        Graph to search, any networkx graph type.
    sources : node or iterable of nodes
        Start node(s) of the paths.
    targets : node or iterable of nodes
        End node(s) of the paths.
    directed : bool
        Follow edge direction (default). If False, search the undirected view of `g`.
        Ignored if `g` is undirected.
    all_shortest : bool
        Return all shortest paths between each source and target (default),
        or only one of them.
    shortest_overall : bool
        If True, for each target keep only the paths that are as short as the shortest
        path from any of the sources (so a target reached in 2 steps from one source
        and 4 from another only keeps the 2-step paths). Default False: keep the
        shortest paths from every source.

    Returns
    -------
    list of list
        Paths as lists of nodes, from source to target. Source/target pairs without a
        path, and sources or targets not in `g`, are skipped. A node that is both a
        source and a target gives no path.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("B", "C"), ("A", "D"), ("D", "C")])
    >>> sorted(get_paths(g, "A", "C"))
    [['A', 'B', 'C'], ['A', 'D', 'C']]
    '''
    search_g = g.to_undirected(as_view=True) if (g.is_directed() and not directed) else g

    sources = [s for s in to_node_list(sources) if s in g]
    targets = [t for t in to_node_list(targets) if t in g]

    paths = []
    for target in targets:
        target_paths = []
        for source in sources:
            if source == target:
                continue
            try:
                if all_shortest:
                    target_paths += list(nx.all_shortest_paths(search_g, source=source, target=target))
                else:
                    target_paths.append(nx.shortest_path(search_g, source=source, target=target))
            except nx.NetworkXNoPath:
                pass

        if shortest_overall and target_paths:
            shortest = min(len(p) for p in target_paths)
            target_paths = [p for p in target_paths if len(p) == shortest]

        paths += target_paths

    return paths


def path_edges(paths, g):
    '''The edges of `g` along each path.

    Parameters
    ----------
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`get_paths`.
    g : networkx.Graph
        Graph the paths are in.

    Returns
    -------
    list of tuple
        Unique edges, in path order: ``(u, v)``, or ``(u, v, key)`` for every parallel
        edge if `g` is a multigraph. For an undirected search on a directed graph
        (``get_paths(..., directed=False)``), a step u-v gives whichever of ``(u, v)``
        and ``(v, u)`` exist in `g`.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("B", "C")])
    >>> path_edges([["A", "B", "C"]], g)
    [('A', 'B'), ('B', 'C')]
    '''
    edges = {}
    for path in paths:
        for s, t in zip(path, path[1:]):
            for u, v in ((s, t), (t, s)) if g.is_directed() else ((s, t),):
                if not g.has_edge(u, v):
                    continue
                if g.is_multigraph():
                    for k in g[u][v]:
                        edges[(u, v, k)] = None
                else:
                    edges[(u, v)] = None
                if (u, v) == (s, t) and g.is_directed():
                    break  # the forward edge exists, don't add the reverse as well

    return list(edges)


def path_subgraph(paths, g):
    '''Edge-induced subgraph of `g` made of the edges along `paths` (a copy).

    Parameters
    ----------
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`get_paths`.
    g : networkx.Graph
        Graph the paths are in.

    Returns
    -------
    networkx.Graph
        Same type as `g`, with only the nodes and edges on the paths (and their attributes).
    '''
    return g.edge_subgraph(path_edges(paths, g)).copy()
