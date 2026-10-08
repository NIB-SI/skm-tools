'''Shortest paths between nodes of interest.'''

import logging

import networkx as nx

from .utils import resolve_nodes

logger = logging.getLogger(__name__)


def _paths_to(target, source, pred):
    '''All shortest paths from `source` to `target` (a generator), backtracking through the
    shortest-path predecessors `pred` of a search from `source`.'''
    stack = [(target, [target])]
    while stack:
        node, path = stack.pop()
        if node == source:
            yield path[::-1]
            continue
        for p in reversed(pred[node]):
            stack.append((p, path + [p]))


def get_paths(g, sources, targets, directed=True, all_shortest=True, shortest_overall=False,
              max_paths=None):
    '''Shortest paths from each source to each target.

    Parameters
    ----------
    g : networkx.Graph
        Graph to search, any networkx graph type. Edge weights are not used: a path's
        length is its number of edges.
    sources : node or iterable of nodes
        Start node(s) of the paths. Nodes not in `g` are skipped, with a warning.
    targets : node or iterable of nodes
        End node(s) of the paths. Nodes not in `g` are skipped, with a warning.
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
    max_paths : int, optional
        Stop after this many paths, with a warning (default: no limit). The number of
        shortest paths can grow very fast in dense networks.

    Returns
    -------
    list of list
        Paths as lists of nodes, from source to target, by target and then by source (in
        the given order). Source/target pairs without a path are skipped. A node that is
        both a source and a target gives no path.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("B", "C"), ("A", "D"), ("D", "C")])
    >>> sorted(get_paths(g, "A", "C"))
    [['A', 'B', 'C'], ['A', 'D', 'C']]
    '''
    if max_paths is not None and (isinstance(max_paths, bool) or not isinstance(max_paths, int)
                                  or max_paths < 1):
        raise ValueError(f"max_paths must be a positive integer or None, not {max_paths!r}.")
    search_g = g.to_undirected(as_view=True) if (g.is_directed() and not directed) else g

    sources, _ = resolve_nodes(g, sources, "sources")
    targets, _ = resolve_nodes(g, targets, "targets")
    sources = list(dict.fromkeys(sources))

    def paths_to(target):
        """The shortest paths from the sources to `target` (a generator)."""
        if not all_shortest:
            # one path per pair: a bidirectional search per pair is fastest
            found = []
            for source in sources:
                if source != target:
                    try:
                        found.append(nx.shortest_path(search_g, source, target))
                    except nx.NetworkXNoPath:
                        pass
            if shortest_overall and found:
                found = [p for p in found if len(p) == min(map(len, found))]
            yield from found
            return
        # all paths: one breadth-first search per source (shortest-path predecessors and
        # distances), backtracking from each target
        reached = [s for s in sources if s != target and target in searches[s][1]]
        if shortest_overall and reached:
            shortest = min(searches[s][1][target] for s in reached)
            reached = [s for s in reached if searches[s][1][target] == shortest]
        for source in reached:
            yield from _paths_to(target, source, searches[source][0])

    searches = ({s: nx.predecessor(search_g, s, return_seen=True) for s in sources}
                if all_shortest else {})

    paths = []
    for target in dict.fromkeys(targets):
        for path in paths_to(target):
            if max_paths is not None and len(paths) >= max_paths:
                logger.warning("Stopped at max_paths=%d paths; there are more.", max_paths)
                return paths
            paths.append(path)

    return paths


def path_edges(g, paths):
    '''The edges of `g` along each path.

    Parameters
    ----------
    g : networkx.Graph
        Graph the paths are in.
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`get_paths`.

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
    >>> path_edges(g, [["A", "B", "C"]])
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


def path_subgraph(g, paths):
    '''Edge-induced subgraph of `g` made of the edges along `paths` (a copy).

    Parameters
    ----------
    g : networkx.Graph
        Graph the paths are in.
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`get_paths`.

    Returns
    -------
    networkx.Graph
        Same type as `g`, with only the nodes and edges on the paths (and their attributes).
    '''
    return g.edge_subgraph(path_edges(g, paths)).copy()
