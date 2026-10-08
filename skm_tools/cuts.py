'''Minimum cuts between sets of nodes.'''

import logging

import networkx as nx

from .utils import resolve_nodes

logger = logging.getLogger(__name__)


def get_cutset(g, sources, targets, return_flow=False):
    '''Minimum edge cut separating `sources` from `targets` (max-flow / min-cut).

    The sources are joined to a super-source and the targets to a super-sink with
    unlimited capacity, and the minimum cut between them is computed with Edmonds-Karp.
    The maximum flow (the total capacity of the cut edges) is logged (INFO).

    Parameters
    ----------
    g : networkx.DiGraph
        Graph with a ``capacity`` edge attribute. Edges without it have infinite capacity
        in networkx, so set it on every edge (e.g. ``nx.set_edge_attributes(g, 1, "capacity")``
        to count edges). `g` is not changed. Multigraphs are not supported (for PSS, use
        :func:`skm_tools.pss.simplify_pss` first).
    sources : node or iterable of nodes
        Source nodes (nodes not in `g` are ignored, with a warning).
    targets : node or iterable of nodes
        Target nodes (nodes not in `g` are ignored, with a warning; targets that are also
        sources are ignored).
    return_flow : bool
        Also return the maximum flow.

    Returns
    -------
    list of tuple
        The cut edges (u, v), sorted. Empty if no target can be reached from the sources.
        With `return_flow`, a tuple ``(cut edges, max flow)``.

    Raises
    ------
    networkx.NetworkXUnbounded
        If a path from a source to a target has only edges without a ``capacity``.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D"), ("D", "E")])
    >>> nx.set_edge_attributes(g, 1, "capacity")
    >>> get_cutset(g, ["A"], ["E"])
    [('D', 'E')]
    >>> get_cutset(g, ["A"], ["D"], return_flow=True)
    ([('B', 'D'), ('C', 'D')], 2)
    '''
    sources = set(resolve_nodes(g, sources, "sources")[0])
    targets = set(resolve_nodes(g, targets, "targets")[0]) - sources

    # super-source and super-sink: new objects, so they can't clash with node names;
    # edges without a capacity have unlimited capacity
    source, sink = object(), object()
    h = g.copy()
    h.add_edges_from((source, n) for n in sources)
    h.add_edges_from((n, sink) for n in targets)
    h.add_nodes_from([source, sink])

    flow, (reachable, non_reachable) = nx.minimum_cut(h, source, sink,
                                                      flow_func=nx.flow.edmonds_karp)
    logger.info("max_flow = %s", flow)

    cut = sorted((u, v) for u in reachable for v in h[u] if v in non_reachable)
    return (cut, flow) if return_flow else cut
