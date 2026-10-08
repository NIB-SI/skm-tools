'''Convert networks and analysis results to JSON-safe Python data (dicts and lists).

The analysis functions return networkx objects; these helpers are the separate, explicit
step for when plain data is needed (JSON files, web APIs, other languages).
'''

import math

import networkx as nx
import pandas as pd


def to_json_safe(x):
    '''Recursively convert a value to JSON-safe Python types.

    dict -> dict (keys as str), list/tuple/set -> list (sets sorted when possible),
    numpy scalars and arrays -> Python values and lists, NaN/inf and pandas' ``NA``/``NaT``
    -> None, and anything else that isn't str/int/float/bool/None -> ``str(x)``.

    Parameters
    ----------
    x : object
        Value to convert.

    Returns
    -------
    object
        A value that :func:`json.dumps` accepts.

    Raises
    ------
    ValueError
        If two keys of a dict are the same as strings (e.g. ``1`` and ``"1"``).
    '''
    if x is None or x is pd.NA or x is pd.NaT:
        return None
    if isinstance(x, (str, bool, int)):
        return x
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        safe = {str(k): to_json_safe(v) for k, v in x.items()}
        if len(safe) < len(x):
            raise ValueError(f"Keys that are the same as strings: {sorted(map(repr, x))}")
        return safe
    if isinstance(x, (set, frozenset)):
        try:
            x = sorted(x)
        except TypeError:
            x = list(x)
        return [to_json_safe(v) for v in x]
    if isinstance(x, (list, tuple)):
        return [to_json_safe(v) for v in x]
    if hasattr(x, "tolist"):  # numpy arrays and scalars
        return to_json_safe(x.tolist())
    return str(x)


def graph_to_dict(g):
    '''A networkx graph as a JSON-safe dict.

    Parameters
    ----------
    g : networkx.Graph
        Any networkx graph type.

    Returns
    -------
    dict
        ``{"directed": bool, "multigraph": bool, "graph": {...},
        "nodes": [{"id": ..., <attributes>}, ...],
        "edges": [{"source": ..., "target": ..., ["key": ...,] <attributes>}, ...]}``.
        Edge keys are only included for multigraphs. Attribute values are converted with
        :func:`to_json_safe`. Read back with :func:`graph_from_dict`.

    Raises
    ------
    ValueError
        If a node id is not a string or an integer (e.g. a tuple, which JSON can't keep), or
        an attribute has the name of a field of the format: ``id`` (nodes), ``source``,
        ``target`` or ``key`` (edges). Rename such attributes first.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B", {"weight": 1})])
    >>> graph_to_dict(g)["edges"]
    [{'source': 'A', 'target': 'B', 'weight': 1}]
    '''
    bad_ids = [n for n in g if not isinstance(n, (str, int))]
    if bad_ids:
        raise ValueError(f"Node ids must be strings or integers for JSON, not e.g. {bad_ids[0]!r}.")
    _check_attribute_names((d for _, d in g.nodes(data=True)), ("id",), "node")
    _check_attribute_names((d for *_, d in g.edges(data=True)), ("source", "target", "key"), "edge")

    nodes = [{**to_json_safe(data), "id": to_json_safe(n)} for n, data in g.nodes(data=True)]

    if g.is_multigraph():
        edges = [
            {**to_json_safe(data), "source": to_json_safe(u), "target": to_json_safe(v), "key": to_json_safe(k)}
            for u, v, k, data in g.edges(keys=True, data=True)
        ]
    else:
        edges = [
            {**to_json_safe(data), "source": to_json_safe(u), "target": to_json_safe(v)}
            for u, v, data in g.edges(data=True)
        ]

    # put the identifying fields first, for readability
    nodes = [{"id": d.pop("id"), **d} for d in nodes]
    edges = [{k: d.pop(k) for k in ("source", "target", "key") if k in d} | d for d in edges]

    return {
        "directed": g.is_directed(),
        "multigraph": g.is_multigraph(),
        "graph": to_json_safe(g.graph),
        "nodes": nodes,
        "edges": edges,
    }


def _check_attribute_names(attribute_dicts, reserved, what):
    clashes = sorted({k for d in attribute_dicts for k in d if k in reserved})
    if clashes:
        raise ValueError(f"{what.capitalize()} attributes {clashes} clash with the fields of the "
                         f"JSON format ({', '.join(reserved)}); rename them first.")


def graph_from_dict(data, create_using=None):
    '''Rebuild a networkx graph from :func:`graph_to_dict` output.

    Parameters
    ----------
    data : dict
        As returned by :func:`graph_to_dict` (or read from its JSON).
    create_using : networkx graph class, optional
        Graph class to create. Default: chosen from ``data["directed"]`` and
        ``data["multigraph"]``.

    Returns
    -------
    networkx.Graph
        Attribute values come back as their JSON types (e.g. tuples become lists, sets
        become sorted lists).
    '''
    if create_using is None:
        create_using = {
            (False, False): nx.Graph, (True, False): nx.DiGraph,
            (False, True): nx.MultiGraph, (True, True): nx.MultiDiGraph,
        }[(bool(data["directed"]), bool(data["multigraph"]))]

    g = create_using()
    g.graph.update(data.get("graph", {}))
    for node in data["nodes"]:
        node = dict(node)
        g.add_node(node.pop("id"), **node)
    for edge in data["edges"]:
        edge = dict(edge)
        u, v = edge.pop("source"), edge.pop("target")
        if g.is_multigraph() and "key" in edge:
            g.add_edge(u, v, key=edge.pop("key"), **edge)
        else:
            edge.pop("key", None)
            g.add_edge(u, v, **edge)
    return g


def paths_to_dict(paths, g=None, node_attrs=("display_label",), edge_attrs=("interaction",)):
    '''Paths as JSON-safe data, optionally with node and edge attributes.

    Parameters
    ----------
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`skm_tools.paths.get_paths`.
    g : networkx.Graph, optional
        The graph of the paths: to add the `node_attrs` of the nodes and the `edge_attrs` of
        each step.
    node_attrs : iterable of str
        Node attributes to add (with `g`; default ``display_label``).
    edge_attrs : iterable of str
        Edge attributes to add (with `g`; default ``interaction``).

    Returns
    -------
    dict
        ``{"paths": [[node, ...], ...], "lengths": [int, ...]}``, where a path's length
        is its number of edges. With `g`, also ``"nodes"``: ``{node: {attribute: value}}`` for
        the nodes on the paths, and ``"steps"``: per path, a list of
        ``{"source": u, "target": v, attribute: value}`` per edge. A step against the edge
        direction (paths from ``get_paths(..., directed=False)``) has the attributes of the
        reverse edge, and ``"reversed": True``. In a multigraph, an attribute's value is the
        list of the values of the parallel edges.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B", {"interaction": "positive-influence"})])
    >>> paths_to_dict([["A", "B"]], g, node_attrs=())["steps"]
    [[{'source': 'A', 'target': 'B', 'interaction': 'positive-influence'}]]
    '''
    paths = [list(p) for p in paths]
    result = {"paths": to_json_safe(paths), "lengths": [len(p) - 1 for p in paths]}
    if g is None:
        return result

    def step(u, v):
        d = {"source": u, "target": v}
        if not g.has_edge(u, v):
            u, v = v, u
            d["reversed"] = True
        edges = list(g[u][v].values()) if g.is_multigraph() else [g[u][v]]
        for a in edge_attrs:
            values = [e.get(a) for e in edges]
            d[a] = values if g.is_multigraph() else values[0]
        return d

    on_paths = dict.fromkeys(n for p in paths for n in p)
    result["nodes"] = to_json_safe({n: {a: g.nodes[n].get(a) for a in node_attrs} for n in on_paths})
    result["steps"] = to_json_safe([[step(u, v) for u, v in zip(p, p[1:])] for p in paths])
    return result
