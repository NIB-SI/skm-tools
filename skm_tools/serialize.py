'''Convert networks and analysis results to JSON-safe Python data (dicts and lists).

The analysis functions return networkx objects; these helpers are the separate, explicit
step for when plain data is needed (JSON files, web APIs, other languages).
'''

import math

import networkx as nx


def to_json_safe(x):
    '''Recursively convert a value to JSON-safe Python types.

    dict -> dict (keys as str), list/tuple/set -> list (sets sorted when possible),
    numpy scalars and arrays -> Python values and lists, NaN/inf -> None, and anything
    else that isn't str/int/float/bool/None -> ``str(x)``.

    Parameters
    ----------
    x : object
        Value to convert.

    Returns
    -------
    object
        A value that :func:`json.dumps` accepts.
    '''
    if x is None or isinstance(x, (str, bool, int)):
        return x
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {str(k): to_json_safe(v) for k, v in x.items()}
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

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.DiGraph([("A", "B", {"weight": 1})])
    >>> graph_to_dict(g)["edges"]
    [{'source': 'A', 'target': 'B', 'weight': 1}]
    '''
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
        Node ids and attribute values come back as their JSON types (e.g. tuples become
        lists, sets become sorted lists).
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


def path_to_dict(paths):
    '''Paths as JSON-safe data.

    Parameters
    ----------
    paths : list of list
        Paths as lists of nodes, e.g. from :func:`skm_tools.paths.get_paths`.

    Returns
    -------
    dict
        ``{"paths": [[node, ...], ...], "lengths": [int, ...]}``, where a path's length
        is its number of edges.
    '''
    paths = [to_json_safe(list(p)) for p in paths]
    return {"paths": paths, "lengths": [len(p) - 1 for p in paths]}
