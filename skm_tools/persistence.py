'''Save networks to disk and load them again.

This module only reads and writes files; it keeps no state between calls.
'''

import json
import pickle
from pathlib import Path

import networkx as nx

from .serialize import graph_from_dict, graph_to_dict


FORMATS = ("pickle", "json", "graphml")


def save_graph(g, path, format="pickle"):
    '''Save a network to a file.

    Parameters
    ----------
    g : networkx.Graph
        Network to save.
    path : str or pathlib.Path
        File to write (overwritten if it exists).
    format : {"pickle", "json", "graphml"}
        - ``"pickle"`` (default): fastest, keeps everything (graph type, attribute types),
          but Python-only, and unsafe to load from untrusted sources.
        - ``"json"``: portable and safe, keeps list attributes; tuples and sets become
          lists (see :func:`skm_tools.serialize.graph_to_dict`).
        - ``"graphml"``: for other network tools (Cytoscape, Gephi, ...). GraphML has no
          lists or nulls, so list attributes are joined with ``";"`` and None values are
          left out; loading does not split them again.
    '''
    path = Path(path)
    _check_format(format)

    if format == "pickle":
        with open(path, "wb") as handle:
            pickle.dump(g, handle, protocol=pickle.HIGHEST_PROTOCOL)
    elif format == "json":
        with open(path, "w") as handle:
            json.dump(graph_to_dict(g), handle)
    else:
        nx.write_graphml(_graphml_safe(g), path)


def load_graph(path, format="pickle"):
    '''Load a network saved with :func:`save_graph`.

    Parameters
    ----------
    path : str or pathlib.Path
        File to read.
    format : {"pickle", "json", "graphml"}
        Format the file was saved in. Only load pickle files you trust: unpickling can
        run arbitrary code.

    Returns
    -------
    networkx.Graph
        The network, of the type it was saved as (for graphml: as declared in the file).
    '''
    path = Path(path)
    _check_format(format)

    if format == "pickle":
        with open(path, "rb") as handle:
            return pickle.load(handle)
    if format == "json":
        with open(path) as handle:
            return graph_from_dict(json.load(handle))
    multigraph = _graphml_is_multigraph(path)
    return nx.read_graphml(path, force_multigraph=multigraph)


def _check_format(format):
    if format not in FORMATS:
        raise ValueError(f"format must be one of {FORMATS}, not {format!r}.")


def _graphml_value(v):
    if isinstance(v, (list, tuple, set, frozenset)):
        return ";".join(str(x) for x in v if x is not None)
    if isinstance(v, (str, int, float, bool)):
        return v
    return str(v)


def _graphml_safe(g):
    '''Copy of g with only GraphML-compatible attribute values.'''
    h = g.copy()
    for _, data in h.nodes(data=True):
        for k in list(data):
            if data[k] is None:
                del data[k]
            else:
                data[k] = _graphml_value(data[k])
    edges = h.edges(keys=True, data=True) if h.is_multigraph() else h.edges(data=True)
    for *_, data in edges:
        for k in list(data):
            if data[k] is None:
                del data[k]
            else:
                data[k] = _graphml_value(data[k])
    return h


def _graphml_is_multigraph(path):
    '''Edge ids are written for multigraph edge keys; networkx doesn't record the type otherwise.'''
    with open(path) as handle:
        for line in handle:
            if "<edge " in line:
                return ' id="' in line
    return False
