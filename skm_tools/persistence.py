'''Save networks to disk and load them again.

This module only reads and writes files; it keeps no state between calls.
'''

import json
import pickle
from pathlib import Path

from .serialize import graph_from_dict, graph_to_dict


FORMATS = ("pickle", "json")


def save_graph(g, path, format="pickle"):
    '''Save a network to a file.

    Parameters
    ----------
    g : networkx.Graph
        Network to save.
    path : str or pathlib.Path
        File to write (overwritten if it exists).
    format : {"pickle", "json"}
        - ``"pickle"`` (default): fastest, keeps everything (graph type, attribute types),
          but Python-only, and unsafe to load from untrusted sources.
        - ``"json"``: portable and safe, keeps list attributes; tuples and sets become
          lists (see :func:`skm_tools.serialize.graph_to_dict`, also for its limits).

        On the full CKN (26k nodes, 899k edges), pickle took 2 s to save, 2 s to load and
        101 MB; JSON 15 s, 6 s and 193 MB. For other network tools (e.g. Cytoscape), see
        :mod:`skm_tools.cytoscape_utils`, or networkx's own writers.
    '''
    path = Path(path)
    _check_format(format)

    if format == "pickle":
        with open(path, "wb") as handle:
            pickle.dump(g, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, "w") as handle:
            json.dump(graph_to_dict(g), handle)


def load_graph(path, format="pickle"):
    '''Load a network saved with :func:`save_graph`.

    Parameters
    ----------
    path : str or pathlib.Path
        File to read.
    format : {"pickle", "json"}
        Format the file was saved in. Only load pickle files you trust: unpickling can
        run arbitrary code.

    Returns
    -------
    networkx.Graph
        The network, of the type it was saved as.
    '''
    path = Path(path)
    _check_format(format)

    if format == "pickle":
        with open(path, "rb") as handle:
            return pickle.load(handle)
    with open(path) as handle:
        return graph_from_dict(json.load(handle))


def _check_format(format):
    if format not in FORMATS:
        raise ValueError(f"format must be one of {FORMATS}, not {format!r}.")
