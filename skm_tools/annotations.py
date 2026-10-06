'''MapMan annotations of PSS and CKN nodes.

Both networks have a ``mapman`` node attribute: a list of MapMan bins (GoMapMan 2, MapMan4),
each as ``<bin code>_<full bin name>``, e.g.
``26.11.3.2.1_External stimuli response.pathogen.defense mechanisms...``.
'''

import re

from .utils import to_node_list


def get_all_annotations(g, key="mapman"):
    '''All values of a list-valued node attribute.

    Parameters
    ----------
    g : networkx.Graph
    key : str
        Node attribute holding a list (or None), e.g. ``"mapman"`` (default) or ``"tissue"``.

    Returns
    -------
    set
    '''
    return {x for _, d in g.nodes(data=True) if d.get(key) is not None for x in d[key]}


def get_nodes_by_mapman(g, bins, children=True):
    '''Nodes annotated with any of the given MapMan bins.

    Parameters
    ----------
    g : networkx.Graph
        PSS or CKN, with list-valued ``mapman`` node attributes.
    bins : str or list of str
        MapMan bin codes, e.g. ``"26.11"`` or ``["26.11", "15.5"]``; values with the bin name
        (``"26.11_External stimuli response.pathogen"``) work too: only the code is used.
    children : bool
        Also match the sub-bins of each bin (default True): ``"26.11"`` matches ``26.11``,
        ``26.11.3``, ``26.11.3.2.1``, ... but not ``26.110``.

    Returns
    -------
    list
        Matching nodes.

    Examples
    --------
    >>> import networkx as nx
    >>> g = nx.Graph()
    >>> g.add_node("AT2G38470", mapman=["26.11.3.2.1_External stimuli response.pathogen..."])
    >>> g.add_node("ABA", mapman=None)
    >>> get_nodes_by_mapman(g, "26.11")
    ['AT2G38470']
    '''
    codes = [b.split("_")[0] for b in to_node_list(bins)]
    if children:
        pattern = re.compile("^(" + "|".join(re.escape(c) for c in codes) + r")[._]")
    else:
        pattern = re.compile("^(" + "|".join(re.escape(c) for c in codes) + r")_")

    return [
        n for n, d in g.nodes(data=True)
        if any(pattern.match(a) for a in (d.get("mapman") or []))
    ]
