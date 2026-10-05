'''Shared helpers.'''

from collections.abc import Iterable, Sequence

import networkx as nx


def lists_intersect(l1, l2):
    '''Whether two lists share any item (False if either is None).

    Parameters
    ----------
    l1, l2 : list or None

    Returns
    -------
    bool
    '''
    if (l1 is None) or (l2 is None):
        return False
    return any(i in l1 for i in l2)


def is_listlike(obj):
    '''Whether `obj` is a sequence other than a string.

    Parameters
    ----------
    obj : object

    Returns
    -------
    bool
    '''
    return isinstance(obj, Sequence) and not isinstance(obj, str)


def to_list(x):
    '''Split a comma-joined string into a list (None for anything else).

    Parameters
    ----------
    x : str or object

    Returns
    -------
    list of str or None
    '''
    if isinstance(x, str):
        return x.split(',')
    return None


def to_node_list(nodes):
    '''A single node or an iterable of nodes, as a list of nodes.

    Strings (and other non-iterables) are a single node; lists, sets, generators etc.
    are several. Tuples are treated as several nodes, so pass a tuple node id in a list.

    Parameters
    ----------
    nodes : node or iterable of nodes

    Returns
    -------
    list
    '''
    if isinstance(nodes, (str, bytes)) or not isinstance(nodes, Iterable):
        return [nodes]
    return list(nodes)


def remove_isolate_nodes(g):
    '''Remove nodes without edges from `g`, in place.

    Parameters
    ----------
    g : networkx.Graph

    Returns
    -------
    dict
        Removed node -> ``"isolate"`` (the removal reason, as in the filter functions).
    '''
    isolates = list(nx.isolates(g))
    g.remove_nodes_from(isolates)
    reasons = {n:"isolate" for n in isolates}

    return reasons


def unique_item(x):
    '''One item from a list, and a warning message if the items aren't all the same.

    Parameters
    ----------
    x : list

    Returns
    -------
    item : object
        One of the items (arbitrary if they differ).
    message : str or None
        None if all items are equal, otherwise a description of the values.
    '''
    m = None
    s = set(x)
    if len(s) > 1:
        m = f"Multiple values in {x}."
    return next(iter(s)), m
