from collections.abc import Sequence
import networkx as nx

def lists_intersect(l1, l2):
    if (l1 is None) or (l2 is None):
        return False
    return any(i in l1 for i in l2)


def is_listlike(obj):
    return isinstance(obj, Sequence) and not isinstance(obj, str)

def to_list(x):
    if isinstance(x, str):
        return x.split(',')
    return None

def remove_isolate_nodes(g):
    '''remove isolate nodes from g'''

    isolates = list(nx.isolates(g))
    g.remove_nodes_from(isolates)
    reasons = {n:"isolate" for n in isolates}

    return reasons

# helper
def unique_item(x):
    ''' Return one unique item from a list, and a warning if non unique'''
    m = None
    s = set(x)
    if len(s) > 1:
        m = f"Multiple values in {x}."
    return next(iter(s)), m
