from collections.abc import Sequence


def lists_intersect(l1, l2):
    if (l1 is None) or (l2 is None):
        return False
    return any(i in l1 for i in l2)


def is_listlike(obj):
    return isinstance(obj, Sequence) and not isinstance(obj, str)
