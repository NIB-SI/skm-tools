'''MapMan annotations of PSS and CKN nodes, and MapMan gene sets from Mercator files.

Both networks have a ``mapman`` node attribute: a list of MapMan bins (GoMapMan 2, MapMan4),
each as ``<bin code>_<full bin name>``, e.g.
``26.11.3.2.1_External stimuli response.pathogen.defense mechanisms...``.
'''

import csv
import gzip
import re
from collections import defaultdict

from .utils import as_list


def get_all_annotations(g, key="mapman"):
    '''All values of a list-valued node attribute.

    Parameters
    ----------
    g : networkx.Graph
    key : str
        Node attribute holding a list, a single value, or None, e.g. ``"mapman"`` (default)
        or ``"tissue"``. A single value (e.g. a string) counts as a one-item list.

    Returns
    -------
    set
    '''
    return {x for _, d in g.nodes(data=True) for x in as_list(d.get(key)) or []}


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
    codes = [b.split("_")[0] for b in as_list(bins)]
    if children:
        pattern = re.compile("^(" + "|".join(re.escape(c) for c in codes) + r")[._]")
    else:
        pattern = re.compile("^(" + "|".join(re.escape(c) for c in codes) + r")_")

    return [
        n for n, d in g.nodes(data=True)
        if any(pattern.match(a) for a in (d.get("mapman") or []))
    ]


def mercator_gene_sets(path, levels=None, min_size=1, max_size=None, exclude_bins=("35",),
                       id_to_gene=None):
    '''MapMan gene sets from a Mercator4 result file, e.g. for an enrichment analysis.

    A gene counts in the bins it is annotated with and in all their parent bins: ``1.1.1``
    also counts in ``1.1`` and ``1``.

    Parameters
    ----------
    path : str or pathlib.Path
        Mercator4 result file (tab-separated, with the columns ``BINCODE``, ``NAME`` and
        ``IDENTIFIER``, values in single quotes), optionally gzipped (``.gz``).
    levels : int or list of int, optional
        Bin levels to keep, e.g. ``[1, 2, 3]`` (``1`` is the top level). Default: all.
    min_size, max_size : int, optional
        Keep bins with at least `min_size` (default 1) and at most `max_size` genes.
    exclude_bins : list of str
        Bins to leave out, with their sub-bins. Default ``("35",)``: MapMan4's
        *not assigned*.
    id_to_gene : callable, optional
        Converts a Mercator identifier (usually a lower-case protein or transcript id, e.g.
        ``"pruarm.3g005300.t1.p1_1"``) to a gene id. Default: the identifier as it is.

    Returns
    -------
    dict
        ``<bin code>_<full bin name>`` (as in the ``mapman`` node attribute) -> sorted list
        of gene ids, ordered by bin code.

    Examples
    --------
    >>> import re
    >>> sets = mercator_gene_sets(  # doctest: +SKIP
    ...     "parm_Mercator4v7_results.txt", levels=[1, 2, 3], min_size=5,
    ...     id_to_gene=lambda i: re.sub(r"^pruarm\\.(\\d+)g(\\d+)\\..*", r"PruarM.\\1G\\2", i))
    '''
    levels = set(as_list(levels)) if levels is not None else None
    excluded = [str(b) for b in as_list(exclude_bins) or []]
    opener = gzip.open if str(path).endswith(".gz") else open

    names, genes = {}, defaultdict(set)
    with opener(path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE):
            code, name, identifier = (row[k].strip().strip("'") for k in ("BINCODE", "NAME", "IDENTIFIER"))
            if not code:
                continue
            names.setdefault(code, name)
            if identifier:
                gene = id_to_gene(identifier) if id_to_gene else identifier
                parts = code.split(".")
                for i in range(1, len(parts) + 1):
                    genes[".".join(parts[:i])].add(gene)

    def keep(code, members):
        return (not any(code == b or code.startswith(b + ".") for b in excluded)
                and (levels is None or code.count(".") + 1 in levels)
                and len(members) >= min_size
                and (max_size is None or len(members) <= max_size))

    def order(code):
        return [int(p) if p.isdigit() else p for p in code.split(".")]

    return {f"{code}_{names.get(code, '')}": sorted(genes[code])
            for code in sorted(genes, key=order) if keep(code, genes[code])}
