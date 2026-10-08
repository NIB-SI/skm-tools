'''Shared helpers.'''

import csv
import gzip
import logging
import os
import shutil
import tempfile
from collections.abc import Iterable
from pathlib import Path
from urllib.request import urlopen

import networkx as nx
import pandas as pd

logger = logging.getLogger(__name__)


def as_list(x):
    '''A single value or an iterable of values, as a list.

    Strings (and other non-iterables) are a single value; lists, sets, tuples, generators
    etc. are several. None stays None (e.g. "no filter"). Tuples are several values, so
    pass a tuple node id in a list.

    Parameters
    ----------
    x : object or iterable or None

    Returns
    -------
    list or None

    Examples
    --------
    >>> as_list("ath")
    ['ath']
    >>> as_list(("ath", "stu"))
    ['ath', 'stu']
    >>> as_list(None) is None
    True
    '''
    if x is None:
        return None
    if isinstance(x, (str, bytes)) or not isinstance(x, Iterable):
        return [x]
    return list(x)


def read_skm_table(path, list_columns=(), bool_columns=(), int_columns=()):
    '''Read an SKM network export table (PSS or CKN) to a DataFrame.

    The format: tab-separated (optionally gzipped), a header, no quoting, an empty cell for no
    value, lists joined with ``|`` (as in TAIR's files and GAF: names contain ``,``, e.g.
    AHK2,3,4, and gene symbols ``;``, e.g. PIP1;3).

    Parameters
    ----------
    path : str or pathlib.Path
    list_columns, bool_columns, int_columns : collection of str
        Columns to split into lists (empty entries become None, so lists of different
        columns stay aligned), to read as ``True``/``False``, and as integers.

    Returns
    -------
    pandas.DataFrame
        Of str, list, bool and int values, with None for empty cells.
    '''
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""],
                     quoting=csv.QUOTE_NONE)
    df = df.astype(object).where(df.notna(), None)

    def to_bool(x):
        if x is None:
            return None
        return {"True": True, "False": False}[x]

    for c in df.columns:
        if c in list_columns:
            df[c] = df[c].map(lambda x: [v if v else None for v in x.split("|")] if x is not None else None)
        elif c in bool_columns:
            df[c] = df[c].map(to_bool)
        elif c in int_columns:
            df[c] = df[c].map(lambda x: int(x) if x is not None else None)
    return df


def resolve_nodes(g, nodes, what="nodes"):
    '''Split `nodes` into those in `g` and those not, logging a warning for the missing ones.

    Parameters
    ----------
    g : networkx.Graph
    nodes : node or iterable of nodes
        As for :func:`as_list`.
    what : str
        What the nodes are, for the warning (e.g. ``"sources"``).

    Returns
    -------
    found : list
        The nodes in `g`, in the given order.
    missing : list
        The nodes not in `g`.
    '''
    nodes = as_list(nodes) or []
    found = [n for n in nodes if n in g]
    missing = [n for n in nodes if n not in g]
    if missing:
        shown = ", ".join(map(str, missing[:10])) + (", ..." if len(missing) > 10 else "")
        logger.warning("%d of %d %s not in the graph: %s", len(missing), len(nodes), what, shown)
    return found, missing


def remove_isolate_nodes(g):
    '''Remove nodes without edges from `g`, in place.

    All nodes without edges are removed, including those that had none to begin with.

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
    return {n: "isolate" for n in isolates}


def merge_values(attribute, values):
    '''The value to keep for an edge attribute when merging edges.

    Missing values (None) are ignored. The rule, for the other values:

    - all equal (``==``): that value;
    - ``interaction``: ``"unknown-influence"`` if they differ;
    - ``directed``: True if any is True;
    - ``rank``: the lowest (best supported) value;
    - anything else: the first value.

    Parameters
    ----------
    attribute : str
        The attribute name.
    values : list
        One value per merged edge, in a fixed order (e.g. by reaction id).

    Returns
    -------
    value : object
        The value to keep (None if all are None).
    differ : bool
        Whether the values (other than None) differed.

    Examples
    --------
    >>> merge_values("interaction", ["positive-influence", "negative-influence"])
    ('unknown-influence', True)
    >>> merge_values("rank", [2, 0])
    (0, True)
    '''
    values = [v for v in values if v is not None]
    if not values:
        return None, False
    if all(v == values[0] for v in values[1:]):
        return values[0], False
    if attribute == "interaction":
        return "unknown-influence", True
    if attribute == "directed":
        return any(v is True for v in values), True
    if attribute == "rank":
        return min(values), True
    return values[0], True


def download(url, path):
    '''Download `url` to `path`, compressed or not as the file name says.

    A name ending in ``.gz`` is saved gzipped, any other name uncompressed: the download is
    decompressed or compressed as needed (gzipped downloads are recognised by their
    ``Content-Type``, ``application/gzip``). The file is written under a temporary name
    and renamed once complete, so an interrupted download leaves no partial file.

    Parameters
    ----------
    url : str
    path : str or pathlib.Path

    Returns
    -------
    pathlib.Path
        `path`.
    '''
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s to %s", url, path)

    with urlopen(url) as response:
        gzipped = response.headers.get_content_type() in ("application/gzip", "application/x-gzip")
        want_gzip = path.suffix == ".gz"
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".part")
        try:
            with os.fdopen(fd, "wb") as raw:
                if gzipped == want_gzip:
                    shutil.copyfileobj(response, raw)
                elif gzipped:
                    with gzip.GzipFile(fileobj=response) as unzipped:
                        shutil.copyfileobj(unzipped, raw)
                else:
                    with gzip.GzipFile(fileobj=raw, mode="wb") as zipped:
                        shutil.copyfileobj(response, zipped)
            os.replace(tmp, path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
    return path


def download_if_missing(path, url):
    '''Download `url` to `path` (see :func:`download`) unless the file exists.'''
    path = Path(path)
    if not path.exists():
        download(url, path)
    return path


def contrast_color(color):
    '''The complementary colour of a hex colour (e.g. for a label on that colour).

    Parameters
    ----------
    color : str
        ``"#RRGGBB"``.

    Returns
    -------
    str
        ``"#RRGGBB"``, each channel inverted.

    Examples
    --------
    >>> contrast_color("#000000")
    '#FFFFFF'
    '''
    r, g, b = (int(color[i:i + 2], 16) for i in (1, 3, 5))
    return f"#{255 - r:02X}{255 - g:02X}{255 - b:02X}"


def to_graph_tool(g):
    """Convert a networkx graph to a graph-tool graph (structure only, no attributes).

    For algorithms that are too slow in networkx on CKN-sized graphs. Requires graph-tool
    (https://graph-tool.skewed.de), which is not a dependency of skm-tools. Based on
    https://bbengfort.github.io/2016/06/graph-tool-from-networkx/.

    Parameters
    ----------
    g : networkx.Graph

    Returns
    -------
    gtG : graph_tool.Graph
        With a vertex property ``id`` holding the networkx node id as a string.
    vertices : dict
        networkx node -> graph-tool vertex.
    """
    import graph_tool as gt

    gtG = gt.Graph(directed=g.is_directed())
    # graph-tool vertices are indices: keep the networkx node ids in a property
    gtG.vertex_properties['id'] = gtG.new_vertex_property('string')

    vertices = {}
    for node in g.nodes():
        v = gtG.add_vertex()
        vertices[node] = v
        gtG.vp["id"][v] = str(node)

    for src, dst in g.edges():
        gtG.add_edge(vertices[src], vertices[dst])

    return gtG, vertices
