'''Add experimental data (e.g. logFC and p-values) to network nodes.'''

import logging
import math

import pandas as pd

logger = logging.getLogger(__name__)


def add_experimental_data(g, df, logfc_col, pvalue_col=None, cutoff=0.05,
                          match_attribute=None, prefix=None):
    '''Copy of `g` with experimental values added as node attributes.

    Unlike the filter functions (which change the graph in place), this returns a new
    graph and leaves `g` unchanged: it's often useful to keep the unannotated graph, or to
    add several experiments to separate copies.

    Parameters
    ----------
    g : networkx.Graph
        Network to annotate.
    df : pandas.DataFrame
        Experimental results, indexed by identifier (e.g. gene id), with a logFC column
        and optionally a p-value column.
    logfc_col : str
        Column of `df` with the (log) fold change.
    pvalue_col : str, optional
        Column of `df` with the (adjusted) p-value.
    cutoff : float
        p-value at or below which a value is marked significant (default 0.05).
        Only used with `pvalue_col`.
    match_attribute : str, optional
        Match `df` index values against this node attribute instead of the node id.
        The attribute can be a single value or a list, e.g. ``"ath_homologues"`` to put
        gene-level data on PSS functional cluster nodes, or ``"TAIR"`` for CKN.
    prefix : str, optional
        Prefix for the new attribute names, to keep several experiments apart, e.g.
        ``"H2O2 10min "``. Default: no prefix.

    Returns
    -------
    networkx.Graph
        A copy of `g`. Each node with at least one matching row gets:

        - ``{prefix}logFC``: the logFC,
        - ``{prefix}pvalue``: the p-value (with `pvalue_col`),
        - ``{prefix}significant``: p-value <= `cutoff` (with `pvalue_col`),
        - ``{prefix}matched_ids``: the `df` identifiers that matched the node.

        If several rows match a node (e.g. several genes of one functional cluster, or an
        identifier that is in `df` more than once), the values are those of the most
        significant row: lowest p-value, then largest absolute logFC. Nodes without a match
        get none of these attributes: values from an earlier call with the same `prefix` are
        removed. Missing (NaN) values in `df` are stored as None.

    Notes
    -----
    Identifiers are matched exactly. Normalise the index of `df` first if it differs from
    the node ids, e.g. for CKN (upper-case TAIR loci, ``AT1G01010``) and a table with
    ``At1g01010`` or isoforms (``AT1G01010.1``)::

        df.index = df.index.str.upper().str.replace(r"\\.\\d+$", "", regex=True)

    A warning is logged if no node matches.

    Examples
    --------
    >>> import networkx as nx, pandas as pd
    >>> g = nx.DiGraph([("AT1G01010", "AT1G01020")])
    >>> df = pd.DataFrame({"logFC": [1.5], "padj": [0.01]}, index=["AT1G01010"])
    >>> h = add_experimental_data(g, df, "logFC", "padj")
    >>> h.nodes["AT1G01010"]["significant"]
    True
    '''
    prefix = prefix or ""

    missing = [c for c in (logfc_col, pvalue_col) if c is not None and c not in df.columns]
    if missing:
        raise KeyError(f"Columns not in df: {missing}")

    def rank(values):
        # sort key: most significant first (lowest p-value, then largest |logFC|); missing
        # values last
        logfc, pvalue = values
        return (pvalue is None, pvalue if pvalue is not None else 0,
                logfc is None, -abs(logfc) if logfc is not None else 0)

    # identifier -> its most significant row's (logFC, p-value)
    pvalues = df[pvalue_col] if pvalue_col is not None else [None] * len(df)
    rows = {}
    duplicated = set()
    for idx, logfc, pvalue in zip(df.index, df[logfc_col], pvalues):
        values = (_value(logfc), _value(pvalue))
        if idx in rows:
            duplicated.add(idx)
            values = min(rows[idx], values, key=rank)
        rows[idx] = values
    if duplicated:
        logger.info("%d identifiers are in df more than once; the most significant row is used.",
                    len(duplicated))

    attributes = [f"{prefix}{a}" for a in ("logFC", "pvalue", "significant", "matched_ids")]
    h = g.copy()
    n_matched = 0
    for n, data in h.nodes(data=True):
        for a in attributes:
            data.pop(a, None)
        keys = data.get(match_attribute) if match_attribute else n
        if keys is None:
            continue
        if isinstance(keys, (list, tuple, set)):
            matched = sorted({k for k in keys if k in rows}, key=str)
        else:
            matched = [keys] if keys in rows else []
        if not matched:
            continue

        n_matched += 1
        logfc, pvalue = min((rows[k] for k in matched), key=rank)
        data[f"{prefix}logFC"] = logfc
        if pvalue_col is not None:
            data[f"{prefix}pvalue"] = pvalue
            data[f"{prefix}significant"] = pvalue is not None and pvalue <= cutoff
        data[f"{prefix}matched_ids"] = matched

    if n_matched == 0:
        example = next((d.get(match_attribute) if match_attribute else n
                        for n, d in g.nodes(data=True)), None)
        logger.warning("No node matched the df index (%s, e.g. %r; df index e.g. %r).",
                       match_attribute or "node ids", example, next(iter(df.index), None))
    else:
        logger.info("Added data to %d of %d nodes.", n_matched, g.number_of_nodes())
    return h


def _value(x):
    '''NaN -> None, numpy scalars -> Python'''
    if x is None or (isinstance(x, float) and math.isnan(x)) or pd.isna(x):
        return None
    return x.item() if hasattr(x, "item") else x
