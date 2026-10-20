'''Add experimental data (e.g. logFC and p-values, or any table) to network nodes.'''

import logging
import math

import pandas as pd

logger = logging.getLogger(__name__)


def add_experimental_data(g, df, logfc_col, pvalue_col=None, cutoff=0.05,
                          match_attribute=None, prefix=None, copy=False):
    '''Add experimental values to the nodes of `g` as node attributes.

    Like the filter functions, this changes `g` in place (and returns it). With
    ``copy=True`` it annotates a copy and leaves `g` unchanged, e.g. to keep the
    unannotated graph, or to add several experiments to separate copies.

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
    copy : bool
        Annotate a copy of `g` instead of `g` itself (default False).

    Returns
    -------
    networkx.Graph
        `g` (with ``copy=True``: a copy of `g`). Each node with at least one matching row gets:

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
    h = g.copy() if copy else g
    n_matched = 0
    for n, data in h.nodes(data=True):
        for a in attributes:
            data.pop(a, None)
        matched = _matched_ids(n, data, match_attribute, rows)
        if not matched:
            continue

        n_matched += 1
        logfc, pvalue = min((rows[k] for k in matched), key=rank)
        data[f"{prefix}logFC"] = logfc
        if pvalue_col is not None:
            data[f"{prefix}pvalue"] = pvalue
            data[f"{prefix}significant"] = pvalue is not None and pvalue <= cutoff
        data[f"{prefix}matched_ids"] = matched

    _log_matches(g, df, match_attribute, n_matched)
    return h


def add_node_data(g, df, columns=None, match_attribute=None, prefix=None, rank_by=None,
                  ascending=False, copy=False):
    '''Add the columns of a table to the nodes of `g` as node attributes.

    The general version of :func:`add_experimental_data`, for data other than logFC and
    p-values, e.g. GWAS results summarised per gene. Matching works the same way: on the
    node id, or on a node attribute (a single value or a list).

    Parameters
    ----------
    g : networkx.Graph
        Network to annotate.
    df : pandas.DataFrame
        Data indexed by identifier (e.g. gene id), one row per identifier.
    columns : list of str, optional
        Columns of `df` to add. Default: all columns.
    match_attribute : str, optional
        Match `df` index values against this node attribute instead of the node id, e.g.
        ``"translation"`` for a CKN translated with
        :func:`~skm_tools.translate.integrate_translation_ckn` (the target-species gene).
        Several nodes can match the same identifier; each gets the same values.
    prefix : str, optional
        Prefix for the new attribute names, e.g. ``"gwas "``. Default: no prefix.
    rank_by : str, optional
        If a node matches several identifiers (a list attribute), use the row with the
        highest value in this column (lowest with ``ascending=True``). Default: the first
        matching row in `df` order.
    ascending : bool
        See `rank_by` (default False).
    copy : bool
        Annotate a copy of `g` instead of `g` itself (default False).

    Returns
    -------
    networkx.Graph
        `g` (with ``copy=True``: a copy of `g`). Each node with a matching row gets
        ``{prefix}{column}`` for each column, and ``{prefix}matched_ids``: the `df`
        identifiers that matched the node. Nodes without a match get none of these
        attributes: values from an earlier call with the same `prefix` and columns are
        removed. Missing (NaN) values are stored as None.

    Raises
    ------
    KeyError
        If a column (or `rank_by`) is not in `df`.
    ValueError
        If the index of `df` has duplicates: aggregate to one row per identifier first.

    Examples
    --------
    >>> import networkx as nx, pandas as pd
    >>> g = nx.DiGraph([("AT1G01010_PruarM.1G000100", "AT1G01020")])
    >>> g.nodes["AT1G01010_PruarM.1G000100"]["translation"] = "PruarM.1G000100"
    >>> df = pd.DataFrame({"maxLog10P": [5.2], "trait": ["Rust"]}, index=["PruarM.1G000100"])
    >>> h = add_node_data(g, df, match_attribute="translation", prefix="gwas ")
    >>> h.nodes["AT1G01010_PruarM.1G000100"]["gwas trait"]
    'Rust'
    '''
    prefix = prefix or ""
    columns = list(df.columns) if columns is None else list(columns)
    missing = [c for c in columns + ([rank_by] if rank_by is not None else [])
               if c not in df.columns]
    if missing:
        raise KeyError(f"Columns not in df: {missing}")
    if df.index.has_duplicates:
        dups = df.index[df.index.duplicated()].unique()
        raise ValueError(f"{len(dups)} identifiers are in df more than once (e.g. {dups[0]!r}); "
                         "aggregate to one row per identifier first.")

    rows = {idx: {c: _value(v) for c, v in zip(columns, values)}
            for idx, *values in df[columns].itertuples(name=None)}
    if rank_by is None:
        order = {idx: i for i, idx in enumerate(df.index)}
    else:
        ranked = df.sort_values(rank_by, ascending=ascending, na_position="last", kind="stable")
        order = {idx: i for i, idx in enumerate(ranked.index)}

    attributes = [f"{prefix}{c}" for c in columns] + [f"{prefix}matched_ids"]
    h = g.copy() if copy else g
    n_matched = 0
    for n, data in h.nodes(data=True):
        for a in attributes:
            data.pop(a, None)
        matched = _matched_ids(n, data, match_attribute, rows)
        if not matched:
            continue

        n_matched += 1
        best = min(matched, key=order.__getitem__)
        for c, v in rows[best].items():
            data[f"{prefix}{c}"] = v
        data[f"{prefix}matched_ids"] = matched

    _log_matches(g, df, match_attribute, n_matched)
    return h


def _matched_ids(n, data, match_attribute, rows):
    '''The identifiers in `rows` that match node `n` (by id or `match_attribute`), sorted'''
    keys = data.get(match_attribute) if match_attribute else n
    if keys is None:
        return []
    if isinstance(keys, (list, tuple, set)):
        return sorted({k for k in keys if k in rows}, key=str)
    return [keys] if keys in rows else []


def _log_matches(g, df, match_attribute, n_matched):
    if n_matched == 0:
        example = next((d.get(match_attribute) if match_attribute else n
                        for n, d in g.nodes(data=True)), None)
        logger.warning("No node matched the df index (%s, e.g. %r; df index e.g. %r).",
                       match_attribute or "node ids", example, next(iter(df.index), None))
    else:
        logger.info("Added data to %d of %d nodes.", n_matched, g.number_of_nodes())


def _value(x):
    '''NaN -> None, numpy scalars -> Python; lists and other containers as they are'''
    if isinstance(x, (list, tuple, set, dict)):
        return x
    if x is None or (isinstance(x, float) and math.isnan(x)) or pd.isna(x):
        return None
    return x.item() if hasattr(x, "item") else x
