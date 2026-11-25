'''CKN (Comprehensive Knowledge Network): load and filter.'''

import logging
from collections import Counter
from pathlib import Path

import networkx as nx

from .annotations import get_nodes_by_mapman
from .skm_download_urls import CKN_EDGE_FILE, CKN_EDGE_URL, CKN_NODE_FILE, CKN_NODE_URL
from .utils import as_list, download_if_missing, read_skm_table, remove_isolate_nodes

logger = logging.getLogger(__name__)

# CKN files (v2.1): the format of the PSS exports (see skm_tools.utils.read_skm_table)
_CKN_LIST_COLUMNS = {"synonyms", "mapman", "tissue", "interactionSources"}
_CKN_BOOL_COLUMNS = {"directed"}
_CKN_INT_COLUMNS = {"rank", "isTFregulation"}

CKN_RANKS = (0, 1, 2, 3, 4)
'''The CKN edge ranks, from 0 (curated, from PSS) to 4 (only predicted).'''


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def ckn_to_networkx(
        edge_path=None,
        node_path=None,
        as_directed=True,
        create_using=nx.DiGraph,
        data_dir=".",
    ):
    ''' Load CKN (v2.1) to a networkx directed graph, including node attributes.

    Downloads the CKN files from skm.nib.si if they don't exist yet.

    Parameters
    ----------
    edge_path : str or pathlib.Path, optional
        Path to the edge file (default ``<data_dir>/AtCKN-v2.1-2026.10.tsv.gz``); if the
        file does not exist, it is downloaded from skm.nib.si, gzipped if the name ends in
        ``.gz``.
    node_path : str or pathlib.Path, optional
        Path to the node file (default ``<data_dir>/AtCKN-v2.1-2026.10_node-annot.tsv.gz``);
        downloaded as for `edge_path`.
    as_directed : bool
        Add undirected edges (``directed`` False, e.g. binding) in both directions, so
        directed path searches can use them either way (default True): A -> B (undirected)
        becomes A -> B and B -> A. If False, every edge is added once, as listed in the file.
        To leave out the undirected edges, use
        ``filter_ckn_edges(g, filter_function=lambda d: d["directed"])``.
    create_using : networkx graph class
        Graph class to create (default ``networkx.DiGraph``; ``networkx.MultiDiGraph``
        keeps parallel edges). CKN can have two edges between the same nodes, of different
        types (e.g. ABI4 binds RGL2, and activates its transcription): a DiGraph keeps the
        last one of each direction, a MultiDiGraph all of them.
    data_dir : str or pathlib.Path
        Folder for the default file names (default: the current folder).

    Returns
    -------
    networkx.DiGraph
        CKN (or the `create_using` type): node attributes ``node_type`` (the PSS class,
        e.g. ``PlantCoding``, ``Metabolite``), ``locus_type`` (genes and RNAs: the TAIR locus
        type, e.g. ``protein_coding``, ``mirna``), ``species``, ``TAIR``, ``display_label``,
        ``short_name``, ``synonyms``, ``description``, ``mapman``, ``note`` and ``tissue``; edge
        attributes ``interaction`` (``positive-influence``, ``negative-influence`` or
        ``unknown-influence``), ``directed``, ``rank``, ``effect``, ``type``, ``species``,
        ``isTFregulation`` and ``interactionSources``. ``synonyms``, ``mapman``, ``tissue`` and
        ``interactionSources`` are lists; empty values are None. Nodes of the node file
        without edges are left out.
    '''
    data_dir = Path(data_dir)
    edge_path = Path(edge_path) if edge_path is not None else data_dir / CKN_EDGE_FILE
    node_path = Path(node_path) if node_path is not None else data_dir / CKN_NODE_FILE

    edge_df = read_skm_table(download_if_missing(edge_path, CKN_EDGE_URL),
                             _CKN_LIST_COLUMNS, _CKN_BOOL_COLUMNS, _CKN_INT_COLUMNS)
    node_df = read_skm_table(download_if_missing(node_path, CKN_NODE_URL), _CKN_LIST_COLUMNS)
    if "locus_type" not in node_df.columns:
        raise ValueError(f"{node_path} is not a CKN v2.1 node file (no locus_type column); "
                         "older CKN versions are not supported.")

    g = create_using()
    g.add_edges_from(
        (data.pop("source"), data.pop("target"), data) for data in edge_df.to_dict("records")
    )
    node_df = node_df.set_index("id")
    nx.set_node_attributes(g, node_df[node_df.index.isin(g.nodes)].to_dict('index'))

    if as_directed:
        reverse = [(v, u, data) for u, v, data in g.edges(data=True)
                   if data["directed"] is False and not g.has_edge(v, u)]
        g.add_edges_from(reverse)

    return g


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def rank_counts(g):
    '''The number of edges of each CKN rank.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`.

    Returns
    -------
    dict
        Rank -> number of edges, for every rank in :data:`CKN_RANKS` (and any other rank
        in `g`). Edges without a ``rank`` are not counted.
    '''
    counts = Counter(d["rank"] for *_, d in g.edges(data=True) if d.get("rank") is not None)
    return {r: counts[r] for r in sorted(set(CKN_RANKS) | set(counts))}


TF_REGULATION = "transcription factor regulation"
'''The CKN edge ``type`` of a transcription factor regulating a gene.'''


def transcription_factors(g, mapman_bins=None):
    '''Transcription factors in CKN.

    The sources of the network's transcription factor regulation edges (``type`` ``"transcription factor
    regulation"``), and optionally the nodes annotated with MapMan bins, e.g. ``"14.5"`` (RNA
    biosynthesis.DNA-binding transcriptional regulation in MapMan4), for TFs without regulation edges in CKN.
    Only the edges in `g` count: filter the edges first (:func:`filter_ckn_edges`) to restrict to some ranks.
    Load CKN as a MultiDiGraph (``ckn_to_networkx(create_using=nx.MultiDiGraph)``): a DiGraph keeps one edge per
    node pair and direction, so a binding edge can hide a regulation edge between the same genes.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`.
    mapman_bins : str or list of str, optional
        MapMan bins whose nodes (and those of their sub-bins) also count as TFs; see
        :func:`skm_tools.annotations.get_nodes_by_mapman`.

    Returns
    -------
    set
    '''
    tfs = {u for u, _, d in g.edges(data=True) if d.get("type") == TF_REGULATION}
    n_edges = len(tfs)
    if mapman_bins is not None:
        tfs |= set(get_nodes_by_mapman(g, mapman_bins))
    logger.info("%d transcription factors (%d with TF regulation edges)", len(tfs), n_edges)
    return tfs


_SIGN = {"positive-influence": 1, "negative-influence": -1}


def tf_target_enrichment(g, degs, universe=None, min_targets=5):
    '''Transcription factors whose CKN targets are over-represented among differentially expressed genes.

    For each TF (source of transcription factor regulation edges in `g`), a one-sided hypergeometric test of its
    targets among the DEGs. The background is the universe genes with at least one incoming transcription factor
    regulation edge in `g`: genes no TF is known to regulate can't be hits, and would make every TF look enriched.
    Only the edges in `g` count: filter the ranks first (:func:`filter_ckn_edges`) if needed.
    Load CKN as a MultiDiGraph (``ckn_to_networkx(create_using=nx.MultiDiGraph)``): a DiGraph keeps one edge per
    node pair and direction, so a binding edge can hide a regulation edge between the same genes.

    The activity score is computed over the TF's DE targets with a signed edge (``interaction``
    positive-influence / negative-influence; all edges from the TF to the target must agree) and a direction:
    the mean of sign(edge) x direction(target), +1 if all are consistent with a more active TF, -1 with a less
    active one.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx` (DiGraph or MultiDiGraph).
    degs : dict or iterable
        The DEGs: node -> direction (+1 up, -1 down, 0 unknown or mixed), or just the nodes (no direction).
    universe : iterable, optional
        The genes that could be DE, e.g. those measured. Default: every node of `g`.
    min_targets : int
        Test only TFs with at least this many targets in the background (default 5).

    Returns
    -------
    pandas.DataFrame
        One row per tested TF, sorted by p: ``tf``, ``targets`` (in the background), ``hits`` (DE targets),
        ``up_hits``, ``down_hits``, ``expected``, ``fold``, ``p``, ``padj`` (Benjamini-Hochberg over the TFs),
        ``n_signed``, ``activity_score`` (NaN without signed hits), ``hit_genes`` (list). The background size and
        the number of DEGs in it are in ``attrs["background"]`` and ``attrs["n_degs"]``.
    '''
    import pandas as pd
    from .stats import bh_fdr, hypergeom_sf

    direction = dict(degs) if isinstance(degs, dict) else {d: 0 for d in degs}
    signs = {}
    for u, v, d in g.edges(data=True):
        if d.get("type") == TF_REGULATION and u != v:
            signs.setdefault(u, {}).setdefault(v, set()).add(_SIGN.get(d.get("interaction"), 0))
    regulated = {v for t in signs.values() for v in t}
    background = regulated & (set(universe) if universe is not None else set(g))
    hits_all = set(direction) & background
    rows = []
    for tf, tg in signs.items():
        targets = set(tg) & background
        if len(targets) < min_targets:
            continue
        hits = targets & hits_all
        signed = [next(iter(tg[h])) * direction[h] for h in hits if len(tg[h]) == 1 and 0 not in tg[h]
                  and direction[h] != 0]
        expected = len(targets) * len(hits_all) / len(background)
        rows.append(dict(tf=tf, targets=len(targets), hits=len(hits),
                         up_hits=sum(direction[h] > 0 for h in hits), down_hits=sum(direction[h] < 0 for h in hits),
                         expected=expected, fold=len(hits) / expected if expected else float("nan"),
                         p=hypergeom_sf(len(hits), len(background), len(targets), len(hits_all)),
                         n_signed=len(signed), activity_score=sum(signed) / len(signed) if signed else float("nan"),
                         hit_genes=sorted(hits, key=str)))
    cols = ["tf", "targets", "hits", "up_hits", "down_hits", "expected", "fold", "p", "padj", "n_signed",
            "activity_score", "hit_genes"]
    res = pd.DataFrame(rows, columns=[c for c in cols if c != "padj"])
    res["padj"] = bh_fdr(res.p.to_numpy())
    res = res[cols].sort_values(["p", "tf"]).reset_index(drop=True)
    res.attrs.update(background=len(background), n_degs=len(hits_all))
    logger.info("%d TFs tested (background %d genes, %d DEGs)", len(res), len(background), len(hits_all))
    return res


def filter_ckn_edges(g,
                     keep_edge_ranks=None,
                     keep_edge_types=None,
                     filter_function=None,
                     remove_isolates=True):
    '''Remove CKN edges, in place.

    The filters are combined: an edge is kept only if it passes all of them.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`. Changed in place.
    keep_edge_ranks : int or iterable of int, optional
        Keep only edges of these ranks (0: best supported, to 4). Edges without a
        ``rank`` attribute are kept.
    keep_edge_types : str or iterable of str, optional
        Keep only edges of these ``type`` values (e.g. ``"binding"``).
    filter_function : callable, optional
        Called with each edge's attribute dict; return True to keep the edge.
    remove_isolates : bool
        Also remove nodes without edges afterwards (default True); this includes nodes that
        had no edges before.

    Returns
    -------
    removed_edges : list
        The removed edges, ``(u, v)`` or, for a multigraph, ``(u, v, key)``.
    removed_nodes : dict
        Removed node -> ``"isolate"``.
    '''
    keep_edge_ranks = as_list(keep_edge_ranks)
    keep_edge_types = as_list(keep_edge_types)

    def keep(d):
        if keep_edge_ranks is not None and d.get("rank") is not None and d["rank"] not in keep_edge_ranks:
            return False
        if keep_edge_types is not None and d.get("type") not in keep_edge_types:
            return False
        return filter_function is None or filter_function(d)

    if g.is_multigraph():
        removed_edges = [(u, v, k) for u, v, k, d in g.edges(keys=True, data=True) if not keep(d)]
    else:
        removed_edges = [(u, v) for u, v, d in g.edges(data=True) if not keep(d)]
    g.remove_edges_from(removed_edges)

    removed_nodes = remove_isolate_nodes(g) if remove_isolates else {}

    logger.info("Removed %d edges and %d nodes from the network.",
                len(removed_edges), len(removed_nodes))
    return removed_edges, removed_nodes


def filter_ckn_nodes(g,
                     node_types=None,
                     species=None,
                     tissues=None,
                     remove_isolates=True):
    '''Remove CKN nodes, in place.

    Complexes with a removed component are removed too: if a component is removed, the
    complex can't exist.

    Parameters
    ----------
    g : networkx.Graph
        CKN, e.g. from :func:`ckn_to_networkx`. Changed in place.
    node_types : str or iterable of str, optional
        Keep only nodes of these ``node_type`` values (the PSS classes, e.g. ``"PlantCoding"``,
        ``"Metabolite"``).
    species : str or iterable of str, optional
        Keep only nodes of these species (e.g. ``"ath"``); nodes without a species
        (e.g. metabolites) are kept. Complexes of an Arabidopsis protein and a pathogen
        (species ``"ath/foreign"``, e.g. ``RISC - vsiRNA``) are removed with
        ``species="ath"``, as their foreign component is.
    tissues : str or iterable of str, optional
        Keep only nodes annotated with at least one of these tissues.
    remove_isolates : bool
        Also remove nodes without edges afterwards (default True); this includes nodes that
        had no edges before.

    Returns
    -------
    dict
        Removed node -> reason (``"wrong species"``, ``"wrong node type"``,
        ``"wrong tissue type"``, ``"complex component removed"`` or ``"isolate"``).
    '''
    node_types = as_list(node_types)
    species = as_list(species)
    tissues = as_list(tissues)
    og_size = g.number_of_nodes()
    reasons = {}

    if species:
        # nodes without a species (e.g. metabolites) are kept
        reasons.update({n: "wrong species" for n, data in g.nodes(data=True)
                        if data.get("species") is not None and data["species"] not in species})

    if node_types:
        reasons.update({n: "wrong node type" for n, data in g.nodes(data=True)
                        if n not in reasons and data.get("node_type") not in node_types})

    if tissues:
        reasons.update({n: "wrong tissue type" for n, data in g.nodes(data=True)
                        if n not in reasons
                        and not any(t in tissues for t in (data.get("tissue") or []))})

    # complexes containing a removed node: CKN complex ids are their components' short
    # names (or ids), joined with "|"
    removed_names = {g.nodes[n].get("short_name") or n for n in reasons}
    reasons.update({n: "complex component removed" for n in g.nodes()
                    if n not in reasons and isinstance(n, str)
                    and removed_names.intersection(n.split("|"))})

    g.remove_nodes_from(list(reasons))

    if remove_isolates:
        reasons.update({n: r for n, r in remove_isolate_nodes(g).items() if n not in reasons})

    logger.info("Removed %d nodes from the network.", og_size - g.number_of_nodes())
    return reasons
