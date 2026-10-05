'''Load the SKM networks (PSS and CKN) from file, or download them, as networkx graphs.'''

from urllib.request import urlretrieve
import csv
import gzip
from pathlib import Path
import networkx as nx
import pandas as pd
from .skm_download_urls import *


# PSS exports (pss-export): tab-separated, header, no quoting, empty = no value,
# lists joined with ";" (never split on ",": names such as AHK2,3,4 contain commas).
_PSS_LIST_SEPARATOR = ";"
_PSS_LIST_COLUMNS = {
    "synonyms",
    "all_pathways",
    "external_links",
    "components",
    "genes",
}
# Gene network only: a gene in several functional clusters has one entry per cluster
# (in the same order) in each of these. Loaded as lists for every node, so the type
# doesn't depend on the row. (display_label stays a string: it's what to show.)
_PSS_GENE_CLUSTER_COLUMNS = {
    "entity",
    "short_name",
    "pathway",
    "functional_cluster_id",
}
_PSS_BOOL_COLUMNS = {
    "directed",
    "location_putative",
    "source_location_putative",
    "target_location_putative",
}


def _split_pss_list(x):
    if x is None:
        return None
    # keep empty entries (as None) so per-cluster lists stay aligned
    return [v if v else None for v in x.split(_PSS_LIST_SEPARATOR)]


def _read_pss_table(path, list_columns=()):
    '''Read a PSS export table to a DataFrame of str/list/bool values, with None for empty.'''
    path = Path(path)
    if not path.exists():
        # TODO: download once the new exports are published on skm.nib.si
        raise FileNotFoundError(f"{path} not found. Download URLs for the PSS network exports "
                                "are not available yet; pass the path to a local export file.")

    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=[""],
                     quoting=csv.QUOTE_NONE)
    df = df.astype(object).where(df.notna(), None)

    list_columns = _PSS_LIST_COLUMNS | set(list_columns)
    for c in df.columns:
        if c in list_columns or c.endswith("_homologues"):
            df[c] = df[c].map(_split_pss_list)
        elif c in _PSS_BOOL_COLUMNS:
            df[c] = df[c].map(lambda x: {"True": True, "False": False}[x] if x is not None else None)

    return df


def _pss_export_to_networkx(edge_path, node_path, edge_key, node_list_columns=()):
    '''Build a MultiDiGraph from a PSS export: nodes (with attributes) from the node file,
    edges from the edge file, keyed by the `edge_key` column.'''
    edge_df = _read_pss_table(edge_path)
    node_df = _read_pss_table(node_path, node_list_columns)

    g = nx.MultiDiGraph()
    g.add_nodes_from(
        (data.pop("id"), data) for data in node_df.to_dict("records")
    )
    g.add_edges_from(
        (data.pop("source"), data.pop("target"), data[edge_key], data)
        for data in edge_df.to_dict("records")
    )
    return g


def pss_reaction_graph_to_networkx(edge_path, node_path):
    ''' Load the PSS reaction graph export to a networkx directed multigraph,
    including node attributes.

    Entities and reactions are both nodes (reactions have `node_type` == "reaction"), with one edge
    per reaction participant: participant -> reaction for inputs and modifiers,
    reaction -> participant for products. This is the lossless form of PSS.
    Edges are keyed by the participant's `role`, since an entity can take part in the same
    reaction twice (e.g. as template and stimulator).

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge file (pss-reaction-graph-edges-*.tsv)

    node_path : str or pathlib.Path
        Path to the node file (pss-reaction-graph-nodes-*.tsv)
    '''
    return _pss_export_to_networkx(edge_path, node_path, edge_key="role")


def pss_interaction_network_to_networkx(edge_path, node_path):
    ''' Load the PSS interaction network export to a networkx directed multigraph,
    including node attributes.

    Entities are nodes, and edges are entity -> entity influences through reactions
    (`interaction`: positive-influence, negative-influence or unknown-influence).
    Edges are keyed by `reaction_id`, as several reactions can link the same node pair.
    Mutual influences (`directed` == False, e.g. binding partners) are already listed in
    both directions.

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge file (pss-interaction-network-edges-*.tsv)

    node_path : str or pathlib.Path
        Path to the node file (pss-interaction-network-nodes-*.tsv)
    '''
    return _pss_export_to_networkx(edge_path, node_path, edge_key="reaction_id")


def pss_gene_network_to_networkx(edge_path, node_path):
    ''' Load a PSS gene network export (one species) to a networkx directed multigraph,
    including node attributes.

    As the interaction network, but with functional clusters expanded into their genes
    of one species (nodes with `node_type` == "gene").
    Edges are keyed by `reaction_id`.

    A gene can be in several functional clusters, so `entity`, `short_name`, `pathway` and
    `functional_cluster_id` are lists for every node (one entry per cluster, in the same order;
    a single entry for nodes that aren't genes). `components` of complexes are entity ids,
    i.e. match them against `entity`, not against the gene node ids.

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge file (pss-gene-network-<species>-edges-*.tsv)

    node_path : str or pathlib.Path
        Path to the node file (pss-gene-network-<species>-nodes-*.tsv)
    '''
    return _pss_export_to_networkx(edge_path, node_path, edge_key="reaction_id",
                                  node_list_columns=_PSS_GENE_CLUSTER_COLUMNS)


def ckn_to_networkx(
        edge_path=None,
        node_path=None,
        add_reciprocal_edges=True,
        directed=False,
        create_using=nx.DiGraph
    ):
    ''' Load CKN to a networkx directed graph, including node attributes.

    Downloads the CKN files from skm.nib.si if they don't exist yet.

    Parameters
    ----------
    edge_path : str or pathlib.Path
        Path to the edge list file (tab-separated, optionally gzipped);
        if the file does not exist, it is downloaded from skm.nib.si (gzipped, with a
        ``.tsv.gz`` suffix).
    node_path : str or pathlib.Path
        Path to the node annotation file; downloaded as for `edge_path`.
    add_reciprocal_edges : bool
        Add the reverse of undirected edges (``isDirected`` == 0), so directed path
        searches can use them in both directions (default True). Not meant to be used
        together with `directed`. For example A -> B (undirected) becomes A -> B and B -> A.
    directed : bool
        Remove undirected edges, and the nodes left without edges (default False).
    create_using : networkx graph class
        Graph class to create (default ``networkx.DiGraph``; ``networkx.MultiDiGraph``
        keeps parallel edges).

    Returns
    -------
    networkx.DiGraph
        CKN (or the `create_using` type). Node attributes ``GMM``, ``synonyms`` and
        ``tissue`` are lists (or None); edges have ``effect``, ``type``, ``rank``,
        ``species``, ``isDirected``, ``isTFregulation`` and ``interactionSources``.
    '''
    edge_path = Path(edge_path)
    node_path = Path(node_path)

    edge_compressed = False
    if not edge_path.exists() or (edge_path.suffix == '.gz'):
        edge_compressed = True
        if edge_path.suffix != '.gz':
            edge_path = edge_path.with_suffix(".tsv.gz")

    if not edge_path.exists():
        print(f"Attempting to download the edge list to {edge_path}.", end=" ")
        url = CKN_EDGE_URL
        urlretrieve(url, edge_path)
        print("Success.")

    node_compressed = False
    if not node_path.exists() or (node_path.suffix == '.gz'):
        node_compressed = True
        if node_path.suffix != '.gz':
            node_path = node_path.with_suffix(".tsv.gz")

    if not node_path.exists():
        print(f"Attempting to download the node annotations to {node_path}.", end=" ")
        urlretrieve(CKN_NODE_URL, node_path)
        print("Success.")

    if edge_compressed:
        open_function = gzip.open
        mode  = "tr"
    else:
        open_function = open
        mode = 'rb'

    with open_function(edge_path, mode) as handle:
        handle.readline()
        g = nx.read_edgelist(handle,
                    delimiter="\t",
                    create_using=create_using,
                    data=[
                        ('effect', str),
                        ('type', str),
                        ('rank', int),
                        ('species', str),
                        ('isDirected', int),
                        ('isTFregulation', int),
                        ('interactionSources', str)
                    ])

    if node_compressed:
        node_df = pd.read_csv(node_path, na_values=[''], keep_default_na=False, sep="\t", compression="gzip")
    else:
        node_df = pd.read_csv(node_path, na_values=[''], keep_default_na=False, sep="\t")

    node_df.set_index("node_ID", inplace=True)

    clean_list = lambda x, delim: [y.strip() for y in x.split(delim)] if not pd.isna(x) else None
    for attr, delim in [("GMM", "|"), ("synonyms", "|"), ("tissue", ",")]:
        node_df[attr] = node_df[attr].apply(clean_list, delim=delim)

    nx.set_node_attributes(g, node_df.to_dict('index'))

    if add_reciprocal_edges:
        edges_to_add = []
        for u, v, data in g.edges(data=True):
            if (data["isDirected"] == 0) and ( not g.has_edge(v, u) ):
                edges_to_add.append((v, u, data))
        _ = g.add_edges_from(edges_to_add)

    if directed:
        to_remove = [(u,v) for u, v, d in g.edges(data=True,) if d["isDirected"]==0]
        g.remove_edges_from(to_remove)

        # remove isolates resulting from filtering
        isolates = list(nx.isolates(g))
        g.remove_nodes_from(isolates)

    return g
