'''Load from file or server to nx object'''

from urllib.request import urlretrieve
import gzip
from pathlib import Path
import networkx as nx
import pandas as pd
from .skm_download_urls import *

def pss_to_networkx(edge_path=None, node_path=None):
    ''' Load PSS from the rxn bipartite projection SIF format to a
    networkx directed multigraph format, including node attributes

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge list file,
        if file does not exist, download from skm.nib.si

    node_path : str or pathlib.Path
        Path to the node annotation file,
        if file does not exist, download from skm.nib.si

    '''
    edge_path = Path(edge_path)
    node_path = Path(node_path)

    if not edge_path.exists():
        print(f"Attempting to download the edge list to {edge_path}.", end=" ")
        urlretrieve(PSS_RXN_EDGE_URL, edge_path)
        print("Success.")

    if not node_path.exists():
        print(f"Attempting to download the node annotations to {node_path}.", end=" ")
        urlretrieve(PSS_RXN_NODE_URL, node_path)
        print("Success.")

    with open(edge_path, "rb") as handle:
        handle.readline()
        g = nx.read_edgelist(handle,
                    delimiter="\t", create_using=nx.MultiDiGraph,
                        data=list({
                            'interaction_type':str,
                            'directed': str,
                            'reaction_type': str,
                            'reaction_effect': str,
                            'reaction_id': str,
                            'source_edge_type': str,
                            'source_location': str,
                            'source_form': str,
                            'target_edge_type': str,
                            'target_location': str,
                            'target_form': str,
                        }.items()))

    node_df = pd.read_csv(node_path, sep="\t")
    node_df.set_index("name", inplace=True, drop=False)

    def to_list(x):
        if isinstance(x, str):
            return x.split(',')
        return None

    for c in [
        'ath_homologues',
        'osa_homologues',
        'stu_homologues',
        'sly_homologues'
        ]:
        node_df[c] = node_df[c].apply(to_list)

    nx.set_node_attributes(g, node_df.to_dict('index'))

    for _, data in g.nodes(data=True):
        if data["node_type"] in ("PlantCoding", "PlantAbstract", "PlantNonCoding"):
            data["display_label"] = data["short_name"]
        else:
            data["display_label"] = data["name"]

    return g



def pss_model_to_networkx(edge_path=None, node_path=None):
    ''' Load PSS from the main SIF format to a
    networkx directed multigraph format, including node attributes

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge list file,
        if file does not exist, download from skm.nib.si

    node_path : str or pathlib.Path
        Path to the node annotation file,
        if file does not exist, download from skm.nib.si

    '''
    edge_path = Path(edge_path)
    node_path = Path(node_path)

    if not edge_path.exists():
        print(f"Attempting to download the edge list to {edge_path}.", end=" ")
        urlretrieve(PSS_EDGE_URL, edge_path)
        print("Success.")

    if not node_path.exists():
        print(f"Attempting to download the node annotations to {node_path}.", end=" ")
        urlretrieve(PSS_NODE_URL, node_path)
        print("Success.")

    with open(edge_path, "rb") as handle:
        handle.readline()
        g = nx.read_edgelist(handle,
                    delimiter="\t", create_using=nx.MultiDiGraph,
                        data=list({
                            'interaction_type':str,
                            'directed': str,
                            'reaction_type': str,
                            'reaction_effect': str,
                            'reaction_id': str,
                            'source_edge_type': str,
                            'source_location': str,
                            'source_form': str,
                            'target_edge_type': str,
                            'target_location': str,
                            'target_form': str,
                        }.items()))

    node_df = pd.read_csv(node_path, sep="\t")
    node_df.set_index("name", inplace=True, drop=False)

    def to_list(x):
        if isinstance(x, str):
            return x.split(',')
        return None

    for c in [
        'ath_homologues',
        'osa_homologues',
        'stu_homologues',
        'sly_homologues'
        ]:
        node_df[c] = node_df[c].apply(to_list)

    nx.set_node_attributes(g, node_df.to_dict('index'))

    for _, data in g.nodes(data=True):
        if data["node_type"] in ("PlantCoding", "PlantAbstract", "PlantNonCoding"):
            data["display_label"] = data["short_name"]
        else:
            data["display_label"] = data["name"]

    return g

























def ckn_to_networkx(
        edge_path=None,
        node_path=None,
        add_reciprocal_edges=True,
        directed=False,
        create_using=nx.DiGraph
    ):
    ''' Load CKN from the compressed edge list format to a
    networkx directed multigraph format, including node attributes

    Parameters
    ----------

    edge_path : str or pathlib.Path
        Path to the edge list file,
        if file does not exist, download from skm.nib.si

    node_path : str or pathlib.Path
        Path to the node annotation file,
        if file does not exist, download from skm.nib.si

    expanded : bool
        Whether to use the expanded CKN. If False, use CKN
        collapsed to a single edge between any pair of nodes.
        Ignored if edge_path already exists

    add_reciprocal_edges : bool
        Whether to add reciprocal edges of undirected edges
        (necessary for doing directed path analysis).
        e.g.
            A-->B (undirected)
        becomes
            A-->B (undirected)
            B-->A (undirected)

    directed : bool
        Whether to remove undirected edges from CKN
        (in contrast to add_reciprocal_edges, and not meant to be used together)
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
