'''Translate Arabidopsis (ath) networks to other species with SKM gene translation files.'''
from pathlib import Path

import networkx as nx
import pandas as pd

from .skm_download_urls import GENE_TRANSLATION_FILE, GENE_TRANSLATION_URL
from .utils import download_if_missing


def load_translation_file(species_code, translation_path=None, data_dir="."):
    """Load a SKM gene translation file, downloading it from skm.nib.si if missing.

    Parameters
    ----------
    species_code : str
        SKM species code of the translation (e.g. ``"stu"`` for potato, ``"parm"`` for apricot).
    translation_path : str or pathlib.Path, optional
        Local file (default ``<data_dir>/translation_ath_to_<species_code>.tsv.gz``);
        downloaded if it doesn't exist (so it can be reused offline), gzipped if the name
        ends in ``.gz``.
    data_dir : str or pathlib.Path
        Folder for the default file name (default: the current folder).

    Returns
    -------
    pandas.DataFrame
        The translation table: the Arabidopsis genes (``ath_source``), the genes of the
        species (a column named after the species, e.g. ``potato``), and the evidence. All
        values are strings (empty: None).
    """
    if translation_path is None:
        translation_path = Path(data_dir) / GENE_TRANSLATION_FILE.format(species_code)
    translation_path = download_if_missing(translation_path, GENE_TRANSLATION_URL.format(species_code))

    df = pd.read_csv(translation_path, sep='\t', dtype=str, keep_default_na=False, na_values=[""])
    return df.astype(object).where(df.notna(), None)


def integrate_translation_ckn(g,
                              translation_df,
                              t_target_col,
                              g_source_attribute='TAIR',
                              t_source_col='ath_source',
                              method='full_replacement',
                              edges_within_translation=True,
                              edges_across_translation=True,
                              keep_unmapped=True):
    """Translate the Arabidopsis genes of CKN to another species' gene identifiers.

    Each ``ath`` node (other than complexes) is replaced by one node per translation of
    its gene; translations can be many-to-many, and all are kept. Translated node ids are
    ``{original}_{translation}``, so one original can have several translations. Edges of
    the original node are copied to every new node. Nodes of other species are copied as is.

    Translation edges (see `edges_within_translation`, `edges_across_translation`) have
    ``type`` and ``interaction`` ``"homology"`` and ``directed`` False; in a directed graph,
    they are added in both directions.

    Parameters
    ----------
    g : networkx.Graph
        CKN with Arabidopsis gene identifiers (``species`` == ``"ath"``).
    translation_df : pandas.DataFrame
        Translation table, e.g. from :func:`load_translation_file`, with a source
        (Arabidopsis) and a target identifier column.
    t_target_col : str
        Column of `translation_df` with the target species' identifiers (e.g. ``"potato"``).
    g_source_attribute : str
        Node attribute of `g` with the Arabidopsis identifier (default ``"TAIR"``).
    t_source_col : str
        Column of `translation_df` with the Arabidopsis identifiers (default ``"ath_source"``).
    method : {"full_replacement"}
        How to add the translations. Only ``"full_replacement"`` (replace the original nodes)
        is implemented; ``"pendant_nodes"`` (add translations as extra nodes linked to the
        originals) is planned.
    edges_within_translation : bool
        Link the translations of the same original node to each other
        (``translation_relation`` == ``"same_translation_source"``).
    edges_across_translation : bool
        Link the nodes of different originals that translate to the same gene
        (``translation_relation`` == ``"same_translation_target"``).
    keep_unmapped : bool
        Keep ``ath`` nodes without a translation (with ``translated`` False), or drop them
        and their edges.

    Returns
    -------
    networkx.Graph
        A new graph: a (Multi)DiGraph if `g` is directed, a (Multi)Graph otherwise; a
        multigraph if either kind of translation edge is added. Translated nodes get the
        attributes ``translation``, ``translated_from`` and ``translated``.

    Raises
    ------
    NotImplementedError
        For a `method` other than ``"full_replacement"``.
    """

    if method != 'full_replacement':
        raise NotImplementedError(f"Method {method} not implemented yet. Only 'full_replacement' is currently supported.")

    # Create a mapping from the source attribute to the target attribute
    # each source identifier can have multiple target identifiers, so we need to handle that
    mapping = translation_df.groupby(t_source_col)[t_target_col].apply(list).to_dict()


    # only want to map nodes from "ath": "species" attribute is "ath"
    # if species is "ath" but there is no mapping, we can choose to either keep the original node or drop it (keep_unmapped)
    # if keep, add an annotation that is is unmapped

    # non "ath" nodes, are copied over as is

    # first add all the nodes (and track ids), then add all the edges

    if g.is_directed():
        if (edges_across_translation or edges_within_translation):
            g_translated = nx.MultiDiGraph()  # use MultiDiGraph to allow for multiple edges between the same nodes (e.g. from original and from translation)
        else:
            g_translated = nx.DiGraph()
    else:
        if (edges_across_translation or edges_within_translation):
            g_translated = nx.MultiGraph()  # use MultiGraph to allow for multiple edges between the same nodes (e.g. from original and from translation)
        else:
            g_translated = nx.Graph()

    id_tracking = {}  # mapping from original node id to list of new node ids (after translation)
    translation_to_new_nodes = {}

    for node, data in g.nodes(data=True):
        source_id = data.get(g_source_attribute)
        species = data.get('species')
        node_type = data.get('node_type')

        if (species == 'ath') and (node_type != 'Complex'):
            if source_id in mapping:
                target_ids = mapping[source_id]
                new_node_ids = []
                for target_id in target_ids:
                    new_node_id = f"{node}_{target_id}"
                    g_translated.add_node(new_node_id, **data, translation=target_id, translated_from=node, translated=True)
                    new_node_ids.append(new_node_id)
                    translation_to_new_nodes.setdefault(target_id, []).append(new_node_id)
                id_tracking[node] = new_node_ids
            else:
                if keep_unmapped:
                    g_translated.add_node(node, **data, translation=None, translated_from=node, translated=False)
                    id_tracking[node] = [node]
        else:
            g_translated.add_node(node, **data)
            id_tracking[node] = [node]

    # add edges
    for u, v, data in g.edges(data=True):
        # nodes dropped (keep_unmapped=False) have no new ids: their edges are dropped too
        u_new_ids = id_tracking.get(u, [])
        v_new_ids = id_tracking.get(v, [])

        for u_new in u_new_ids:
            for v_new in v_new_ids:
                g_translated.add_edge(u_new, v_new, **data)


    # supplemental edges
    def add_translation_edges(groups, relation):
        for new_nodes in groups:
            for i in range(len(new_nodes)):
                for j in range(i + 1, len(new_nodes)):
                    pairs = [(new_nodes[i], new_nodes[j])]
                    if g_translated.is_directed():
                        pairs.append((new_nodes[j], new_nodes[i]))
                    for u, v in pairs:
                        g_translated.add_edge(u, v, translation_relation=relation, type="homology",
                                              interaction="homology", directed=False)

    if edges_within_translation:
        # translations of the same original node
        add_translation_edges(id_tracking.values(), "same_translation_source")

    if edges_across_translation:
        # nodes of different originals that translate to the same gene
        add_translation_edges(translation_to_new_nodes.values(), "same_translation_target")

    return g_translated
