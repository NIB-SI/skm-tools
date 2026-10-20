import networkx as nx

from skm_tools.annotations import get_all_annotations, get_nodes_by_mapman, mercator_gene_sets
from skm_tools.ckn import ckn_to_networkx
from skm_tools.pss import pss_interaction_network_to_networkx


def _graph():
    g = nx.Graph()
    g.add_node("A", mapman=["26.11.3.2.1_External stimuli response.pathogen...", "15.5_RNA biosynthesis..."])
    g.add_node("B", mapman=["26.110_not a child of 26.11"])
    g.add_node("C", mapman=["26.11_External stimuli response.pathogen"])
    g.add_node("D", mapman=None)
    return g


def test_children():
    assert sorted(get_nodes_by_mapman(_graph(), "26.11")) == ["A", "C"]


def test_without_children():
    assert get_nodes_by_mapman(_graph(), "26.11", children=False) == ["C"]


def test_several_bins_and_names():
    assert sorted(get_nodes_by_mapman(_graph(), ["15.5_RNA biosynthesis", "26.110"])) == ["A", "B"]


def test_get_all_annotations():
    assert len(get_all_annotations(_graph())) == 4


def test_get_all_annotations_single_values():
    g = nx.Graph()
    g.add_node("A", tissue="leaf")
    g.add_node("B", tissue=["root", "leaf"])
    assert get_all_annotations(g, "tissue") == {"leaf", "root"}


def test_string_bin():
    assert get_nodes_by_mapman(_graph(), "15.5") == ["A"]


def test_same_bins_work_on_pss_and_ckn(pss_interaction_network_edge_path, pss_interaction_network_node_path,
                                       ckn_edge_path, ckn_node_path):
    pss = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
    pss_hits = get_nodes_by_mapman(pss, [a.split(".")[0] for a in get_all_annotations(pss)])
    assert pss_hits
    ckn = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    bins = {a.split(".")[0] for a in get_all_annotations(ckn)}
    assert get_nodes_by_mapman(ckn, list(bins))


MERCATOR = """BINCODE\tNAME\tIDENTIFIER\tDESCRIPTION\tTYPE
'1'\t'Photosynthesis'\t''\t''\t
'1.1'\t'Photosynthesis.photophosphorylation'\t''\t''\t
'1.1.1'\t'Photosynthesis.photophosphorylation.photosystem II'\t''\t''\t
'1.1.1'\t'Photosynthesis.photophosphorylation.photosystem II'\t'pruarm.1g000100.t1.p1_1'\t'desc'\tT
'1.1.1'\t'Photosynthesis.photophosphorylation.photosystem II'\t'pruarm.1g000100.t2.p1_1'\t'desc'\tT
'1.2'\t'Photosynthesis.calvin cycle'\t'pruarm.1g000200.t1.p1_1'\t'desc'\tT
'10'\t'Redox homeostasis'\t'pruarm.2g000300.t1.p1_1'\t'desc'\tT
'35'\t'not assigned'\t''\t''\t
'35.2'\t'not assigned.not annotated'\t'pruarm.2g000400.t1.p1_1'\t'desc'\tT
"""


def _gene(identifier):
    return identifier.split(".t")[0]


def test_mercator_gene_sets(tmp_path):
    path = tmp_path / "mercator.txt"
    path.write_text(MERCATOR)
    sets = mercator_gene_sets(path, id_to_gene=_gene)
    assert sets == {
        "1_Photosynthesis": ["pruarm.1g000100", "pruarm.1g000200"],     # with the sub-bins' genes
        "1.1_Photosynthesis.photophosphorylation": ["pruarm.1g000100"],
        "1.1.1_Photosynthesis.photophosphorylation.photosystem II": ["pruarm.1g000100"],   # once
        "1.2_Photosynthesis.calvin cycle": ["pruarm.1g000200"],
        "10_Redox homeostasis": ["pruarm.2g000300"],                    # after 1.x, not after 1
    }


def test_mercator_gene_sets_filters(tmp_path):
    path = tmp_path / "mercator.txt.gz"
    import gzip
    path.write_bytes(gzip.compress(MERCATOR.encode()))
    assert list(mercator_gene_sets(path, levels=[2], id_to_gene=_gene)) == [
        "1.1_Photosynthesis.photophosphorylation", "1.2_Photosynthesis.calvin cycle"]
    assert list(mercator_gene_sets(path, min_size=2, id_to_gene=_gene)) == ["1_Photosynthesis"]
    assert list(mercator_gene_sets(path, max_size=1, levels=1, id_to_gene=_gene)) == ["10_Redox homeostasis"]
    assert "35_not assigned" in mercator_gene_sets(path, exclude_bins=[], id_to_gene=_gene)
    # without id_to_gene: the identifiers
    assert mercator_gene_sets(path, levels=[3])["1.1.1_Photosynthesis.photophosphorylation.photosystem II"] == [
        "pruarm.1g000100.t1.p1_1", "pruarm.1g000100.t2.p1_1"]
