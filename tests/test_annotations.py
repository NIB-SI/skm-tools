import networkx as nx

from skm_tools.annotations import get_all_annotations, get_nodes_by_mapman
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
