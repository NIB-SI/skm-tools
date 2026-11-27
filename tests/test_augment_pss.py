import logging
from unittest.mock import patch

import networkx as nx
import pytest

from skm_tools import augment_pss as ap


def _pss():
    '''PSS: A -> B; C without edges to A or B'''
    g = nx.MultiDiGraph(pss_network="gene_network", species="ath")
    for n in "ABC":
        g.add_node(n, node_type="PlantCoding", species="ath")
    g.add_edge("A", "B", key="rx1", reaction_id="rx1", interaction="positive-influence")
    g.add_edge("B", "C", key="rx2", reaction_id="rx2", interaction="positive-influence")
    return g


def _ckn():
    '''CKN around the PSS genes:
    B -> A (PSS has A -> B: skipped), A -> C (new direct, rank 1; also C -> A, rank 0),
    A -> X, X -> B (X: an intermediate of A and B, but they have a PSS edge),
    A -> Y, Y -> C (Y: an intermediate of A and C), X -> Y (between added nodes),
    Z -> W (not touching PSS)'''
    g = nx.DiGraph()
    for n in "ABCXYZW":
        g.add_node(n, node_type="PlantCoding", display_label=f"ckn {n}")
    edges = [("B", "A", 0), ("A", "C", 1), ("C", "A", 0), ("A", "X", 2), ("X", "B", 1),
             ("A", "Y", 1), ("Y", "C", 3), ("X", "Y", 0), ("Z", "W", 0)]
    for u, v, rank in edges:
        g.add_edge(u, v, rank=rank, interaction="unknown-influence", interactionSources=["x"])
    return g


@pytest.fixture
def augmented():
    return ap.build_augmented_pss(_pss(), _ckn())


def _edges(g, origin):
    return {(u, v) for u, v, d in g.edges(data=True) if d["origin"] == origin}


def test_build_adds_the_ckn_edges_touching_pss(augmented):
    assert _edges(augmented, "PSS") == {("A", "B"), ("B", "C")}
    assert _edges(augmented, "CKN") == {("A", "C"), ("C", "A"), ("A", "X"), ("X", "B"), ("A", "Y"),
                                        ("Y", "C"), ("X", "Y")}
    assert {n: d["origin"] for n, d in augmented.nodes(data=True)} == {
        "A": "PSS", "B": "PSS", "C": "PSS", "X": "CKN", "Y": "CKN"}
    # PSS keeps its keys, attributes and graph attributes; CKN attributes are copied, typed
    assert augmented["A"]["B"]["rx1"]["reaction_id"] == "rx1"
    assert augmented.graph == {"pss_network": "gene_network", "species": "ath"}
    assert augmented.nodes["A"]["species"] == "ath" and "display_label" not in augmented.nodes["A"]
    assert augmented.nodes["X"]["display_label"] == "ckn X"
    ckn_edge = next(iter(augmented["A"]["C"].values()))
    assert ckn_edge["rank"] == 1 and ckn_edge["interactionSources"] == ["x"]


def test_build_steps(augmented):
    steps = {(u, v): d["step"] for u, v, d in augmented.edges(data=True)}
    assert steps[("A", "B")] == "input"
    assert steps[("A", "X")] == "step1"
    assert steps[("X", "Y")] == "step1+step2"
    assert augmented.nodes["X"]["step"] == "step1" and augmented.nodes["A"]["step"] == "input"


def test_build_without_edges_between_added_nodes():
    g = ap.build_augmented_pss(_pss(), _ckn(), between_added=False)
    assert ("X", "Y") not in _edges(g, "CKN")


def test_build_leaves_the_inputs_unchanged():
    pss, ckn = _pss(), _ckn()
    ap.build_augmented_pss(pss, ckn)
    assert "origin" not in pss.nodes["A"] and "origin" not in ckn["A"]["C"]


def test_build_from_a_digraph():
    pss = nx.DiGraph([("A", "B")])
    g = ap.build_augmented_pss(pss, _ckn())
    assert _edges(g, "PSS") == {("A", "B")}


def test_build_warns_without_shared_nodes(caplog):
    with caplog.at_level(logging.WARNING, logger="skm_tools"):
        ap.build_augmented_pss(nx.DiGraph([("MYC2", "JAZ")]), _ckn())
    assert "No PSS node is in CKN" in caplog.text


def test_build_rejects_the_reaction_graph():
    with pytest.raises(ValueError, match="reaction graph"):
        ap.build_augmented_pss(nx.MultiDiGraph(pss_network="reaction_graph"), _ckn())


def test_check_unaugmented_nodes(augmented):
    pss = _pss()
    pss.add_node("D")
    assert ap.check_unaugmented_nodes(pss, augmented) == {"D"}


def test_summary_and_diff(augmented):
    assert ap.get_augmentation_summary(_pss(), augmented) == {
        "pss_nodes": 3, "pss_edges": 2, "augmented_nodes": 5, "augmented_edges": 9,
        "ckn_nodes_added": 2, "ckn_edges_added": 7}
    diff = ap.augmentation_diff(augmented)
    assert diff.kind.value_counts().to_dict() == {"edge": 7, "node": 2}
    assert set(diff[diff.kind == "node"].source) == {"X", "Y"}
    assert ap.augmentation_diff(ap.build_augmented_pss(_pss(), nx.DiGraph())).empty


def test_find_new_direct_edges(augmented):
    augmented.add_edge("B", "B", origin="CKN", rank=0)  # a self-loop is not a pair
    df = ap.find_new_direct_edges(augmented, _pss())
    assert df.to_dict("records") == [{"u": "A", "v": "C", "rank": 0, "direction": "both"}]
    assert ap.find_new_direct_edges(augmented, _pss(), max_rank=-1).empty


def test_find_new_intermediate_pairs(augmented):
    df = ap.find_new_intermediate_pairs(augmented, _pss())
    # through X: A and B have a PSS edge; through Y: A and C
    assert df.to_dict("records") == [
        {"u": "A", "v": "C", "intermediate": "Y", "rank_ux": 1, "rank_xv": 3, "path_rank": 3}]
    assert ap.find_new_intermediate_pairs(augmented, _pss(), max_rank=2).empty


def test_extract_pair_subnetwork(augmented):
    sub = ap.extract_pair_subnetwork(augmented, _pss(), "A", "C", intermediate="Y")
    assert dict(sub.nodes(data="role")) == {"A": "focus", "C": "focus", "Y": "intermediate",
                                           "B": "context"}
    status = {(u, v): d["edge_status"] for u, v, d in sub.edges(data=True)}
    assert status[("A", "C")] == status[("A", "Y")] == "new_from_ckn"
    assert status[("A", "B")] == "pss_original"
    assert "role" not in augmented.nodes["A"]  # a copy
    assert set(ap.extract_pair_subnetwork(augmented, _pss(), "A", "C", max_neighbours=0)) == {"A", "C"}


def test_pick_random_node(augmented):
    assert ap.pick_random_node(augmented, origin="CKN", seed=1) in {"X", "Y"}
    assert ap.pick_random_node(augmented, seed=3) == ap.pick_random_node(augmented, seed=3)
    with pytest.raises(ValueError):
        ap.pick_random_node(augmented, origin="other")


def test_extract_node_neighbourhood(augmented):
    sub = ap.extract_node_neighbourhood(augmented, "X")
    assert set(sub) == {"X", "A", "B", "Y"} and sub.nodes["X"]["role"] == "focus"
    capped = ap.extract_node_neighbourhood(augmented, "A", radius=2, max_neighbours=1)
    assert "A" in capped and len(capped) <= 3
    with pytest.raises(ValueError):
        ap.extract_node_neighbourhood(augmented, "nope")


def test_get_node_neighbourhood_augmentation(augmented):
    stats = ap.get_node_neighbourhood_augmentation(_pss(), augmented, "A")
    assert stats == {"node": "A", "radius": 1, "pss_neighbours": 1, "augmented_neighbours": 4,
                     "ckn_added_neighbours": 2, "pct_new": 50.0}


def test_show_augmented(augmented):
    from skm_tools import cytoscape_utils as cu
    sub = ap.extract_pair_subnetwork(augmented, _pss(), "A", "C", intermediate="Y")
    with patch.object(cu, "p4c") as p4c, patch.object(cu, "load_network", return_value=7) as load, \
            patch.object(cu, "apply_builtin_style", return_value="SKM"), \
            patch.object(cu, "copy_style") as copy_style:
        assert ap.show_augmented(sub, "pair") == 7
    load.assert_called_once_with(sub, "pair", collection=None)
    copy_style.assert_called_once_with("SKM", "SKM-augmented", networks=[7], overwrite=True)
    assert p4c.set_edge_color_mapping.call_args.args == ("origin", ["PSS", "CKN"], ["#2171B5", "#E6550D"])
    assert p4c.set_node_border_width_mapping.call_args.args[0] == "role"
    assert p4c.set_edge_line_width_mapping.call_args.args[0] == "edge_status"


def test_show_augmented_whole_network_borders_by_origin(augmented):
    from skm_tools import cytoscape_utils as cu
    with patch.object(cu, "p4c") as p4c, patch.object(cu, "load_network", return_value=7), \
            patch.object(cu, "apply_builtin_style", return_value="SKM"), patch.object(cu, "copy_style"):
        ap.show_augmented(augmented, "all")
    assert p4c.set_node_border_width_mapping.call_args.args[0] == "origin"
    p4c.set_edge_line_width_mapping.assert_not_called()
