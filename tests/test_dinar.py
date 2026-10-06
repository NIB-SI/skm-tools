"""Tests for skm_tools.dinar"""

import networkx as nx
import pandas as pd
import pytest

from skm_tools.dinar import to_dinar, write_dinar
from skm_tools.pss import pss_gene_network_to_networkx


def _graph():
    g = nx.MultiDiGraph()
    g.add_node("AT2G38470", display_label="WRKY33", description="WRKY DNA-binding protein 33",
               mapman=["14.5.7.5.1_RNA biosynthesis.WRKY", "26.11.3.2.1_External stimuli"])
    g.add_node("camalexin", display_label=None, description=None, mapman=None)
    g.add_edge("AT2G38470", "camalexin", key="rx1", interaction="positive-influence")
    g.add_edge("AT2G38470", "camalexin", key="rx2", interaction="positive-influence")
    g.add_edge("AT2G38470", "camalexin", key="rx3", interaction="negative-influence")
    g.add_edge("camalexin", "AT2G38470", key="rx4")
    return g


TEXT = ["shortDescription", "shortName", "MapManBin"]


def test_node_columns():
    nodes, _ = to_dinar(_graph())
    assert list(nodes.columns) == [
        "geneID", "shortDescription", "shortName", "MapManBin", "clusterID", "clusterName",
        "x", "y", "clusterSimplifiedNodeDegree", "expressed",
    ]
    assert nodes.set_index("geneID").loc["AT2G38470", TEXT].tolist() == [
        "WRKY DNA-binding protein 33", "WRKY33",
        "14.5.7.5.1_RNA biosynthesis.WRKY | 26.11.3.2.1_External stimuli",
    ]


def test_missing_values_are_dashes_and_name_falls_back_to_id():
    nodes, _ = to_dinar(_graph())
    assert nodes.set_index("geneID").loc["camalexin", TEXT].tolist() == ["-", "camalexin", "-"]


def test_edges_one_row_per_interaction():
    _, edges = to_dinar(_graph())
    assert list(edges.columns) == [
        "geneID1", "geneID2", "reactionType", "clusterID_geneID1", "clusterID_geneID2",
        "clusterSimplifiedNodeDegree_geneID1", "clusterSimplifiedNodeDegree_geneID2", "exists",
    ]
    assert sorted(map(tuple, edges[["geneID1", "geneID2", "reactionType"]].values.tolist())) == [
        ("AT2G38470", "camalexin", "negative-influence"),
        ("AT2G38470", "camalexin", "positive-influence"),
        ("camalexin", "AT2G38470", "-"),
    ]


def test_write_dinar(tmp_path):
    write_dinar(_graph(), tmp_path / "nodes.txt", tmp_path / "edges.txt")
    nodes = pd.read_csv(tmp_path / "nodes.txt", sep="\t", dtype=str, keep_default_na=False)
    edges = pd.read_csv(tmp_path / "edges.txt", sep="\t", dtype=str, keep_default_na=False)
    assert len(nodes) == 2 and len(edges) == 3
    assert (nodes.values != "").all() and (edges.values != "").all()


def test_gene_network(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
    g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
    nodes, edges = to_dinar(g)
    assert set(nodes["geneID"]) == set(g.nodes())
    assert set(edges["reactionType"]) <= {"positive-influence", "negative-influence", "unknown-influence"}
    assert set(edges["geneID1"]) | set(edges["geneID2"]) <= set(nodes["geneID"])


def _clustered_graph():
    """two triangles (pathways A and B) joined by one edge, plus an unclustered node"""
    g = nx.DiGraph([("a1", "a2"), ("a2", "a3"), ("a3", "a1"), ("a1", "a2"),
                    ("b1", "b2"), ("b2", "b3"), ("b3", "b1"), ("b4", "b1"),
                    ("a1", "b1"), ("x", "a1"), ("a1", "a1")])
    for n in g:
        g.nodes[n]["pathway"] = {"a": ["A", "B"], "b": "B"}.get(n[0])
    return g


def test_one_cluster_by_default():
    nodes, edges = to_dinar(_clustered_graph())
    assert set(nodes["clusterID"]) == {1}
    assert (edges["exists"] == 1).all() and (nodes["expressed"] == 1).all()


def test_clusters_from_attribute():
    nodes, _ = to_dinar(_clustered_graph(), clusters="pathway")
    n = nodes.set_index("geneID")
    # B (4 nodes) is the largest cluster: 1; lists use the first value; no pathway: 0
    assert n.loc["b1", ["clusterID", "clusterName"]].tolist() == [1, "B"]
    assert n.loc["a1", ["clusterID", "clusterName"]].tolist() == [2, "A"]
    assert n.loc["x", ["clusterID", "clusterName"]].tolist() == [0, "-"]


def test_cluster_degree_ignores_direction_loops_and_other_clusters():
    nodes, edges = to_dinar(_clustered_graph(), clusters="pathway")
    n = nodes.set_index("geneID")["clusterSimplifiedNodeDegree"]
    assert n["a1"] == 2   # a2, a3; not a1 itself, b1 or x
    assert n["b1"] == 3   # b2, b3, b4
    e = edges.set_index(["geneID1", "geneID2"]).loc[("a1", "b1")]
    assert e[["clusterID_geneID1", "clusterID_geneID2"]].tolist() == [2, 1]
    assert e[["clusterSimplifiedNodeDegree_geneID1", "clusterSimplifiedNodeDegree_geneID2"]].tolist() == [2, 3]


def test_clusters_from_mapping():
    clusters = {"a1": "c", "a2": "c", "b1": "d"}
    nodes, _ = to_dinar(_clustered_graph(), clusters=clusters)
    n = nodes.set_index("geneID")["clusterID"]
    assert n["a1"] == n["a2"] == 1 and n["b1"] == 2 and n["b2"] == 0


def test_positions_from_mapping_and_pos_attribute():
    g = _clustered_graph()
    pos = {n: (i, -i) for i, n in enumerate(g)}
    nodes, _ = to_dinar(g, positions=pos)
    assert nodes.set_index("geneID").loc["a2", ["x", "y"]].tolist() == [1.0, -1.0]

    nx.set_node_attributes(g, pos, "pos")
    assert to_dinar(g)[0].equals(nodes)


def test_default_layout_is_reproducible_and_groups_clusters():
    g = _clustered_graph()
    nodes, _ = to_dinar(g, clusters="pathway")
    assert nodes.equals(to_dinar(g, clusters="pathway")[0])
    # each cluster is laid out in its own 2 x 2 box
    for _, c in nodes.groupby("clusterID"):
        assert c["x"].max() - c["x"].min() <= 2 and c["y"].max() - c["y"].min() <= 2


def test_missing_positions_raise():
    with pytest.raises(ValueError, match="No position"):
        to_dinar(_clustered_graph(), positions={"a1": (0, 0)})
