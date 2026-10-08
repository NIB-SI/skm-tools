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
        "geneID", "shortDescription", "shortName", "MapManBin", "clusterID",
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
        ("AT2G38470", "camalexin", "act_influence"),
        ("AT2G38470", "camalexin", "inh_influence"),
        ("camalexin", "AT2G38470", "unk_influence"),
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
    # as DiNAR's PSS tables: the sign (act/inh/unk, which DiNAR draws), then the reaction type
    signs = {"positive-influence": "act", "negative-influence": "inh", "unknown-influence": "unk"}
    expected = {f"{signs[d['interaction']]}_{d['reaction_type']}" for *_, d in g.edges(data=True)}
    assert set(edges["reactionType"]) == expected
    assert "act_protein activation" in expected
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
    # B (4 nodes) is the largest cluster: 1; lists use the first value; no pathway: one more
    # cluster, last
    assert n.loc["b1", "clusterID"] == 1
    assert n.loc["a1", "clusterID"] == 2
    assert n.loc["x", "clusterID"] == 3


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
    assert n["a1"] == n["a2"] == 1 and n["b1"] == 2 and n["b2"] == n["x"] == 3


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


def test_nodes_without_edges_are_cluster_0():
    g = _clustered_graph()
    g.add_node("lonely", pathway="A")
    nodes, _ = to_dinar(g, clusters="pathway")
    n = nodes.set_index("geneID")["clusterID"]
    assert n["lonely"] == 0
    assert sorted(set(n)) == [0, 1, 2, 3]


def test_nan_and_set_cluster_labels():
    g = _clustered_graph()
    clusters = {n: float("nan") for n in g} | {"a1": {"Z", "A"}, "a2": ["A"]}
    n = to_dinar(g, clusters=clusters)[0].set_index("geneID")["clusterID"]
    # NaN: no cluster (all together), not a cluster each; sets: the first in sorted order
    assert n["a1"] == n["a2"] == 1
    assert set(n.drop(["a1", "a2"])) == {2}  # no cluster: together, last


def test_text_values():
    g = nx.DiGraph([("A", "B")])
    g.nodes["A"]["description"] = "line 1\nline\t2"
    g.nodes["B"]["description"] = float("nan")
    nodes, _ = to_dinar(g)
    assert nodes["shortDescription"].tolist() == ["line 1 line 2", "-"]


def test_colliding_ids_raise():
    with pytest.raises(ValueError, match="same as text"):
        to_dinar(nx.Graph([(1, "1")]))


def test_cluster_edge_limit():
    g = nx.gnm_random_graph(200, 2001, seed=1, directed=True)
    with pytest.raises(ValueError, match="more than 2000 edges"):
        to_dinar(g, positions={n: (0, 0) for n in g})
    clusters = {n: n % 2 for n in g}
    to_dinar(g, clusters=clusters, positions={n: (0, 0) for n in g})  # two smaller clusters


def test_written_without_quoting(tmp_path):
    g = nx.DiGraph([("A", "B")])
    g.nodes["A"]["description"] = 'U2B"-LIKE'
    write_dinar(g, tmp_path / "nodes.txt", tmp_path / "edges.txt")
    assert '\tU2B"-LIKE\t' in (tmp_path / "nodes.txt").read_text()
