import networkx as nx

from skm_tools.paths import get_paths, path_edges, path_subgraph


def _diamond():
    """A -> B -> D, A -> C -> D, D -> E, X -> A"""
    return nx.DiGraph([("A", "B"), ("B", "D"), ("A", "C"), ("C", "D"), ("D", "E"), ("X", "A")])


def test_all_shortest_paths():
    assert sorted(get_paths(_diamond(), "A", "D")) == [["A", "B", "D"], ["A", "C", "D"]]


def test_single_shortest_path():
    paths = get_paths(_diamond(), "A", "D", all_shortest=False)
    assert len(paths) == 1
    assert paths[0][0] == "A" and paths[0][-1] == "D" and len(paths[0]) == 3


def test_multiple_sources_and_targets():
    paths = get_paths(_diamond(), ["A", "B"], ["D", "E"])
    assert ["B", "D"] in paths
    assert ["B", "D", "E"] in paths
    assert ["A", "B", "D", "E"] in paths


def test_shortest_overall_keeps_shortest_per_target():
    paths = get_paths(_diamond(), ["A", "B"], "D", shortest_overall=True)
    assert paths == [["B", "D"]]


def test_missing_nodes_and_no_path_are_skipped():
    g = _diamond()
    assert get_paths(g, "not-a-node", "D") == []
    assert get_paths(g, "E", "A") == []  # no directed path
    assert get_paths(g, "A", "A") == []


def test_undirected_search():
    assert get_paths(_diamond(), "E", "A") == []
    assert sorted(get_paths(_diamond(), "E", "A", directed=False)) == [["E", "D", "B", "A"], ["E", "D", "C", "A"]]


def test_path_edges_digraph():
    assert path_edges([["A", "B", "D"], ["A", "B"]], _diamond()) == [("A", "B"), ("B", "D")]


def test_path_edges_undirected_search_uses_existing_direction():
    assert path_edges([["E", "D"]], _diamond()) == [("D", "E")]


def test_path_edges_multigraph_includes_parallel_edges():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    assert path_edges([["A", "B", "C"]], g) == [("A", "B", "rx1"), ("A", "B", "rx2"), ("B", "C", "rx3")]


def test_path_subgraph():
    g = _diamond()
    g["A"]["B"]["weight"] = 3
    h = path_subgraph([["A", "B", "D"]], g)
    assert set(h.edges()) == {("A", "B"), ("B", "D")}
    assert h["A"]["B"]["weight"] == 3
    h["A"]["B"]["weight"] = 0
    assert g["A"]["B"]["weight"] == 3  # a copy
