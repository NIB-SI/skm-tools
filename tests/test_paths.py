import networkx as nx
import pytest

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


def test_missing_nodes_and_no_path_are_skipped(caplog):
    g = _diamond()
    assert get_paths(g, "not-a-node", "D") == []
    assert "1 of 1 sources not in the graph: not-a-node" in caplog.text
    assert get_paths(g, "E", "A") == []  # no directed path
    assert get_paths(g, "A", "A") == []


def test_undirected_search():
    assert get_paths(_diamond(), "E", "A") == []
    assert sorted(get_paths(_diamond(), "E", "A", directed=False)) == [["E", "D", "B", "A"], ["E", "D", "C", "A"]]


def test_path_edges_digraph():
    assert path_edges(_diamond(), [["A", "B", "D"], ["A", "B"]]) == [("A", "B"), ("B", "D")]


def test_path_edges_undirected_search_uses_existing_direction():
    assert path_edges(_diamond(), [["E", "D"]]) == [("D", "E")]


def test_path_edges_multigraph_includes_parallel_edges():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    assert path_edges(g, [["A", "B", "C"]]) == [("A", "B", "rx1"), ("A", "B", "rx2"), ("B", "C", "rx3")]


def test_path_subgraph():
    g = _diamond()
    g["A"]["B"]["weight"] = 3
    h = path_subgraph(g, [["A", "B", "D"]])
    assert set(h.edges()) == {("A", "B"), ("B", "D")}
    assert h["A"]["B"]["weight"] == 3
    h["A"]["B"]["weight"] = 0
    assert g["A"]["B"]["weight"] == 3  # a copy


def test_undirected_search_on_multidigraph():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("C", "B", key="rx3")
    assert get_paths(g, "A", "C") == []
    assert get_paths(g, "A", "C", directed=False) == [["A", "B", "C"]]


def test_undirected_graph():
    g = nx.Graph([("A", "B"), ("B", "C"), ("A", "D"), ("D", "C")])
    assert sorted(get_paths(g, "C", "A")) == [["C", "B", "A"], ["C", "D", "A"]]


def test_single_shortest_path_several_sources_shortest_overall():
    paths = get_paths(_diamond(), ["X", "A", "B"], ["D", "E"], all_shortest=False, shortest_overall=True)
    assert paths == [["B", "D"], ["B", "D", "E"]]


def test_order_by_target_then_source():
    g = nx.DiGraph([("S1", "T1"), ("S2", "T1"), ("S1", "T2")])
    assert get_paths(g, ["S2", "S1"], ["T1", "T2"]) == [["S2", "T1"], ["S1", "T1"], ["S1", "T2"]]


def test_max_paths(caplog):
    # a 4 x 4 grid has 20 shortest paths from corner to corner
    g = nx.DiGraph(nx.grid_2d_graph(4, 4))
    assert len(get_paths(g, [(0, 0)], [(3, 3)])) == 20  # tuple node ids: in a list
    assert len(get_paths(g, [(0, 0)], [(3, 3)], max_paths=5)) == 5
    assert "Stopped at max_paths=5" in caplog.text


@pytest.mark.parametrize("max_paths", [0, -1, 1.5, True])
def test_bad_max_paths(max_paths):
    with pytest.raises(ValueError):
        get_paths(_diamond(), "A", "D", max_paths=max_paths)


def test_path_subgraph_multigraph():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1", interaction="positive-influence")
    g.add_edge("A", "B", key="rx2", interaction="negative-influence")
    g.add_edge("B", "C", key="rx3")
    g.add_edge("A", "C", key="rx4")
    h = path_subgraph(g, [["A", "B", "C"]])
    assert type(h) is nx.MultiDiGraph
    assert set(h.edges(keys=True)) == {("A", "B", "rx1"), ("A", "B", "rx2"), ("B", "C", "rx3")}
