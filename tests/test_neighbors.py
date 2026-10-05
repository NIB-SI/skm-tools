import networkx as nx
import pytest

from skm_tools.neighbors import get_neighborhood, neighborhood_nodes


def _chain():
    """U -> A -> B -> C -> D, and A -> X, X -> B"""
    return nx.DiGraph([("U", "A"), ("A", "B"), ("B", "C"), ("C", "D"), ("A", "X"), ("X", "B")])


def test_first_neighbours_both_directions():
    assert neighborhood_nodes(_chain(), "A") == {"A": 0, "U": 1, "B": 1, "X": 1}


def test_depth():
    assert neighborhood_nodes(_chain(), "A", depth=2, direction="out") == {"A": 0, "B": 1, "X": 1, "C": 2}
    assert neighborhood_nodes(_chain(), "A", depth=0) == {"A": 0}


def test_direction_in():
    assert neighborhood_nodes(_chain(), "B", direction="in") == {"B": 0, "A": 1, "X": 1}


def test_several_start_nodes_take_closest_distance():
    d = neighborhood_nodes(_chain(), ["A", "D"], depth=2, direction="out")
    assert d["A"] == 0 and d["D"] == 0 and d["C"] == 2


def test_unknown_nodes_ignored_and_bad_direction_raises():
    assert neighborhood_nodes(_chain(), ["nope"]) == {}
    with pytest.raises(ValueError):
        neighborhood_nodes(_chain(), "A", direction="sideways")


def test_undirected_graph():
    g = nx.Graph([("A", "B"), ("B", "C")])
    assert neighborhood_nodes(g, "A", depth=2, direction="out") == {"A": 0, "B": 1, "C": 2}


def test_neighbourhood_is_node_induced_copy_with_distance():
    g = _chain()
    h = get_neighborhood(g, "A")
    assert set(h) == {"U", "A", "B", "X"}
    assert h.has_edge("X", "B")  # between two first neighbours
    assert h.nodes["B"]["distance"] == 1
    assert "distance" not in g.nodes["B"]
    assert type(h) is nx.DiGraph


def test_neighbourhood_not_induced_only_followed_edges():
    h = get_neighborhood(_chain(), "A", induced=False)
    assert set(h.edges()) == {("U", "A"), ("A", "B"), ("A", "X")}


def test_neighbourhood_multigraph_keeps_parallel_edges():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    h = get_neighborhood(g, "A", induced=False)
    assert set(h.edges(keys=True)) == {("A", "B", "rx1"), ("A", "B", "rx2")}


def test_isolated_start_node_kept():
    g = _chain()
    g.add_node("lonely")
    assert set(get_neighborhood(g, "lonely", induced=False)) == {"lonely"}
