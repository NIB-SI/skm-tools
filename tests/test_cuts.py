"""Tests for skm_tools.cuts"""

import networkx as nx
import pytest

from skm_tools.cuts import get_cutset


def _graph(edges, capacity=1):
    g = nx.DiGraph(edges)
    nx.set_edge_attributes(g, capacity, "capacity")
    return g


def _assert_min_cut(g, sources, targets, cut, size):
    '''`cut` has `size` edges of `g`, and disconnects every source from every target
    (with ties, which of the minimum cuts is returned is up to networkx).'''
    assert len(cut) == size and all(g.has_edge(u, v) for u, v in cut)
    h = g.copy()
    h.remove_edges_from(cut)
    assert not any(nx.has_path(h, s, t) for s in sources for t in targets)


def test_bottleneck_edge():
    g = _graph([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D"), ("D", "E")])
    assert get_cutset(["A"], ["E"], g) == [("D", "E")]


def test_parallel_paths_all_cut(capsys):
    g = _graph([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")])
    _assert_min_cut(g, ["A"], ["D"], get_cutset(["A"], ["D"], g), 2)
    assert "max_flow = 2" in capsys.readouterr().out


def test_capacity_is_used():
    g = _graph([("A", "B"), ("B", "C")])
    g["B"]["C"]["capacity"] = 5
    assert get_cutset(["A"], ["C"], g) == [("A", "B")]


def test_several_sources_and_targets():
    g = _graph([("S1", "X"), ("S2", "X"), ("X", "T1"), ("X", "T2"), ("S2", "T2")])
    cut = get_cutset(["S1", "S2"], ["T1", "T2"], g)
    # max flow 3: S1 -> X -> T1, S2 -> X -> T2, S2 -> T2
    _assert_min_cut(g, ["S1", "S2"], ["T1", "T2"], cut, 3)


def test_cut_disconnects_sources_from_targets():
    g = _graph([("A", "B"), ("B", "C"), ("A", "C"), ("C", "D"), ("B", "D"), ("D", "E"), ("C", "E")])
    _assert_min_cut(g, ["A"], ["E"], get_cutset(["A"], ["E"], g), 2)


def test_no_path_gives_empty_cut(capsys):
    g = _graph([("A", "B"), ("C", "D")])
    assert get_cutset(["A"], ["D"], g) == []
    assert "max_flow = 0" in capsys.readouterr().out


def test_missing_nodes_and_source_target_overlap_ignored():
    g = _graph([("A", "B"), ("B", "C")])
    g["B"]["C"]["capacity"] = 5
    assert get_cutset(["A", "missing"], ["C", "A", "also missing"], g) == [("A", "B")]


def test_node_names_like_source_and_sink():
    g = _graph([("source", "x"), ("x", "sink")])
    g["x"]["sink"]["capacity"] = 5
    assert get_cutset(["source"], ["sink"], g) == [("source", "x")]


def test_graph_not_changed():
    g = _graph([("A", "B"), ("B", "C")])
    before = (set(g.nodes()), set(g.edges()))
    get_cutset(["A"], ["C"], g)
    assert (set(g.nodes()), set(g.edges())) == before


def test_missing_capacity_raises():
    g = nx.DiGraph([("A", "B"), ("B", "C")])
    with pytest.raises(nx.NetworkXUnbounded):
        get_cutset(["A"], ["C"], g)
