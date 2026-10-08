import json

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from skm_tools.persistence import load_graph, save_graph
from skm_tools.serialize import graph_from_dict, graph_to_dict, paths_to_dict, to_json_safe


def _multigraph():
    g = nx.MultiDiGraph(name="test")
    g.add_node("A", node_type="PlantCoding", synonyms=["a", "a1"], pathway=None)
    g.add_node("B", node_type="Metabolite", score=np.float64(0.5))
    g.add_edge("A", "B", key="rx1", interaction="positive-influence", directed=True)
    g.add_edge("A", "B", key="rx2", interaction="negative-influence", directed=True)
    return g


def test_to_json_safe():
    value = {"a": {1, 3, 2}, "b": (1, np.int64(2)), "c": float("nan"), 4: np.array([1.0, 2.0]), "d": object}
    safe = to_json_safe(value)
    json.dumps(safe)
    assert safe["a"] == [1, 2, 3]
    assert safe["b"] == [1, 2]
    assert safe["c"] is None
    assert safe["4"] == [1.0, 2.0]


def test_graph_to_dict_is_json_and_keeps_keys():
    d = graph_to_dict(_multigraph())
    json.dumps(d)
    assert d["directed"] and d["multigraph"]
    assert d["graph"] == {"name": "test"}
    assert {e["key"] for e in d["edges"]} == {"rx1", "rx2"}
    assert list(d["edges"][0])[:3] == ["source", "target", "key"]
    assert d["nodes"][0]["id"] == "A"


def test_graph_to_dict_simple_graph_has_no_keys():
    d = graph_to_dict(nx.DiGraph([("A", "B")]))
    assert d["edges"] == [{"source": "A", "target": "B"}]
    assert not d["multigraph"]


def test_dict_round_trip():
    g = _multigraph()
    h = graph_from_dict(json.loads(json.dumps(graph_to_dict(g))))
    assert type(h) is nx.MultiDiGraph
    assert set(h.edges(keys=True)) == set(g.edges(keys=True))
    assert h.nodes["A"]["synonyms"] == ["a", "a1"]
    assert h.graph["name"] == "test"


def test_to_json_safe_pandas_missing_values():
    assert to_json_safe([pd.NA, pd.NaT, np.nan]) == [None, None, None]


def test_to_json_safe_key_collision():
    with pytest.raises(ValueError, match="same as strings"):
        to_json_safe({1: "a", "1": "b"})


@pytest.mark.parametrize("graph_type", [nx.Graph, nx.DiGraph, nx.MultiGraph])
def test_round_trip_graph_types_and_int_ids(graph_type):
    g = graph_type()
    g.add_edge(1, 2, weight=0.5)
    h = graph_from_dict(json.loads(json.dumps(graph_to_dict(g))))
    assert type(h) is graph_type
    assert set(h.nodes) == {1, 2}
    assert h.number_of_edges() == 1


def test_round_trip_create_using():
    h = graph_from_dict(graph_to_dict(_multigraph()), create_using=nx.MultiGraph)
    assert type(h) is nx.MultiGraph


def test_graph_to_dict_rejects_tuple_ids():
    with pytest.raises(ValueError, match="strings or integers"):
        graph_to_dict(nx.Graph([(("A", 1), "B")]))


@pytest.mark.parametrize("node_attrs, edge_attrs, match", [
    ({"id": "TAIR:1"}, {}, "Node attributes \\['id'\\]"),
    ({}, {"source": "lit"}, "Edge attributes \\['source'\\]"),
])
def test_graph_to_dict_rejects_attribute_name_clashes(node_attrs, edge_attrs, match):
    g = nx.DiGraph()
    g.add_node("A", **node_attrs)
    g.add_edge("A", "B", **edge_attrs)
    with pytest.raises(ValueError, match=match):
        graph_to_dict(g)


def test_paths_to_dict():
    assert paths_to_dict([["A", "B", "C"], ("A", "C")]) == {"paths": [["A", "B", "C"], ["A", "C"]], "lengths": [2, 1]}


def test_paths_to_dict_with_graph():
    g = nx.MultiDiGraph()
    g.add_node("A", display_label="a")
    g.add_node("B", display_label="b")
    g.add_edge("A", "B", interaction="positive-influence")
    g.add_edge("A", "B", interaction="negative-influence")
    g.add_edge("C", "B", interaction="positive-influence")
    d = paths_to_dict([["A", "B", "C"]], g)
    assert d["nodes"] == {"A": {"display_label": "a"}, "B": {"display_label": "b"},
                          "C": {"display_label": None}}
    assert d["steps"] == [[
        {"source": "A", "target": "B", "interaction": ["positive-influence", "negative-influence"]},
        {"source": "B", "target": "C", "reversed": True, "interaction": ["positive-influence"]},
    ]]
    json.dumps(d)


@pytest.mark.parametrize("format", ["pickle", "json"])
def test_save_load_round_trip(tmp_path, format):
    g = _multigraph()
    path = tmp_path / f"g.{format}"
    save_graph(g, path, format=format)
    h = load_graph(str(path), format=format)
    assert type(h) is nx.MultiDiGraph
    assert set(h.edges(keys=True)) == set(g.edges(keys=True))
    assert h.nodes["A"]["synonyms"] == ["a", "a1"]
    assert h.nodes["A"]["pathway"] is None


@pytest.mark.parametrize("format", ["csv", "graphml"])
def test_bad_format(tmp_path, format):
    with pytest.raises(ValueError):
        save_graph(nx.DiGraph(), tmp_path / "g", format=format)
