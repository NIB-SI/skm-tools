import json

import networkx as nx
import numpy as np
import pytest

from skm_tools.persistence import load_graph, save_graph
from skm_tools.serialize import graph_from_dict, graph_to_dict, path_to_dict, to_json_safe


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


def test_path_to_dict():
    assert path_to_dict([["A", "B", "C"], ("A", "C")]) == {"paths": [["A", "B", "C"], ["A", "C"]], "lengths": [2, 1]}


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


def test_graphml_flattens_lists_and_drops_none(tmp_path):
    g = _multigraph()
    path = tmp_path / "g.graphml"
    save_graph(g, path, format="graphml")
    h = load_graph(path, format="graphml")
    assert h.is_multigraph() and h.is_directed()
    assert h.number_of_edges() == 2
    assert h.nodes["A"]["synonyms"] == "a;a1"
    assert "pathway" not in h.nodes["A"]
    assert "synonyms" in g.nodes["A"] and g.nodes["A"]["synonyms"] == ["a", "a1"]  # g unchanged


def test_graphml_digraph(tmp_path):
    path = tmp_path / "g.graphml"
    save_graph(nx.DiGraph([("A", "B")]), path, format="graphml")
    h = load_graph(path, format="graphml")
    assert type(h) is nx.DiGraph


def test_bad_format(tmp_path):
    with pytest.raises(ValueError):
        save_graph(nx.DiGraph(), tmp_path / "g", format="csv")
