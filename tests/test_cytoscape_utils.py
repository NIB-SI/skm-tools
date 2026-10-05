"""Cytoscape helpers, with py4cytoscape mocked: no running Cytoscape needed."""

from unittest.mock import patch

import networkx as nx
import pandas as pd
import pytest

pytest.importorskip("py4cytoscape")

from skm_tools import cytoscape_utils as cu


NODE_TABLE = pd.DataFrame({"name": ["A", "B", "C", "EIN3(like)", "X (1) Y"]}, index=[1, 2, 3, 4, 5])
EDGE_TABLE = pd.DataFrame(
    {"name": [
        "A (positive-influence) B",
        "A (negative-influence) B",
        "B (interacts with) C",
        "EIN3(like) (unknown-influence) X (1) Y",
        "C (positive-influence) A",
    ]},
    index=[11, 12, 13, 14, 15],
)


def _get_table_columns(table="node", columns=None, network=None):
    return NODE_TABLE if table == "node" else EDGE_TABLE


@pytest.fixture
def p4c():
    with patch.object(cu, "p4c") as mock:
        mock.tables.get_table_columns.side_effect = _get_table_columns
        mock.styles.get_current_style.return_value = "PSS-default"
        yield mock


def test_match_edge_names_any_interaction_and_parallel_edges():
    names = EDGE_TABLE["name"].to_dict()
    assert cu._match_edge_names(names, [("A", "B")]) == {11: ("A", "B"), 12: ("A", "B")}
    assert cu._match_edge_names(names, [("B", "C", "rx1")]) == {13: ("B", "C")}
    assert cu._match_edge_names(names, [("B", "A")]) == {}


def test_match_edge_names_with_brackets_in_node_names():
    names = EDGE_TABLE["name"].to_dict()
    assert cu._match_edge_names(names, [("EIN3(like)", "X (1) Y")]) == {14: ("EIN3(like)", "X (1) Y")}


def test_highlight_edges_uses_suids_and_skips(p4c):
    edges = cu.highlight_edges([("A", "B"), ("B", "C"), ("C", "A")], "#FF0000", skip_edges=[("C", "A")], network=7)
    assert edges == [("A", "B"), ("B", "C")]
    suids = p4c.style_bypasses.set_edge_color_bypass.call_args.args[0]
    assert sorted(suids) == [11, 12, 13]


def test_highlight_path_returns_nodes_and_edges(p4c):
    nodes, edges = cu.highlight_path(["A", "B", "C"], "#00FF00", skip_nodes=["A"])
    assert nodes == ["A", "B", "C"]
    assert edges == [("A", "B"), ("B", "C")]
    assert p4c.style_bypasses.set_node_color_bypass.call_args.args[0] == [2, 3]


def test_highlight_nodes_only_sets_given_properties(p4c):
    assert cu.highlight_nodes("C", colour="#000000") == [3]
    p4c.style_bypasses.set_node_color_bypass.assert_called_once()
    p4c.style_bypasses.set_node_border_width_bypass.assert_not_called()


def test_subnetwork_edge_induced_from_paths(p4c):
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    cu.subnetwork_edge_induced_from_paths([["A", "B", "C"]], g, parent_suid=99, name="paths")
    kwargs = p4c.networks.create_subnetwork.call_args.kwargs
    assert sorted(kwargs["nodes"]) == [1, 2, 3]
    assert sorted(kwargs["edges"]) == [11, 12, 13]
    assert kwargs["exclude_edges"] is True
    assert kwargs["network"] == 99


def test_subnetwork_node_induced(p4c):
    cu.subnetwork_node_induced(["A", "C", "not-there"], parent_suid=99)
    assert p4c.networks.create_subnetwork.call_args.kwargs["nodes"] == [1, 3]


def test_contrast_colour():
    assert cu.contrast_colour("#000000") == "#FFFFFF"
    assert cu.contrast_colour("#123456") == "#EDCBA9"
