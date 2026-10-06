"""Tests for skm_tools.ckn using small real-data fixtures."""

from unittest.mock import patch

import networkx as nx
import pytest

from skm_tools.ckn import ckn_to_networkx, filter_ckn_nodes

from .conftest import FIXTURES


# ---------------------------------------------------------------------------
# CKN
# ---------------------------------------------------------------------------

class TestCKNToNetworkx:

    def test_returns_digraph(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        assert isinstance(g, nx.DiGraph)

    def test_nodes_loaded(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        assert g.number_of_nodes() > 0

    def test_edges_loaded(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        assert g.number_of_edges() > 0

    def test_node_attributes_set(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        node = next(iter(g.nodes()))
        attrs = g.nodes[node]
        assert "node_type" in attrs
        assert "species" in attrs

    def test_edge_attributes_set(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        u, v = next(iter(g.edges()))
        attrs = g[u][v]
        assert "effect" in attrs
        assert isinstance(attrs["rank"], int)
        assert isinstance(attrs["directed"], bool)
        assert attrs["interaction"] in ("positive-influence", "negative-influence", "unknown-influence")

    def test_same_attributes_in_both_formats(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        node_attrs = {k for _, d in g.nodes(data=True) for k in d}
        edge_attrs = {k for *_, d in g.edges(data=True) for k in d}
        assert node_attrs == {"node_type", "species", "TAIR", "display_label", "short_name", "synonyms",
                              "description", "mapman", "note", "tissue"}
        assert edge_attrs == {"interaction", "directed", "rank", "effect", "type", "species",
                              "isTFregulation", "interactionSources"}

    def test_list_attributes(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        for _, data in g.nodes(data=True):
            for k in ("synonyms", "mapman", "tissue"):
                assert data[k] is None or isinstance(data[k], list)
        assert any(len(d["tissue"] or []) > 1 for _, d in g.nodes(data=True))
        for *_, data in g.edges(data=True):
            assert isinstance(data["interactionSources"], list)

    def test_display_label(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        for n, d in g.nodes(data=True):
            assert d["display_label"] == (d["short_name"] or n)

    def test_add_reciprocal_edges_default(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=True)
        # every undirected edge should have a reciprocal
        for u, v, data in g.edges(data=True):
            if not data["directed"]:
                assert g.has_edge(v, u), f"Missing reciprocal edge for undirected ({u}, {v})"

    def test_add_reciprocal_edges_false(self, ckn_edge_path, ckn_node_path):
        g_with = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=True)
        g_without = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=False)
        assert g_with.number_of_edges() >= g_without.number_of_edges()

    def test_directed_removes_undirected_edges(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path, directed=True)
        for _, _, data in g.edges(data=True):
            assert data["directed"]

    def test_directed_removes_isolates(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path, directed=True)
        assert list(nx.isolates(g)) == []

    def test_no_download_when_files_exist(self, ckn_edge_path, ckn_node_path):
        with patch("skm_tools.ckn.urlretrieve") as mock_dl:
            ckn_to_networkx(ckn_edge_path, ckn_node_path)
            mock_dl.assert_not_called()

    def test_accepts_string_paths(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(str(ckn_edge_path), str(ckn_node_path))
        assert isinstance(g, nx.DiGraph)


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def test_metabolite_species_is_missing(ckn_edge_path, ckn_node_path):
    # the CKN v2 node file has species "N/A" for metabolites, v2.0.1 an empty value
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    metabolites = [d for _, d in g.nodes(data=True) if d["node_type"] == "metabolite"]
    assert metabolites
    assert all(d["species"] is None for d in metabolites)


def test_empty_values_are_none(ckn_edge_path, ckn_node_path):
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    values = [v for _, d in g.nodes(data=True) for v in d.values()]
    assert None in values
    assert not any(isinstance(v, float) and v != v for v in values)  # no NaN


def test_filter_ckn_nodes_species_keeps_nodes_without_species(ckn_edge_path, ckn_node_path):
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    metabolites = {n for n, d in g.nodes(data=True) if d["node_type"] == "metabolite"}

    reasons = filter_ckn_nodes(g, species=["ath"], remove_isolates=False)

    assert metabolites <= set(g)
    assert "wrong species" not in reasons.values()


# ---------------------------------------------------------------------------
# Format specifics
# ---------------------------------------------------------------------------

def test_gene_symbols_with_semicolon_stay_whole():
    g = ckn_to_networkx(FIXTURES / "ckn_v2.0.1_edges.tsv", FIXTURES / "ckn_v2.0.1_nodes.tsv.gz")
    assert "PIP1;3" in g.nodes["AT1G01620"]["synonyms"]


def test_v2_effect_converted_to_interaction():
    g = ckn_to_networkx(FIXTURES / "ckn_v2_edges.tsv", FIXTURES / "ckn_v2_nodes.tsv.gz", add_reciprocal_edges=False)
    expected = {"act": "positive-influence", "inh": "negative-influence"}
    for *_, d in g.edges(data=True):
        assert d["interaction"] == expected.get(d["effect"], "unknown-influence")
