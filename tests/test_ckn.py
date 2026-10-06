"""Tests for skm_tools.ckn using small real-data fixtures."""

from unittest.mock import patch

import networkx as nx
import pytest

import pandas as pd

from skm_tools.ckn import ckn_to_networkx, filter_ckn_nodes


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
        assert "rank" in attrs
        assert "isDirected" in attrs

    def test_gmm_attribute_is_list_or_none(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        for _, data in g.nodes(data=True):
            assert data["GMM"] is None or isinstance(data["GMM"], list)

    def test_add_reciprocal_edges_default(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=True)
        # every undirected edge should have a reciprocal
        for u, v, data in g.edges(data=True):
            if data["isDirected"] == 0:
                assert g.has_edge(v, u), f"Missing reciprocal edge for undirected ({u}, {v})"

    def test_add_reciprocal_edges_false(self, ckn_edge_path, ckn_node_path):
        g_with = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=True)
        g_without = ckn_to_networkx(ckn_edge_path, ckn_node_path, add_reciprocal_edges=False)
        assert g_with.number_of_edges() >= g_without.number_of_edges()

    def test_directed_removes_undirected_edges(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path, directed=True)
        for _, _, data in g.edges(data=True):
            assert data["isDirected"] == 1

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

def test_metabolite_species_na_is_missing(ckn_edge_path, ckn_node_path):
    # the CKN node file has species "N/A" for metabolites
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    metabolites = [d for _, d in g.nodes(data=True) if d["node_type"] == "metabolite"]
    assert metabolites
    assert all(pd.isna(d["species"]) for d in metabolites)


def test_filter_ckn_nodes_species_keeps_nodes_without_species(ckn_edge_path, ckn_node_path):
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    metabolites = {n for n, d in g.nodes(data=True) if d["node_type"] == "metabolite"}

    reasons = filter_ckn_nodes(g, species=["ath"], remove_isolates=False)

    assert metabolites <= set(g)
    assert "wrong species" not in reasons.values()
