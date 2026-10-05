"""Tests for skm_tools.load_networks using small real-data fixtures."""

import gzip
from unittest.mock import patch

import networkx as nx
import pytest

from skm_tools.load_networks import (
    ckn_to_networkx,
    pss_gene_network_to_networkx,
    pss_interaction_network_to_networkx,
    pss_reaction_graph_to_networkx,
)


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
        with patch("skm_tools.load_networks.urlretrieve") as mock_dl:
            ckn_to_networkx(ckn_edge_path, ckn_node_path)
            mock_dl.assert_not_called()

    def test_accepts_string_paths(self, ckn_edge_path, ckn_node_path):
        g = ckn_to_networkx(str(ckn_edge_path), str(ckn_node_path))
        assert isinstance(g, nx.DiGraph)


# ---------------------------------------------------------------------------
# PSS exports (pss-export): shared behaviour of all three loaders
# ---------------------------------------------------------------------------

PSS_LOADERS = {
    "reaction_graph": pss_reaction_graph_to_networkx,
    "interaction_network": pss_interaction_network_to_networkx,
    "gene_network_ath": pss_gene_network_to_networkx,
}


@pytest.fixture(params=list(PSS_LOADERS))
def pss_export(request):
    """(loader, edge path, node path) for each PSS export"""
    name = request.param
    return (
        PSS_LOADERS[name],
        request.getfixturevalue(f"pss_{name}_edge_path"),
        request.getfixturevalue(f"pss_{name}_node_path"),
    )


def _rows(path):
    with open(path) as handle:
        return [line.rstrip("\n").split("\t") for line in handle][1:]


class TestPSSExports:

    def test_returns_multidigraph(self, pss_export):
        loader, edge_path, node_path = pss_export
        assert isinstance(loader(edge_path, node_path), nx.MultiDiGraph)

    def test_all_nodes_and_edges_loaded(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        assert g.number_of_nodes() == len(_rows(node_path))
        assert g.number_of_edges() == len(_rows(edge_path))

    def test_accepts_string_paths(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(str(edge_path), str(node_path))
        assert g.number_of_edges() > 0

    def test_node_attributes_set(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, data in g.nodes(data=True):
            assert data["node_type"]
            assert isinstance(data["display_label"], str) and data["display_label"]
            assert "id" not in data

    def test_edge_endpoints_not_kept_as_attributes(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, _, data in g.edges(data=True):
            assert "source" not in data and "target" not in data

    def test_empty_values_are_none(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        values = [v for _, d in g.nodes(data=True) for v in d.values()]
        values += [v for _, _, d in g.edges(data=True) for v in d.values()]
        assert "" not in values
        assert None in values

    def test_list_columns_are_list_or_none(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, data in g.nodes(data=True):
            for k, v in data.items():
                if k in ("synonyms", "all_pathways", "external_links", "components") or k.endswith("_homologues"):
                    assert v is None or isinstance(v, list)

    def test_lists_split_on_semicolon_not_comma(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        synonyms = [s for _, d in g.nodes(data=True) for s in (d["synonyms"] or [])]
        assert len(synonyms) > g.number_of_nodes() / 2
        assert not any(";" in s for s in synonyms)
        # names with commas stay whole (10,11-EHT, AHK2,3,4, ...)
        assert any("," in s for s in synonyms)

    def test_complex_components_are_lists(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        complexes = [d for _, d in g.nodes(data=True) if d["node_type"] == "Complex"]
        assert complexes
        assert all(isinstance(d["components"], list) for d in complexes)

    def test_booleans_parsed(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, _, data in g.edges(data=True):
            assert isinstance(data["directed"], bool)

    def test_missing_file_raises(self, pss_export, tmp_path):
        loader, edge_path, node_path = pss_export
        with pytest.raises(FileNotFoundError):
            loader(tmp_path / "missing-edges.tsv", node_path)

    def test_no_download_attempted(self, pss_export, tmp_path):
        loader, edge_path, node_path = pss_export
        with patch("skm_tools.load_networks.urlretrieve") as mock_dl:
            with pytest.raises(FileNotFoundError):
                loader(tmp_path / "missing-edges.tsv", node_path)
            mock_dl.assert_not_called()


# ---------------------------------------------------------------------------
# PSS exports: per-format specifics
# ---------------------------------------------------------------------------

class TestPSSReactionGraph:

    def test_reactions_are_nodes(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        reactions = [n for n, d in g.nodes(data=True) if d["node_type"] == "reaction"]
        assert reactions
        for r in reactions:
            assert g.nodes[r]["reaction_type"]

    def test_edges_keyed_by_role(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["role"]

    def test_same_entity_twice_in_one_reaction(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        # WRKY33 is both template and stimulator of its own expression (rx00160)
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        assert set(g["WRKY33[fc00166]"]["rx00160"]) == {"template", "stimulator"}


class TestPSSInteractionNetwork:

    def test_edges_keyed_by_reaction_id(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["reaction_id"]

    def test_mutual_edges_listed_both_ways(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        mutual = [(u, v, k) for u, v, k, d in g.edges(keys=True, data=True) if not d["directed"]]
        assert mutual
        for u, v, k in mutual:
            assert g.has_edge(v, u, key=k)

    def test_homologues_loaded(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        plant = [d for _, d in g.nodes(data=True) if d["node_type"] == "PlantCoding"]
        assert plant
        assert any(d["ath_homologues"] for d in plant)


class TestPSSGeneNetwork:

    def test_genes_have_species_and_entity(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        genes = [d for _, d in g.nodes(data=True) if d["node_type"] == "gene"]
        assert genes
        for d in genes:
            assert d["species"] == "ath"
            assert d["entity"]

    def test_cluster_columns_are_aligned_lists(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        multi = [d for _, d in g.nodes(data=True) if len(d["functional_cluster_id"] or []) > 1]
        assert multi
        for _, d in g.nodes(data=True):
            n = len(d["entity"])
            for c in ("short_name", "pathway", "functional_cluster_id"):
                assert d[c] is None or len(d[c]) == n

    def test_edges_keyed_by_reaction_id(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["reaction_id"]
