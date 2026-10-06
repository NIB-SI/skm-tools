"""Tests for skm_tools.ckn using small real-data fixtures."""

from unittest.mock import patch

import networkx as nx
import pandas as pd
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
        assert node_attrs == {"node_type", "locus_type", "species", "TAIR", "display_label", "short_name",
                              "synonyms", "description", "mapman", "note", "tissue"}
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
    metabolites = [d for _, d in g.nodes(data=True) if d["node_type"] == "Metabolite"]
    assert metabolites
    assert all(d["species"] is None for d in metabolites)


def test_empty_values_are_none(ckn_edge_path, ckn_node_path):
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    values = [v for _, d in g.nodes(data=True) for v in d.values()]
    assert None in values
    assert not any(isinstance(v, float) and v != v for v in values)  # no NaN


def test_filter_ckn_nodes_species_keeps_nodes_without_species(ckn_edge_path, ckn_node_path):
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    g_species = dict(g.nodes(data="species"))
    metabolites = {n for n, d in g.nodes(data=True) if d["node_type"] == "Metabolite"}

    reasons = filter_ckn_nodes(g, species=["ath"], remove_isolates=False)

    assert metabolites <= set(g)
    # only nodes of other species are removed (e.g. "foreign": pathogens, abiotic stresses)
    for n, reason in reasons.items():
        assert reason != "wrong species" or g_species[n] not in (None, "ath")


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


# ---------------------------------------------------------------------------
# node_type: PSS classes, with the TAIR locus type of genes and RNAs in locus_type
# ---------------------------------------------------------------------------

PSS_CLASSES = {"PlantCoding", "PlantNonCoding", "PlantPseudogene", "Metabolite", "Complex", "Process",
               "ForeignCoding", "ForeignNonCoding", "ForeignEntity", "ForeignAbiotic"}


def test_node_types_are_pss_classes(ckn_edge_path, ckn_node_path):
    # v2.0.1 files have them; older files (v2) are converted when loaded
    g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
    for n, d in g.nodes(data=True):
        assert d["node_type"] in PSS_CLASSES, n
        is_locus = d["node_type"] in ("PlantCoding", "PlantNonCoding", "PlantPseudogene")
        assert (d["locus_type"] is not None) == is_locus, n


def test_v2_0_1_fixture_has_every_node_type():
    g = ckn_to_networkx(FIXTURES / "ckn_v2.0.1_edges.tsv", FIXTURES / "ckn_v2.0.1_nodes.tsv.gz")
    assert {d["node_type"] for _, d in g.nodes(data=True)} == PSS_CLASSES


def test_older_node_types_converted():
    from skm_tools.ckn import _ckn_pss_node_types
    df = pd.DataFrame({
        "id": ["AT1G01010", "MIR165A", "AT1G16140", "ABA", "bacteria", "virus_CP", "abiotic_heat"],
        "node_type": ["protein_coding", "mirna", "pseudogene", "metabolite", "biotic", "biotic", "abiotic"],
        "species": ["ath", "ath", "ath", None, "foreign", "foreign", "foreign"],
    })
    out = _ckn_pss_node_types(df)
    assert list(out.columns[:3]) == ["id", "node_type", "locus_type"]
    assert out["node_type"].tolist() == ["PlantCoding", "PlantNonCoding", "PlantPseudogene", "Metabolite",
                                         "ForeignEntity", "ForeignCoding", "ForeignAbiotic"]
    assert out["locus_type"].tolist() == ["protein_coding", "mirna", "pseudogene", None, None, None, None]


def test_older_node_types_unknown_biotic_node_raises():
    from skm_tools.ckn import _ckn_pss_node_types
    df = pd.DataFrame({"id": ["virus_new"], "node_type": ["biotic"]})
    with pytest.raises(ValueError, match="virus_new"):
        _ckn_pss_node_types(df)
