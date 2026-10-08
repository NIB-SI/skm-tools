"""Tests for skm_tools.ckn using small real-data fixtures."""

import gzip
import shutil
from unittest.mock import patch

import networkx as nx
import pytest

from skm_tools.ckn import CKN_RANKS, ckn_to_networkx, filter_ckn_edges, filter_ckn_nodes, rank_counts
from skm_tools.skm_download_urls import CKN_EDGE_FILE, CKN_EDGE_URL, CKN_NODE_FILE, CKN_NODE_URL

from .conftest import FIXTURES

EDGES = FIXTURES / "ckn_v2.0.1_edges.tsv"
NODES = FIXTURES / "ckn_v2.0.1_nodes.tsv.gz"


@pytest.fixture
def ckn():
    return ckn_to_networkx(EDGES, NODES)


def _write_ckn(tmp_path, edge_rows, node_rows):
    '''Small CKN files: rows as lists of values, in the v2.0.1 columns.'''
    edge_header = ["source", "interaction", "target", "directed", "rank", "effect", "type", "species",
                   "isTFregulation", "interactionSources"]
    node_header = ["id", "node_type", "locus_type", "species", "TAIR", "display_label", "short_name",
                   "synonyms", "description", "mapman", "note", "tissue"]
    edges, nodes = tmp_path / "edges.tsv", tmp_path / "nodes.tsv"
    edges.write_text("\n".join("\t".join(r) for r in [edge_header] + edge_rows) + "\n")
    nodes.write_text("\n".join("\t".join(r) for r in [node_header] + node_rows) + "\n")
    return edges, nodes


def _gene(n, **values):
    row = {"id": n, "node_type": "PlantCoding", "locus_type": "protein_coding", "species": "ath",
           "TAIR": n, "display_label": n, "short_name": n, "synonyms": "", "description": "",
           "mapman": "", "note": "", "tissue": "leaf"}
    row.update(values)
    return list(row.values())


def _edge(u, v, directed="True", rank="1"):
    return [u, "positive-influence", v, directed, rank, "act", "binding", "ath", "0", "skm"]


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

class TestCKNToNetworkx:

    def test_returns_digraph(self, ckn):
        assert type(ckn) is nx.DiGraph
        assert ckn.number_of_nodes() > 0 and ckn.number_of_edges() > 0

    def test_attributes(self, ckn):
        node_attrs = {k for _, d in ckn.nodes(data=True) for k in d}
        edge_attrs = {k for *_, d in ckn.edges(data=True) for k in d}
        assert node_attrs == {"node_type", "locus_type", "species", "TAIR", "display_label", "short_name",
                              "synonyms", "description", "mapman", "note", "tissue"}
        assert edge_attrs == {"interaction", "directed", "rank", "effect", "type", "species",
                              "isTFregulation", "interactionSources"}
        for *_, d in ckn.edges(data=True):
            assert isinstance(d["rank"], int)
            assert isinstance(d["directed"], bool)
            assert d["interaction"] in ("positive-influence", "negative-influence", "unknown-influence")

    def test_list_attributes(self, ckn):
        for _, data in ckn.nodes(data=True):
            for k in ("synonyms", "mapman", "tissue"):
                assert data[k] is None or isinstance(data[k], list)
        assert any(len(d["tissue"] or []) > 1 for _, d in ckn.nodes(data=True))
        for *_, data in ckn.edges(data=True):
            assert isinstance(data["interactionSources"], list)

    def test_display_label(self, ckn):
        for n, d in ckn.nodes(data=True):
            assert d["display_label"] == (d["short_name"] or n)

    def test_as_directed_adds_reverse_of_undirected_edges(self, ckn):
        undirected = [(u, v) for u, v, d in ckn.edges(data=True) if not d["directed"]]
        assert undirected
        for u, v in undirected:
            assert ckn.has_edge(v, u)

    def test_as_directed_false_adds_edges_as_listed(self):
        g = ckn_to_networkx(EDGES, NODES, as_directed=False)
        assert g.number_of_edges() == len(EDGES.read_text().splitlines()) - 1
        assert g.number_of_edges() < ckn_to_networkx(EDGES, NODES).number_of_edges()

    def test_multidigraph(self):
        g = ckn_to_networkx(EDGES, NODES, create_using=nx.MultiDiGraph)
        assert type(g) is nx.MultiDiGraph

    def test_accepts_string_paths(self):
        assert isinstance(ckn_to_networkx(str(EDGES), str(NODES)), nx.DiGraph)

    def test_no_download_when_files_exist(self):
        with patch("skm_tools.utils.download") as download:
            ckn_to_networkx(EDGES, NODES)
        download.assert_not_called()

    def test_default_paths_download_to_data_dir(self, tmp_path):
        def fake_download(url, path):
            # as utils.download: gzipped, as the names end in .gz
            source = {CKN_EDGE_URL: EDGES, CKN_NODE_URL: NODES}[url]
            if source.suffix == ".gz":
                shutil.copyfile(source, path)
            else:
                with open(source, "rb") as src, gzip.open(path, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        with patch("skm_tools.utils.download", side_effect=fake_download) as download:
            g = ckn_to_networkx(data_dir=tmp_path)
        download.assert_any_call(CKN_EDGE_URL, tmp_path / CKN_EDGE_FILE)
        download.assert_any_call(CKN_NODE_URL, tmp_path / CKN_NODE_FILE)
        assert g.number_of_nodes() > 0

    def test_older_format_raises(self, tmp_path):
        nodes = tmp_path / "nodes.tsv"
        with gzip.open(NODES, "rt") as handle:
            lines = handle.read().splitlines()
        # drop the locus_type column (the format before CKN v2.0.1)
        nodes.write_text("\n".join("\t".join(l.split("\t")[:2] + l.split("\t")[3:]) for l in lines))
        with pytest.raises(ValueError, match="not a CKN v2.0.1 node file"):
            ckn_to_networkx(EDGES, nodes)


def test_metabolite_species_is_missing(ckn):
    metabolites = [d for _, d in ckn.nodes(data=True) if d["node_type"] == "Metabolite"]
    assert metabolites
    assert all(d["species"] is None for d in metabolites)


def test_empty_values_are_none(ckn):
    values = [v for _, d in ckn.nodes(data=True) for v in d.values()]
    assert None in values
    assert not any(isinstance(v, float) and v != v for v in values)  # no NaN


def test_gene_symbols_with_semicolon_stay_whole(ckn):
    assert "PIP1;3" in ckn.nodes["AT1G01620"]["synonyms"]


def test_quotes_are_part_of_values(tmp_path):
    # CKN v2.0.1 is written without quoting: a value can start with '"'
    edges, nodes = _write_ckn(tmp_path, [_edge("A", "B")],
                              [_gene("A", synonyms='"sunken"|low'), _gene("B", description='"x" and y')])
    g = ckn_to_networkx(edges, nodes)
    assert g.nodes["A"]["synonyms"] == ['"sunken"', "low"]
    assert g.nodes["B"]["description"] == '"x" and y'


def test_missing_directed_is_none(tmp_path):
    edges, nodes = _write_ckn(tmp_path, [_edge("A", "B", directed="")], [_gene("A"), _gene("B")])
    g = ckn_to_networkx(edges, nodes)
    assert g["A"]["B"]["directed"] is None
    assert not g.has_edge("B", "A")  # not added in reverse: not known to be undirected


# ---------------------------------------------------------------------------
# node_type: PSS classes, with the TAIR locus type of genes and RNAs in locus_type
# ---------------------------------------------------------------------------

PSS_CLASSES = {"PlantCoding", "PlantNonCoding", "PlantPseudogene", "Metabolite", "Complex", "Process",
               "ForeignCoding", "ForeignNonCoding", "ForeignEntity", "ForeignAbiotic"}


def test_node_types_are_pss_classes(ckn):
    for n, d in ckn.nodes(data=True):
        assert d["node_type"] in PSS_CLASSES, n
        is_locus = d["node_type"] in ("PlantCoding", "PlantNonCoding", "PlantPseudogene")
        assert (d["locus_type"] is not None) == is_locus, n


def test_fixture_has_every_node_type(ckn):
    assert {d["node_type"] for _, d in ckn.nodes(data=True)} == PSS_CLASSES


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def test_rank_counts(ckn):
    counts = rank_counts(ckn)
    assert list(counts)[:len(CKN_RANKS)] == list(CKN_RANKS)
    assert sum(counts.values()) == ckn.number_of_edges()


def test_filter_ckn_edges_ranks_and_return_value(ckn):
    removed_edges, removed_nodes = filter_ckn_edges(ckn, keep_edge_ranks=[0])
    assert removed_edges
    assert all(d["rank"] == 0 for *_, d in ckn.edges(data=True))
    assert all(r == "isolate" for r in removed_nodes.values())
    assert not set(removed_nodes) & set(ckn)


@pytest.mark.parametrize("ranks", [(0,), {0}, range(1), 0])
def test_filter_ckn_edges_any_iterable(ranks):
    g = nx.DiGraph([("A", "B", {"rank": 0}), ("B", "C", {"rank": 3})])
    filter_ckn_edges(g, keep_edge_ranks=ranks)
    assert list(g.edges()) == [("A", "B")]


def test_filter_ckn_edges_multigraph_removes_the_right_parallel_edge():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", rank=0, type="binding")
    g.add_edge("A", "B", rank=3, type="binding")
    removed_edges, _ = filter_ckn_edges(g, keep_edge_ranks=[3])
    assert [d["rank"] for *_, d in g.edges(data=True)] == [3]
    assert removed_edges == [("A", "B", 0)]


def test_filter_ckn_edges_types_and_function():
    g = nx.DiGraph([("A", "B", {"type": "binding", "directed": False}),
                    ("B", "C", {"type": "binding", "directed": True}),
                    ("C", "D", {"type": "other", "directed": True})])
    filter_ckn_edges(g, keep_edge_types="binding", filter_function=lambda d: d["directed"])
    assert list(g.edges()) == [("B", "C")]
    assert set(g) == {"B", "C"}


def test_filter_ckn_nodes_species_keeps_nodes_without_species(ckn):
    species = dict(ckn.nodes(data="species"))
    metabolites = {n for n, d in ckn.nodes(data=True) if d["node_type"] == "Metabolite"}

    reasons = filter_ckn_nodes(ckn, species=["ath"], remove_isolates=False)

    assert metabolites <= set(ckn)
    # only nodes of other species are removed (e.g. "foreign": pathogens, abiotic stresses)
    for n, reason in reasons.items():
        assert reason != "wrong species" or species[n] not in (None, "ath")


def test_filter_ckn_nodes_string_arguments(ckn):
    h = ckn.copy()
    assert filter_ckn_nodes(ckn, species="ath", tissues="leaf") == \
        filter_ckn_nodes(h, species=["ath"], tissues=["leaf"])
    assert set(ckn) == set(h)


def test_filter_ckn_nodes_mixed_complex_removed_with_its_foreign_component(ckn):
    assert ckn.nodes["RISC|virus_vsiRNA"]["species"] == "ath/foreign"
    reasons = filter_ckn_nodes(ckn, species="ath")
    assert "RISC|virus_vsiRNA" in reasons


def test_filter_ckn_nodes_complex_component_removed():
    g = nx.DiGraph([("A", "A|B"), ("B", "A|B"), ("A|B", "C")])
    nx.set_node_attributes(g, {"A": "PlantCoding", "B": "Metabolite", "A|B": "Complex", "C": "PlantCoding"},
                           "node_type")
    nx.set_node_attributes(g, None, "short_name")
    reasons = filter_ckn_nodes(g, node_types=["PlantCoding", "Complex"])
    assert reasons == {"B": "wrong node type", "A|B": "complex component removed",
                       "A": "isolate", "C": "isolate"}
