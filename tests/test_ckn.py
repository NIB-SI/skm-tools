"""Tests for skm_tools.ckn using small real-data fixtures."""

import gzip
import shutil
from unittest.mock import patch

import networkx as nx
import numpy as np
import pytest

from skm_tools.ckn import CKN_RANKS, ckn_to_networkx, filter_ckn_edges, filter_ckn_nodes, rank_counts
from skm_tools.skm_download_urls import CKN_EDGE_FILE, CKN_EDGE_URL, CKN_NODE_FILE, CKN_NODE_URL

from .conftest import FIXTURES

EDGES = FIXTURES / "ckn_v2.1_edges.tsv"
NODES = FIXTURES / "ckn_v2.1_nodes.tsv.gz"


@pytest.fixture
def ckn():
    return ckn_to_networkx(EDGES, NODES)


def _write_ckn(tmp_path, edge_rows, node_rows):
    '''Small CKN files: rows as lists of values, in the v2.1 columns.'''
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
        # drop the locus_type column (the format before CKN v2.1)
        nodes.write_text("\n".join("\t".join(l.split("\t")[:2] + l.split("\t")[3:]) for l in lines))
        with pytest.raises(ValueError, match="not a CKN v2.1 node file"):
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
    # CKN v2.1 is written without quoting: a value can start with '"'
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

PSS_CLASSES = {"PlantCoding", "PlantNonCoding", "PlantPseudogene", "PlantAbstract", "Metabolite",
               "MetaboliteFamily", "Complex", "Process", "ForeignCoding", "ForeignNonCoding", "ForeignEntity",
               "ForeignAbiotic"}


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
    assert ckn.nodes["RISC - vsiRNA"]["species"] == "ath/foreign"
    reasons = filter_ckn_nodes(ckn, species="ath")
    assert "RISC - vsiRNA" in reasons


def test_filter_ckn_nodes_complex_component_removed():
    g = nx.DiGraph([("A", "A|B"), ("B", "A|B"), ("A|B", "C")])
    nx.set_node_attributes(g, {"A": "PlantCoding", "B": "Metabolite", "A|B": "Complex", "C": "PlantCoding"},
                           "node_type")
    nx.set_node_attributes(g, None, "short_name")
    reasons = filter_ckn_nodes(g, node_types=["PlantCoding", "Complex"])
    assert reasons == {"B": "wrong node type", "A|B": "complex component removed",
                       "A": "isolate", "C": "isolate"}


class TestTranscriptionFactors:
    def test_sources_of_tf_regulation_edges(self):
        from skm_tools.ckn import TF_REGULATION, transcription_factors
        g = nx.DiGraph()
        g.add_edge("TF1", "G1", type=TF_REGULATION)
        g.add_edge("TF1", "G2", type=TF_REGULATION)
        g.add_edge("K1", "TF1", type="post-translational modification")
        g.add_node("TF2", mapman=["14.5.7.1_RNA biosynthesis.DNA-binding transcriptional regulation.NAC"])
        g.add_node("G3", mapman=["14.50_something else"])
        assert transcription_factors(g) == {"TF1"}
        assert transcription_factors(g, mapman_bins="14.5") == {"TF1", "TF2"}

    def test_fixture(self, ckn_edge_path, ckn_node_path):
        from skm_tools.ckn import transcription_factors
        g = ckn_to_networkx(ckn_edge_path, ckn_node_path)
        tfs = transcription_factors(g)
        assert "AT2G43010" in tfs
        assert all(any(d["type"] == "transcription factor regulation" for *_, d in g.out_edges(t, data=True))
                   for t in tfs)


class TestTfTargetEnrichment:
    @staticmethod
    def graph():
        from skm_tools.ckn import TF_REGULATION
        g = nx.MultiDiGraph()
        # TF1 regulates G0-G9 (activating), TF2 regulates G10-G29, TF3 has only 2 targets
        g.add_edges_from((("TF1", f"G{i}", dict(type=TF_REGULATION, interaction="positive-influence"))
                          for i in range(10)))
        g.add_edges_from((("TF2", f"G{i}", dict(type=TF_REGULATION, interaction="unknown-influence"))
                          for i in range(10, 30)))
        g.add_edges_from((("TF3", f"G{i}", dict(type=TF_REGULATION)) for i in (30, 31)))
        g.add_edge("TF1", "G0", type="binding", interaction="unknown-influence")     # not a regulation edge
        g.add_nodes_from(f"N{i}" for i in range(50))                                 # not regulated
        return g

    def test_enriched_tf_and_background(self):
        from skm_tools.ckn import tf_target_enrichment
        degs = {f"G{i}": -1 for i in range(8)} | {"G12": 1, "N1": 1}
        r = tf_target_enrichment(self.graph(), degs)
        assert list(r.tf) == ["TF1", "TF2"]                  # TF3: fewer than 5 targets
        assert r.attrs == {"background": 32, "n_degs": 9}    # N1 is not regulated: not in the background
        tf1 = r.iloc[0]
        assert (tf1.targets, tf1.hits, tf1.down_hits, tf1.up_hits) == (10, 8, 8, 0)
        assert tf1.fold == pytest.approx(8 / (10 * 9 / 32))
        assert tf1.padj < 0.01
        # activating edges to down-regulated targets: a less active TF
        assert (tf1.n_signed, tf1.activity_score) == (8, -1.0)
        assert r.iloc[1].n_signed == 0 and np.isnan(r.iloc[1].activity_score)

    def test_universe_and_plain_node_list(self):
        from skm_tools.ckn import tf_target_enrichment
        r = tf_target_enrichment(self.graph(), [f"G{i}" for i in range(5)], universe=[f"G{i}" for i in range(20)],
                                 min_targets=1)
        assert r.attrs["background"] == 20
        assert set(r.tf) == {"TF1", "TF2"}
        assert r.set_index("tf").at["TF1", "up_hits"] == 0      # no directions given
