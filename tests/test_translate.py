"""Tests for skm_tools.translate"""

import gzip
from unittest.mock import patch

import networkx as nx
import pandas as pd
import pytest

from skm_tools.translate import integrate_translation_ckn, load_translation_file

TABLE = "ath_source\tstu\nAT1G01010\tSoltu.DM.01G000010\n"


def test_load_translation_file(tmp_path):
    path = tmp_path / "stu.tsv"
    path.write_text(TABLE)
    assert load_translation_file("stu", path).to_dict("records") == [
        {"ath_source": "AT1G01010", "stu": "Soltu.DM.01G000010"}]


def test_load_translation_file_gzipped(tmp_path):
    # as downloaded from skm.nib.si
    path = tmp_path / "stu.tsv.gz"
    path.write_bytes(gzip.compress(TABLE.encode()))
    assert load_translation_file("stu", path)["stu"].tolist() == ["Soltu.DM.01G000010"]


def test_load_translation_file_na_like_ids(tmp_path):
    path = tmp_path / "stu.tsv"
    path.write_text("ath_source\tpotato\nAT1G01010\tNA\nAT1G01020\tnull\nAT1G01030\t\n")
    df = load_translation_file("stu", path)
    assert df["potato"].tolist() == ["NA", "null", None]


def test_load_translation_file_default_name(tmp_path):
    with patch("skm_tools.utils.download",
               side_effect=lambda url, path: path.write_bytes(gzip.compress(TABLE.encode()))) as download:
        load_translation_file("stu", data_dir=tmp_path)
    download.assert_called_once_with("https://skm.nib.si/downloads/translations/stu",
                                     tmp_path / "translation_ath_to_stu.tsv.gz")


# ---------------------------------------------------------------------------
# integrate_translation_ckn
# ---------------------------------------------------------------------------

def _ckn():
    """A -> B -> C (Arabidopsis genes), C -> A|B (complex), M (metabolite) -> A"""
    g = nx.DiGraph([("A", "B"), ("B", "C"), ("C", "A|B"), ("M", "A")], interaction="positive-influence")
    for n in "ABC":
        g.add_node(n, species="ath", node_type="PlantCoding", TAIR=f"AT{n}")
    g.add_node("A|B", species="ath", node_type="Complex", TAIR=None)
    g.add_node("M", species=None, node_type="Metabolite", TAIR=None)
    return g


TRANSLATIONS = pd.DataFrame({"ath_source": ["ATA", "ATA", "ATB", "ATB"],
                             "potato": ["St1", "St2", "St2", "St3"]})


def test_integrate_translation_many_to_many():
    h = integrate_translation_ckn(_ckn(), TRANSLATIONS, "potato",
                                  edges_within_translation=False, edges_across_translation=False)
    assert set(h) == {"A_St1", "A_St2", "B_St2", "B_St3", "C", "A|B", "M"}
    # edges copied to every translation
    assert set(h.edges()) == {(a, b) for a in ("A_St1", "A_St2") for b in ("B_St2", "B_St3")} | \
        {("B_St2", "C"), ("B_St3", "C"), ("C", "A|B"), ("M", "A_St1"), ("M", "A_St2")}
    assert h.nodes["A_St1"]["translation"] == "St1" and h.nodes["A_St1"]["translated_from"] == "A"
    assert h.nodes["C"]["translated"] is False


def test_integrate_translation_keeps_complexes():
    h = integrate_translation_ckn(_ckn(), TRANSLATIONS, "potato")
    assert "A|B" in h and "translated" not in h.nodes["A|B"]


def test_integrate_translation_drop_unmapped_with_edges():
    h = integrate_translation_ckn(_ckn(), TRANSLATIONS, "potato", keep_unmapped=False)
    assert "C" not in h
    assert all("node_type" in d for _, d in h.nodes(data=True))


def test_integrate_translation_edges():
    h = integrate_translation_ckn(_ckn(), TRANSLATIONS, "potato")
    homology = [(u, v, d) for u, v, d in h.edges(data=True) if d.get("type") == "homology"]
    relations = {(u, v): d["translation_relation"] for u, v, d in homology}
    # translations of one gene, and genes of different originals with the same translation:
    # both directions, as undirected edges
    assert relations[("A_St1", "A_St2")] == relations[("A_St2", "A_St1")] == "same_translation_source"
    assert relations[("A_St2", "B_St2")] == relations[("B_St2", "A_St2")] == "same_translation_target"
    for *_, d in homology:
        assert d["interaction"] == "homology" and d["directed"] is False
    assert type(h) is nx.MultiDiGraph


def test_integrate_translation_needs_target_column():
    with pytest.raises(TypeError):
        integrate_translation_ckn(_ckn(), TRANSLATIONS)
