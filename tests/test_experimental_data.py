import networkx as nx
import numpy as np
import pandas as pd
import pytest

from skm_tools.experimental_data import overlay_experimental_data


@pytest.fixture
def df():
    return pd.DataFrame(
        {"logFC": [1.5, -0.2, 3.0, np.nan], "padj": [0.01, 0.5, 0.2, 0.03]},
        index=["G1", "G2", "G3", "G4"],
    )


def _graph():
    g = nx.DiGraph([("G1", "G2"), ("G2", "N")])
    g.add_node("FC", ath_homologues=["G2", "G3", "G1"])
    return g


def test_returns_copy(df):
    g = _graph()
    h = overlay_experimental_data(g, df, "logFC", "padj")
    assert h is not g
    assert "logFC" not in g.nodes["G1"]


def test_values_by_node_id(df):
    h = overlay_experimental_data(_graph(), df, "logFC", "padj")
    assert h.nodes["G1"]["logFC"] == 1.5
    assert h.nodes["G1"]["pvalue"] == 0.01
    assert h.nodes["G1"]["significant"] is True
    assert h.nodes["G2"]["significant"] is False
    assert h.nodes["G1"]["matched_ids"] == ["G1"]
    assert "logFC" not in h.nodes["N"]


def test_python_types(df):
    h = overlay_experimental_data(_graph(), df, "logFC", "padj")
    assert type(h.nodes["G1"]["logFC"]) is float


def test_match_attribute_picks_most_significant(df):
    h = overlay_experimental_data(_graph(), df, "logFC", "padj", match_attribute="ath_homologues")
    assert h.nodes["FC"]["logFC"] == 1.5  # G1 has the lowest p-value
    assert h.nodes["FC"]["matched_ids"] == ["G2", "G3", "G1"]
    assert "logFC" not in h.nodes["G1"]  # G1 has no ath_homologues attribute


def test_without_pvalue_picks_largest_abs_logfc(df):
    h = overlay_experimental_data(_graph(), df, "logFC", match_attribute="ath_homologues")
    assert h.nodes["FC"]["logFC"] == 3.0
    assert "pvalue" not in h.nodes["FC"] and "significant" not in h.nodes["FC"]


def test_nan_is_none_and_prefix(df):
    g = _graph()
    g.add_node("G4")
    h = overlay_experimental_data(g, df, "logFC", "padj", prefix="t10 ")
    assert h.nodes["G4"]["t10 logFC"] is None
    assert h.nodes["G4"]["t10 significant"] is True


def test_cutoff(df):
    h = overlay_experimental_data(_graph(), df, "logFC", "padj", cutoff=0.001)
    assert h.nodes["G1"]["significant"] is False


def test_missing_column_raises(df):
    with pytest.raises(KeyError):
        overlay_experimental_data(_graph(), df, "log2FoldChange", "padj")
