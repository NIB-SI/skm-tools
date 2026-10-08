import networkx as nx
import numpy as np
import pandas as pd
import pytest

from skm_tools.experimental_data import add_experimental_data


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
    h = add_experimental_data(g, df, "logFC", "padj")
    assert h is not g
    assert "logFC" not in g.nodes["G1"]


def test_values_by_node_id(df):
    h = add_experimental_data(_graph(), df, "logFC", "padj")
    assert h.nodes["G1"]["logFC"] == 1.5
    assert h.nodes["G1"]["pvalue"] == 0.01
    assert h.nodes["G1"]["significant"] is True
    assert h.nodes["G2"]["significant"] is False
    assert h.nodes["G1"]["matched_ids"] == ["G1"]
    assert "logFC" not in h.nodes["N"]


def test_python_types(df):
    h = add_experimental_data(_graph(), df, "logFC", "padj")
    assert type(h.nodes["G1"]["logFC"]) is float


def test_match_attribute_picks_most_significant(df):
    h = add_experimental_data(_graph(), df, "logFC", "padj", match_attribute="ath_homologues")
    assert h.nodes["FC"]["logFC"] == 1.5  # G1 has the lowest p-value
    assert h.nodes["FC"]["matched_ids"] == ["G1", "G2", "G3"]  # sorted
    assert "logFC" not in h.nodes["G1"]  # G1 has no ath_homologues attribute


def test_without_pvalue_picks_largest_abs_logfc(df):
    h = add_experimental_data(_graph(), df, "logFC", match_attribute="ath_homologues")
    assert h.nodes["FC"]["logFC"] == 3.0
    assert "pvalue" not in h.nodes["FC"] and "significant" not in h.nodes["FC"]


def test_nan_is_none_and_prefix(df):
    g = _graph()
    g.add_node("G4")
    h = add_experimental_data(g, df, "logFC", "padj", prefix="t10 ")
    assert h.nodes["G4"]["t10 logFC"] is None
    assert h.nodes["G4"]["t10 significant"] is True


def test_cutoff(df):
    h = add_experimental_data(_graph(), df, "logFC", "padj", cutoff=0.001)
    assert h.nodes["G1"]["significant"] is False


def test_missing_column_raises(df):
    with pytest.raises(KeyError):
        add_experimental_data(_graph(), df, "log2FoldChange", "padj")


def test_duplicated_ids_use_the_most_significant_row(caplog):
    df = pd.DataFrame({"logFC": [5.0, 1.0, 2.0], "padj": [0.5, 0.001, 0.001]}, index=["G1", "G1", "G1"])
    with caplog.at_level("INFO", logger="skm_tools"):
        h = add_experimental_data(_graph(), df, "logFC", "padj")
    # lowest p-value, then largest |logFC|
    assert (h.nodes["G1"]["logFC"], h.nodes["G1"]["pvalue"]) == (2.0, 0.001)
    assert "1 identifiers are in df more than once" in caplog.text


def test_duplicated_ids_without_pvalue():
    df = pd.DataFrame({"logFC": [1.0, -4.0, np.nan]}, index=["G1", "G1", "G1"])
    assert add_experimental_data(_graph(), df, "logFC").nodes["G1"]["logFC"] == -4.0


def test_reannotating_removes_old_values(df):
    h = add_experimental_data(_graph(), df, "logFC", "padj")
    df2 = pd.DataFrame({"logFC": [0.5], "padj": [0.2]}, index=["G2"])
    h2 = add_experimental_data(h, df2, "logFC", "padj")
    assert h2.nodes["G2"]["logFC"] == 0.5
    assert not {"logFC", "pvalue", "significant", "matched_ids"} & set(h2.nodes["G1"])
    # another prefix is kept
    h3 = add_experimental_data(h, df2, "logFC", "padj", prefix="other ")
    assert h3.nodes["G1"]["logFC"] == 1.5


def test_int_node_ids():
    g = nx.Graph([(1, 2)])
    df = pd.DataFrame({"logFC": [1.0]}, index=[2])
    assert add_experimental_data(g, df, "logFC").nodes[2]["logFC"] == 1.0


def test_nullable_dtypes():
    df = pd.DataFrame({"logFC": pd.array([1.0, None], dtype="Float64"),
                       "padj": pd.array([None, 0.01], dtype="Float64")}, index=["G1", "G2"])
    h = add_experimental_data(_graph(), df, "logFC", "padj")
    assert h.nodes["G1"]["pvalue"] is None and h.nodes["G1"]["significant"] is False
    assert h.nodes["G2"]["logFC"] is None and h.nodes["G2"]["significant"] is True
    assert type(h.nodes["G1"]["logFC"]) is float


def test_no_match_warns(df, caplog):
    df.index = df.index.str.lower()
    add_experimental_data(_graph(), df, "logFC", "padj")
    assert "No node matched the df index" in caplog.text
