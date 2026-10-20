import networkx as nx
import numpy as np
import pandas as pd
import pytest

from skm_tools.experimental_data import add_experimental_data, add_node_data


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


def test_in_place_by_default(df):
    g = _graph()
    h = add_experimental_data(g, df, "logFC", "padj")
    assert h is g
    assert g.nodes["G1"]["logFC"] == 1.5


def test_copy(df):
    g = _graph()
    h = add_experimental_data(g, df, "logFC", "padj", copy=True)
    assert h is not g
    assert h.nodes["G1"]["logFC"] == 1.5
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
    h2 = add_experimental_data(h, df2, "logFC", "padj", copy=True)
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


# add_node_data

@pytest.fixture
def gene_df():
    return pd.DataFrame(
        {"maxLog10P": [5.0, 7.5, np.nan], "trait": ["Rust", "Blight", None],
         "traits": [["Rust"], ["Blight", "Rust"], []]},
        index=["P1", "P2", "P3"],
    )


def _translated():
    # two Arabidopsis genes translated to the same apricot gene P1, as integrate_translation_ckn does
    g = nx.MultiDiGraph([("A1_P1", "A2_P1"), ("A2_P1", "A3_P2"), ("A3_P2", "M")])
    for n in ("A1_P1", "A2_P1", "A3_P2"):
        g.nodes[n]["translation"] = n.split("_")[1]
    g.add_node("FC", homologues=["P2", "P1", "X"])
    return g


def test_node_data_by_node_id(gene_df):
    g = nx.Graph([("P1", "P2"), ("P2", "N")])
    h = add_node_data(g, gene_df)
    assert h is g
    assert h.nodes["P2"] == {"maxLog10P": 7.5, "trait": "Blight", "traits": ["Blight", "Rust"],
                             "matched_ids": ["P2"]}
    assert h.nodes["N"] == {}


def test_node_data_same_id_on_several_nodes(gene_df):
    h = add_node_data(_translated(), gene_df, match_attribute="translation", prefix="gwas ")
    for n in ("A1_P1", "A2_P1"):
        assert h.nodes[n]["gwas maxLog10P"] == 5.0
        assert h.nodes[n]["gwas matched_ids"] == ["P1"]
    assert h.nodes["A3_P2"]["gwas trait"] == "Blight"
    assert "gwas trait" not in h.nodes["M"]


def test_node_data_columns_and_types(gene_df):
    h = add_node_data(_translated(), gene_df, columns=["maxLog10P"], match_attribute="translation")
    assert "trait" not in h.nodes["A1_P1"]
    assert type(h.nodes["A1_P1"]["maxLog10P"]) is float


def test_node_data_nan_is_none(gene_df):
    g = nx.Graph([("P3", "N")])
    h = add_node_data(g, gene_df)
    assert h.nodes["P3"]["maxLog10P"] is None and h.nodes["P3"]["trait"] is None
    assert h.nodes["P3"]["traits"] == []


def test_node_data_list_attribute_rank_by(gene_df):
    g = _translated()
    # default: first matching row in df order
    assert add_node_data(g, gene_df, match_attribute="homologues").nodes["FC"]["trait"] == "Rust"
    h = add_node_data(g, gene_df, match_attribute="homologues", rank_by="maxLog10P")
    assert h.nodes["FC"]["trait"] == "Blight"
    assert h.nodes["FC"]["matched_ids"] == ["P1", "P2"]
    h = add_node_data(g, gene_df, match_attribute="homologues", rank_by="maxLog10P", ascending=True)
    assert h.nodes["FC"]["trait"] == "Rust"


def test_node_data_copy_and_reannotate(gene_df):
    g = _translated()
    h = add_node_data(g, gene_df, match_attribute="translation", copy=True)
    assert "trait" not in g.nodes["A1_P1"]
    h2 = add_node_data(h, gene_df.loc[["P2"]], match_attribute="translation")
    assert not {"maxLog10P", "trait", "matched_ids"} & set(h2.nodes["A1_P1"])
    assert h2.nodes["A3_P2"]["trait"] == "Blight"


def test_node_data_duplicated_index_raises(gene_df):
    with pytest.raises(ValueError, match="aggregate"):
        add_node_data(_translated(), pd.concat([gene_df, gene_df]))


def test_node_data_missing_column_raises(gene_df):
    with pytest.raises(KeyError):
        add_node_data(_translated(), gene_df, columns=["pvalue"])
    with pytest.raises(KeyError):
        add_node_data(_translated(), gene_df, rank_by="pvalue")


def test_node_data_no_match_warns(gene_df, caplog):
    add_node_data(_translated(), gene_df)
    assert "No node matched the df index" in caplog.text
