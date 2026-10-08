"""Tests for skm_tools.node_images (no Cytoscape needed)."""

import logging
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from skm_tools import node_images as ni


@pytest.fixture(autouse=True)
def _info_logs(caplog):
    caplog.set_level(logging.INFO, logger="skm_tools")


def test_node_file_key():
    assert ni.node_file_key("WRKY33[fc00166]") == "WRKY33_fc00166_"
    assert ni.node_file_key("11-/12-OH-JA") == "11-_12-OH-JA"
    assert ni.node_file_key("AT2G38470") == "AT2G38470"


def test_match_files_to_nodes(tmp_path, caplog):
    for f in ["WRKY33_fc00166_.svg", "11-_12-OH-JA.png", "Pro.png", "WRKY33.png", "notes.txt"]:
        (tmp_path / f).write_text("x")
    nodes = ["WRKY33[fc00166]", "WRKY33[fc00999]", "11-/12-OH-JA", "Proline accumulation", "JA"]

    matched = ni.match_files_to_nodes(tmp_path, nodes, aliases={"Pro": "Proline accumulation"})

    assert {n: p.name for n, p in matched.items()} == {
        "WRKY33[fc00166]": "WRKY33_fc00166_.svg",
        "11-/12-OH-JA": "11-_12-OH-JA.png", "Proline accumulation": "Pro.png",
    }
    assert "WRKY33.png" in caplog.text


def test_match_files_to_nodes_custom_key(tmp_path):
    (tmp_path / "WRKY33.png").write_text("x")
    matched = ni.match_files_to_nodes(tmp_path, ["WRKY33[fc00166]"], key=lambda n: n.partition("[")[0])
    assert matched == {"WRKY33[fc00166]": tmp_path / "WRKY33.png"}


def test_unique_image_copies_new_names_every_call(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    (src / "WRKY33.png").write_bytes(b"png")
    images = {"A": src / "WRKY33.png", "B": src / "WRKY33.png"}

    first = ni.unique_image_copies(images, tmp_path / "out")
    second = ni.unique_image_copies(images, tmp_path / "out")

    assert first["A"] == first["B"]   # one copy per image
    assert first["A"] != second["A"]  # new names on every call
    assert first["A"].name.startswith("WRKY33_") and first["A"].suffix == ".png"
    assert first["A"].read_bytes() == b"png"
    assert len(list((tmp_path / "out").iterdir())) == 2


def test_unique_image_copies_converts_svg(tmp_path):
    (tmp_path / "a.svg").write_text("<svg/>")
    cairosvg = MagicMock()
    with patch.dict("sys.modules", {"cairosvg": cairosvg}):
        copies = ni.unique_image_copies({"A": tmp_path / "a.svg"}, tmp_path / "out", to_png=True)
    assert copies["A"].suffix == ".png"
    cairosvg.svg2png.assert_called_once()


def test_chart_column():
    df = pd.DataFrame({"t1": [1.5, None], "t2": [-0.5, 2.0]}, index=["A", "B"])
    charts = ni.chart_column(df, ["t1", "t2"], ["#E41A1C", "#377EB8"], value_range=(-2, 2),
                             labels=["10 min", "30 min"], separation=2)
    assert charts["B"] == ('barchart: colorlist="#E41A1C,#377EB8" valuelist="0,2" '
                           'labellist="10 min,30 min" range="-2,2" separation=2')




def test_chart_column_rejects_commas_and_quotes_in_labels():
    df = pd.DataFrame({"t1": [1.5]}, index=["A"])
    for label in ["10 min, mock", 'the "mock"']:
        with pytest.raises(ValueError, match="can't contain"):
            ni.chart_column(df, ["t1"], "#E41A1C", labels=[label])
