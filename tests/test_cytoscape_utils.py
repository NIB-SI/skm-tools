"""Cytoscape helpers, with py4cytoscape mocked: no running Cytoscape needed."""

from unittest.mock import MagicMock, patch

import networkx as nx
import pandas as pd
import pytest

pytest.importorskip("py4cytoscape")

from skm_tools import cytoscape_utils as cu


NODE_TABLE = pd.DataFrame({"name": ["A", "B", "C", "EIN3(like)", "X (1) Y"]}, index=[1, 2, 3, 4, 5])
EDGE_TABLE = pd.DataFrame(
    {"name": [
        "A (positive-influence) B",
        "A (negative-influence) B",
        "B (interacts with) C",
        "EIN3(like) (unknown-influence) X (1) Y",
        "C (positive-influence) A",
    ]},
    index=[11, 12, 13, 14, 15],
)


def _get_table_columns(table="node", columns=None, network=None):
    return NODE_TABLE if table == "node" else EDGE_TABLE


@pytest.fixture
def p4c():
    with patch.object(cu, "p4c") as mock:
        mock.tables.get_table_columns.side_effect = _get_table_columns
        mock.styles.get_current_style.return_value = "SKM"
        yield mock


def test_match_edge_names_any_interaction_and_parallel_edges():
    names = EDGE_TABLE["name"].to_dict()
    assert cu._match_edge_names(names, [("A", "B")]) == {11: ("A", "B"), 12: ("A", "B")}
    assert cu._match_edge_names(names, [("B", "C", "rx1")]) == {13: ("B", "C")}
    assert cu._match_edge_names(names, [("B", "A")]) == {}


def test_match_edge_names_with_brackets_in_node_names():
    names = EDGE_TABLE["name"].to_dict()
    assert cu._match_edge_names(names, [("EIN3(like)", "X (1) Y")]) == {14: ("EIN3(like)", "X (1) Y")}


def test_highlight_edges_uses_suids_and_skips(p4c):
    edges = cu.highlight_edges([("A", "B"), ("B", "C"), ("C", "A")], "#FF0000", skip_edges=[("C", "A")], network=7)
    assert edges == [("A", "B"), ("B", "C")]
    suids = p4c.style_bypasses.set_edge_color_bypass.call_args.args[0]
    assert sorted(suids) == [11, 12, 13]


def test_highlight_path_returns_nodes_and_edges(p4c):
    nodes, edges = cu.highlight_path(["A", "B", "C"], "#00FF00", skip_nodes=["A"])
    assert nodes == ["A", "B", "C"]
    assert edges == [("A", "B"), ("B", "C")]
    assert p4c.style_bypasses.set_node_color_bypass.call_args.args[0] == [2, 3]


def test_highlight_nodes_only_sets_given_properties(p4c):
    assert cu.highlight_nodes("C", colour="#000000") == [3]
    p4c.style_bypasses.set_node_color_bypass.assert_called_once()
    p4c.style_bypasses.set_node_border_width_bypass.assert_not_called()


def test_subnetwork_edge_induced_from_paths(p4c):
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    cu.subnetwork_edge_induced_from_paths([["A", "B", "C"]], g, parent_suid=99, name="paths")
    kwargs = p4c.networks.create_subnetwork.call_args.kwargs
    assert sorted(kwargs["nodes"]) == [1, 2, 3]
    assert sorted(kwargs["edges"]) == [11, 12, 13]
    assert kwargs["exclude_edges"] is True
    assert kwargs["network"] == 99


def test_subnetwork_node_induced(p4c):
    cu.subnetwork_node_induced(["A", "C", "not-there"], parent_suid=99)
    assert p4c.networks.create_subnetwork.call_args.kwargs["nodes"] == [1, 3]


def test_contrast_colour():
    assert cu.contrast_colour("#000000") == "#FFFFFF"
    assert cu.contrast_colour("#123456") == "#EDCBA9"


# ---------------------------------------------------------------------------
# Networks, collections and styles
# ---------------------------------------------------------------------------

def test_clone_network_renames_collection_then_network(p4c):
    p4c.networks.clone_network.return_value = 50
    p4c.collections.get_collection_suid.return_value = 49

    assert cu.clone_network(7, name="Heat", collection="Collection - Heat") == 50

    op = p4c.commands.cyrest_put.call_args
    assert op.args[0] == "collections/49/tables/default"
    assert op.kwargs["body"]["data"] == [{"SUID": 49, "name": "Collection - Heat"}]
    p4c.rename_network.assert_called_once_with("Heat", network=50)
    names = [c[0] for c in p4c.mock_calls]
    assert names.index("commands.cyrest_put") < names.index("rename_network")


def test_clone_network_without_names(p4c):
    p4c.networks.clone_network.return_value = 50
    cu.clone_network(7)
    p4c.commands.cyrest_put.assert_not_called()
    p4c.rename_network.assert_not_called()


def test_copy_style_applies_to_networks(p4c):
    assert cu.copy_style("SKM", "SKM-heat", networks=[1, 2]) == "SKM-heat"
    p4c.copy_visual_style.assert_called_once_with("SKM", "SKM-heat")
    assert [c.kwargs["network"] for c in p4c.set_visual_style.call_args_list] == [1, 2]


def test_delete_other_networks(p4c):
    p4c.get_network_list.return_value = [{"name": "a", "suid": 1}, {"name": "b", "suid": 2},
                                         {"name": "c", "suid": 3}]
    assert cu.delete_other_networks(2) == [1, 3]
    assert [c.args[0] for c in p4c.delete_network.call_args_list] == [1, 3]
    p4c.set_current_network.assert_called_once_with(2)


def test_get_or_create_subnetwork(p4c):
    p4c.get_network_list.return_value = ["existing"]
    p4c.get_network_suid.return_value = 5
    p4c.networks.create_subnetwork.return_value = 6

    assert cu.get_or_create_subnetwork(["A"], 1, "existing") == 5
    p4c.networks.create_subnetwork.assert_not_called()

    assert cu.get_or_create_subnetwork(["A", "B"], 1, "new") == 6
    assert sorted(p4c.networks.create_subnetwork.call_args.kwargs["nodes"]) == [1, 2]


# ---------------------------------------------------------------------------
# Node images and charts
# ---------------------------------------------------------------------------

def test_node_file_key():
    assert cu.node_file_key("WRKY33[fc00166]") == "WRKY33_fc00166_"
    assert cu.node_file_key("11-/12-OH-JA") == "11-_12-OH-JA"
    assert cu.node_file_key("AT2G38470") == "AT2G38470"


def test_match_files_to_nodes(tmp_path, capsys):
    for f in ["WRKY33_fc00166_.svg", "11-_12-OH-JA.png", "Pro.png", "WRKY33.png", "notes.txt"]:
        (tmp_path / f).write_text("x")
    nodes = ["WRKY33[fc00166]", "WRKY33[fc00999]", "11-/12-OH-JA", "Proline accumulation", "JA"]

    matched = cu.match_files_to_nodes(tmp_path, nodes, aliases={"Pro": "Proline accumulation"})

    assert {n: p.name for n, p in matched.items()} == {
        "WRKY33[fc00166]": "WRKY33_fc00166_.svg",
        "11-/12-OH-JA": "11-_12-OH-JA.png", "Proline accumulation": "Pro.png",
    }
    assert "WRKY33.png" in capsys.readouterr().out


def test_match_files_to_nodes_custom_key(tmp_path):
    (tmp_path / "WRKY33.png").write_text("x")
    matched = cu.match_files_to_nodes(tmp_path, ["WRKY33[fc00166]"], key=lambda n: n.partition("[")[0])
    assert matched == {"WRKY33[fc00166]": tmp_path / "WRKY33.png"}


def test_unique_image_copies_new_names_every_call(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    (src / "WRKY33.png").write_bytes(b"png")
    images = {"A": src / "WRKY33.png", "B": src / "WRKY33.png"}

    first = cu.unique_image_copies(images, tmp_path / "out")
    second = cu.unique_image_copies(images, tmp_path / "out")

    assert first["A"] == first["B"]   # one copy per image
    assert first["A"] != second["A"]  # new names on every call
    assert first["A"].name.startswith("WRKY33_") and first["A"].suffix == ".png"
    assert first["A"].read_bytes() == b"png"
    assert len(list((tmp_path / "out").iterdir())) == 2


def test_unique_image_copies_converts_svg(tmp_path):
    (tmp_path / "a.svg").write_text("<svg/>")
    cairosvg = MagicMock()
    with patch.dict("sys.modules", {"cairosvg": cairosvg}):
        copies = cu.unique_image_copies({"A": tmp_path / "a.svg"}, tmp_path / "out", to_png=True)
    assert copies["A"].suffix == ".png"
    cairosvg.svg2png.assert_called_once()


def test_load_node_images(p4c, tmp_path):
    (tmp_path / "a.png").write_bytes(b"png")
    images = cu.load_node_images({"A": tmp_path / "a.png"}, "image_heat", network=7,
                                 unique_dir=tmp_path / "unique")

    assert images["A"].parent == tmp_path / "unique" and images["A"].name.startswith("a_")
    table = p4c.load_table_data.call_args.args[0]
    assert table.loc["A", "image_heat"] == f"file:{images['A']}"
    assert p4c.load_table_data.call_args.kwargs["network"] == 7


def test_show_node_images(p4c):
    cu.show_node_images("SKM-heat", "image_heat", slot=2, position="above", size=90)

    p4c.style_mappings.map_visual_property.assert_called_once_with(
        visual_prop="NODE_CUSTOMGRAPHICS_2", table_column="image_heat", mapping_type="p")
    defaults = {c.args[0]["visualProperty"]: c.args[0]["value"]
                for c in p4c.style_defaults.set_visual_property_default.call_args_list}
    assert defaults == {"NODE_CUSTOMGRAPHICS_POSITION_2": "N,S,c,0.00,0.00",
                        "NODE_CUSTOMGRAPHICS_SIZE_2": 90}


def test_add_custom_png(p4c, tmp_path):
    (tmp_path / "A.png").write_bytes(b"png")
    p4c.get_all_nodes.return_value = ["A", "B"]

    cu.add_custom_png(7, lambda n: tmp_path / f"{n}.png" if n == "A" else None)

    assert list(p4c.load_table_data.call_args.args[0].index) == ["A"]
    p4c.style_dependencies.sync_node_custom_graphics_size.assert_called_once_with(False, style_name="SKM")


def test_chart_column():
    df = pd.DataFrame({"t1": [1.5, None], "t2": [-0.5, 2.0]}, index=["A", "B"])
    charts = cu.chart_column(df, ["t1", "t2"], ["#E41A1C", "#377EB8"], value_range=(-2, 2),
                             labels=["10 min", "30 min"], separation=2)
    assert charts["B"] == ('barchart: colorlist="#E41A1C,#377EB8" valuelist="0,2" '
                           'labellist="10 min,30 min" range="-2,2" separation=2')


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def test_export_collection(p4c, tmp_path):
    p4c.get_collection_networks.return_value = [12, 11]
    p4c.get_network_name.side_effect = lambda suid: {11: "JA - Heat", 12: "SA/Heat"}[suid]

    with patch.object(cu.time, "sleep") as sleep:
        files = cu.export_collection(11, tmp_path / "out", zoom=300)
    assert sleep.call_count == 2  # waits for Cytoscape before each export

    assert [f.name for f in files] == ["JA_Heat_11.png", "SA_Heat_12.png"]
    call = p4c.network_views.export_image.call_args.kwargs
    assert call["zoom"] == 300 and call["all_graphics_details"] is True and call["type"] == "PNG"


def test_clone_network_retries_rename(p4c):
    import py4cytoscape
    p4c.CyError = py4cytoscape.CyError
    p4c.networks.clone_network.return_value = 50
    p4c.rename_network.side_effect = [py4cytoscape.CyError("unrecognized (table entry)"), None]

    with patch.object(cu, "_RETRY_WAIT", 0):
        assert cu.clone_network(7, name="Heat") == 50
    assert p4c.rename_network.call_count == 2


def test_silence_py4cytoscape(capsys):
    from py4cytoscape.exceptions import CyError

    CyError("handled error")  # silenced on import
    assert capsys.readouterr().err == ""

    try:
        cu.silence_py4cytoscape(False)
        e = CyError("shown error")
        assert "shown error" in capsys.readouterr().err
    finally:
        cu.silence_py4cytoscape(True)

    e = CyError("silenced again")
    assert capsys.readouterr().err == "" and "silenced again" in str(e)


def test_set_style_makes_network_current_first(p4c):
    cu.set_style("SKM", 7)
    names = [c[0] for c in p4c.mock_calls]
    assert names.index("set_current_network") < names.index("set_visual_style")
    p4c.set_current_network.assert_called_once_with(7)
    p4c.set_visual_style.assert_called_once_with("SKM", network=7)


# ---------------------------------------------------------------------------
# Bundled styles
# ---------------------------------------------------------------------------

def _bundled_styles():
    import xml.etree.ElementTree as ET
    from skm_tools import resources
    return {vs.get("name"): vs for vs in ET.parse(resources.get_style_xml_path()).getroot().iter("visualStyle")}


def test_bundled_xml_has_the_builtin_styles():
    from skm_tools import resources
    assert set(resources.BUILTIN_STYLES.values()) <= set(_bundled_styles())


def test_bundled_styles_map_node_type_to_pss_classes_only():
    old = {"gene", "protein_coding", "metabolite", "complex", "biotic", "abiotic", "mirna"}
    for name, vs in _bundled_styles().items():
        for vp in vs.iter("visualProperty"):
            for m in vp.iter("discreteMapping"):
                if m.get("attributeName") == "node_type":
                    values = {e.get("attributeValue") for e in m.iter("discreteMappingEntry")}
                    assert not values & old, (name, vp.get("name"))


@pytest.mark.parametrize("key, expected", [("skm", "SKM"), ("SKM-reactions", "SKM-reactions"),
                                           ("pss", "SKM"), ("ckn", "SKM")])
def test_apply_builtin_style_imports_once_and_applies(p4c, key, expected):
    p4c.styles.get_visual_style_names.return_value = ["default"]
    cu.apply_builtin_style(7, key)
    p4c.import_visual_styles.assert_called_once()
    p4c.set_visual_style.assert_called_once_with(expected, network=7)


def test_apply_builtin_style_does_not_reimport(p4c):
    p4c.styles.get_visual_style_names.return_value = ["default", "SKM"]
    cu.apply_builtin_style(7)
    p4c.import_visual_styles.assert_not_called()


def test_apply_builtin_style_unknown(p4c):
    with pytest.raises(ValueError, match="skm-reactions"):
        cu.apply_builtin_style(7, "fancy")
