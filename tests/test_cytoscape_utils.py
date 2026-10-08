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
    edges = cu.highlight_edges([("A", "B"), ("B", "C"), ("C", "A"), ("B", "A")], "#FF0000",
                               skip_edges=[("C", "A")], network=7)
    assert edges == [("A", "B"), ("B", "C")]  # found, not skipped
    suids = p4c.style_bypasses.set_edge_color_bypass.call_args.args[0]
    assert sorted(suids) == [11, 12, 13]


def test_highlight_edges_undirected(p4c):
    assert cu.highlight_edges([("B", "A")], "#FF0000") == []
    assert cu.highlight_edges([("B", "A")], "#FF0000", directed=False) == [("A", "B")]


def test_highlight_path_returns_highlighted_nodes_and_edges(p4c):
    nodes, edges = cu.highlight_path(["A", "B", "C", "nope"], "#00FF00", skip_nodes=["A"])
    assert nodes == ["B", "C"]
    assert edges == [("A", "B"), ("B", "C")]
    assert p4c.style_bypasses.set_node_color_bypass.call_args.args[0] == [2, 3]


def test_highlight_path_reversed_step(p4c):
    # a path from get_paths(..., directed=False): A - C uses the edge C -> A
    _, edges = cu.highlight_path(["A", "C"], "#00FF00", directed=False)
    assert edges == [("C", "A")]


def test_highlight_path_skip_workaround(p4c):
    # as in the PSS tutorial: later paths don't paint over earlier ones
    done_nodes, done_edges = [], []
    for path, color in [(["A", "B"], "#FF0000"), (["A", "B", "C"], "#0000FF")]:
        nodes, edges = cu.highlight_path(path, color, skip_nodes=done_nodes, skip_edges=done_edges)
        done_nodes += nodes
        done_edges += edges
    assert (nodes, edges) == (["C"], [("B", "C")])
    assert done_nodes == ["A", "B", "C"]


def test_highlight_nodes_only_sets_given_properties(p4c, caplog):
    assert cu.highlight_nodes(["C", "nope"], color="#000000", border_width=0) == ["C"]
    p4c.style_bypasses.set_node_color_bypass.assert_called_once()
    p4c.style_bypasses.set_node_border_width_bypass.assert_called_once_with([3], 0, network=None)
    p4c.style_bypasses.set_node_label_color_bypass.assert_not_called()
    assert "1 of 2 nodes not in the Cytoscape network: nope" in caplog.text


def test_clear_highlights(p4c):
    cu.clear_highlights(nodes=["A"], edges=[("A", "B")], network=7)
    node_calls = p4c.style_bypasses.clear_node_property_bypass.call_args_list
    assert {c.args[1] for c in node_calls} == set(cu._NODE_BYPASSES)
    assert all(c.args[0] == [1] for c in node_calls)
    edge_calls = p4c.style_bypasses.clear_edge_property_bypass.call_args_list
    assert {c.args[1] for c in edge_calls} == set(cu._EDGE_BYPASSES)
    assert all(sorted(c.args[0]) == [11, 12] for c in edge_calls)


def test_clear_highlights_all(p4c):
    cu.clear_highlights()
    assert p4c.style_bypasses.clear_node_property_bypass.call_args.args[0] == [1, 2, 3, 4, 5]


def test_apply_shortest_paths_style(p4c):
    _style_names(p4c, ["SKM"])
    g = nx.DiGraph([("A", "B"), ("B", "C"), ("C", "A")])

    name = cu.apply_shortest_paths_style(g, ["A", "C"], [[["A", "B", "C"]], [["C", "A"]]],
                                         edge_colors=["#FF0000", "#0000FF"], node_colors=["#000000"])

    assert name == "SKM-shortest-paths"
    p4c.copy_visual_style.assert_called_once_with("SKM", "SKM-shortest-paths")
    nodes = p4c.tables.load_table_data.call_args_list[0].args[0].set_index("name")
    # steps to the end of the path; on several paths, the smallest (A: 2 on A -> B -> C, 0 on C -> A)
    assert nodes["distance-to-target"].to_dict() == {"A": 0, "B": 1, "C": 0}
    assert nodes["node-path-source"].to_dict() == {"A": "A", "B": "A", "C": "A"}
    colors = p4c.style_mappings.set_node_color_mapping.call_args.kwargs["colors"]
    assert colors == ["#000000", "#FFFFFF"]  # one colour: to white


def test_apply_shortest_paths_style_rerun_does_not_nest(p4c):
    p4c.styles.get_current_style.return_value = "SKM-shortest-paths"
    g = nx.DiGraph([("A", "B")])
    assert cu.apply_shortest_paths_style(g, ["A"], [[["A", "B"]]]) == "SKM-shortest-paths"
    p4c.copy_visual_style.assert_not_called()


def test_subnetwork_edge_induced_from_paths(p4c):
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1")
    g.add_edge("A", "B", key="rx2")
    g.add_edge("B", "C", key="rx3")
    cu.subnetwork_edge_induced_from_paths(g, [["A", "B", "C"]], parent_suid=99, name="paths")
    kwargs = p4c.networks.create_subnetwork.call_args.kwargs
    assert sorted(kwargs["nodes"]) == [1, 2, 3]
    assert sorted(kwargs["edges"]) == [11, 12, 13]
    assert kwargs["exclude_edges"] is True
    assert kwargs["network"] == 99


def test_subnetwork_node_induced(p4c):
    cu.subnetwork_node_induced(["A", "C", "not-there"], parent_suid=99)
    assert p4c.networks.create_subnetwork.call_args.kwargs["nodes"] == [1, 3]


def test_subnetwork_neighbours_uses_networkx(p4c):
    g = nx.DiGraph([("A", "B"), ("B", "C"), ("X", "A")])
    cu.subnetwork_neighbours(g, "A", parent_suid=99, direction="out")
    assert p4c.networks.create_subnetwork.call_args.kwargs["nodes"] == [1, 2]
    p4c.select_nodes.assert_not_called()


# ---------------------------------------------------------------------------
# Networks, collections and styles
# ---------------------------------------------------------------------------

def test_clone_network_renames_collection_then_network(p4c):
    p4c.networks.clone_network.return_value = 50
    p4c.collections.get_collection_suid.return_value = 49

    assert cu.clone_network(name="Heat", collection="Collection - Heat", network=7) == 50

    op = p4c.commands.cyrest_put.call_args
    assert op.args[0] == "collections/49/tables/default"
    assert op.kwargs["body"]["data"] == [{"SUID": 49, "name": "Collection - Heat"}]
    p4c.rename_network.assert_called_once_with("Heat", network=50)
    names = [c[0] for c in p4c.mock_calls]
    assert names.index("commands.cyrest_put") < names.index("rename_network")


def test_clone_network_without_names(p4c):
    p4c.networks.clone_network.return_value = 50
    cu.clone_network(network=7)
    p4c.networks.clone_network.assert_called_once_with(network=7)
    p4c.commands.cyrest_put.assert_not_called()
    p4c.rename_network.assert_not_called()


def _style_names(p4c, names):
    """get_visual_style_names, with copy_visual_style and delete_visual_style changing it"""
    names = list(names)
    p4c.styles.get_visual_style_names.side_effect = lambda: list(names)
    p4c.copy_visual_style.side_effect = lambda style, new: names.append(new)
    p4c.styles.delete_visual_style.side_effect = names.remove
    return names


def test_copy_style_applies_to_networks(p4c):
    _style_names(p4c, ["SKM"])
    assert cu.copy_style("SKM", "SKM-heat", networks=[1, 2]) == "SKM-heat"
    p4c.copy_visual_style.assert_called_once_with("SKM", "SKM-heat")
    assert [c.kwargs["network"] for c in p4c.set_visual_style.call_args_list] == [1, 2]


def test_copy_style_existing_name(p4c):
    _style_names(p4c, ["SKM", "SKM-heat"])
    with pytest.raises(ValueError, match="overwrite=True"):
        cu.copy_style("SKM", "SKM-heat")
    assert cu.copy_style("SKM", "SKM-heat", overwrite=True) == "SKM-heat"
    p4c.styles.delete_visual_style.assert_called_once_with("SKM-heat")


def test_copy_style_name_not_as_asked(p4c):
    names = _style_names(p4c, ["SKM"])
    p4c.copy_visual_style.side_effect = lambda style, new: names.append(new + "_1")
    with pytest.raises(RuntimeError):
        cu.copy_style("SKM", "SKM-heat")


def test_delete_other_networks(p4c):
    p4c.get_network_list.return_value = [{"name": "a", "suid": 1}, {"name": "b", "suid": 2},
                                         {"name": "c", "suid": 3}]
    assert cu.delete_other_networks(2) == [1, 3]
    assert [c.args[0] for c in p4c.delete_network.call_args_list] == [1, 3]
    p4c.set_current_network.assert_called_once_with(2)


def test_get_or_create_subnetwork(p4c):
    # only the networks of the parent's collection
    p4c.get_collection_networks.return_value = [1, 5]
    p4c.get_network_name.side_effect = {1: "parent", 5: "existing"}.get
    p4c.networks.create_subnetwork.return_value = 6

    assert cu.get_or_create_subnetwork(["A"], 1, "existing") == 5
    p4c.networks.create_subnetwork.assert_not_called()

    assert cu.get_or_create_subnetwork(["A", "B"], 1, "new") == 6
    assert sorted(p4c.networks.create_subnetwork.call_args.kwargs["nodes"]) == [1, 2]


# ---------------------------------------------------------------------------
# Node images and charts
# ---------------------------------------------------------------------------

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

    p4c.styles.get_current_style.return_value = "my style"

    style = cu.add_custom_png(lambda n: tmp_path / f"{n}.png" if n == "A" else None, network=7,
                              position="left")

    assert style == "my style"
    assert list(p4c.load_table_data.call_args.args[0].index) == ["A"]
    p4c.style_dependencies.sync_node_custom_graphics_size.assert_called_once_with(False, style_name="my style")
    positions = [c.args[0]["value"] for c in p4c.style_defaults.set_visual_property_default.call_args_list
                 if c.args[0]["visualProperty"] == "NODE_CUSTOMGRAPHICS_POSITION_1"]
    assert positions == ["W,E,c,0.00,0.00"]


def test_add_custom_png_copies_bundled_style(p4c, tmp_path):
    _style_names(p4c, ["SKM"])
    p4c.get_all_nodes.return_value = ["A"]

    style = cu.add_custom_png(lambda n: tmp_path / "A.png", column="heat", network=7)

    assert style == "SKM-heat"
    p4c.copy_visual_style.assert_called_once_with("SKM", "SKM-heat")
    p4c.set_visual_style.assert_called_with("SKM-heat", network=7)
    p4c.style_dependencies.sync_node_custom_graphics_size.assert_called_once_with(False, style_name="SKM-heat")


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

def test_layout_from_coords(p4c):
    p4c.styles.get_visual_style_names.return_value = ["SKM"]
    cu.layout_from_coords({"A": (1.0, 2.0), "B": (-1.0, 0.5)}, scale=10, network=7)

    table = p4c.load_table_data.call_args.args[0]
    x, y = table.columns
    assert x.startswith("skm_tools_x_") and y.startswith("skm_tools_y_")
    assert table.loc["A"].tolist() == [10.0, -20.0]  # scaled, y flipped
    # the temporary style and columns are removed, the style set back
    tmp_style = p4c.copy_visual_style.call_args.args[1]
    assert p4c.set_visual_style.call_args.args[0] == "SKM"
    assert [c.args[0] for c in p4c.tables.delete_table_column.call_args_list] == [x, y]
    assert tmp_style.startswith("skm-tools-layout-")


def test_layout_from_coords_cleans_up_on_errors(p4c):
    p4c.styles.get_visual_style_names.side_effect = lambda: ["SKM", p4c.copy_visual_style.call_args.args[1]]
    p4c.update_style_mapping.side_effect = RuntimeError("Cytoscape error")
    table = pd.DataFrame({"x": [1.0], "y": [2.0]}, index=["A"])
    with pytest.raises(RuntimeError):
        cu.layout_from_coords(table, flip_y=False)
    assert p4c.load_table_data.call_args.args[0].iloc[0].tolist() == [1.0, 2.0]
    p4c.delete_visual_style.assert_called_once()
    assert p4c.tables.delete_table_column.call_count == 2
    assert p4c.set_visual_style.call_args.args[0] == "SKM"


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

def test_export_network_returns_the_written_file(p4c, tmp_path):
    p4c.network_views.export_image.side_effect = lambda filename, type, **kw: {"file": f"{filename}.pdf"}
    with patch.object(cu.time, "sleep") as sleep:
        path = cu.export_network(tmp_path / "net", network=7, wait=2)
    sleep.assert_called_once_with(2)
    assert path == tmp_path / "net.pdf"
    p4c.network_views.fit_content.assert_called_once_with(network=7)


def test_export_collection(p4c, tmp_path):
    p4c.get_collection_networks.return_value = [12, 11]
    p4c.get_network_name.side_effect = lambda suid: {11: "JA - Heat", 12: "SA/Heat"}[suid]
    # py4cytoscape adds the extension of the format
    p4c.network_views.export_image.side_effect = lambda filename, type, **kw: {"file": f"{filename}.jpeg"}

    with patch.object(cu.time, "sleep") as sleep:
        files = cu.export_collection(tmp_path / "out", format="JPEG", network=11, zoom=300)
    assert sleep.call_count == 2  # waits for Cytoscape before each export

    assert [f.name for f in files] == ["JA_Heat_11.jpeg", "SA_Heat_12.jpeg"]
    call = p4c.network_views.export_image.call_args.kwargs
    assert call["zoom"] == 300 and call["all_graphics_details"] is True and call["type"] == "JPEG"


def test_export_collection_crop_only_pdf(p4c, tmp_path):
    with pytest.raises(ValueError, match="PDF"):
        cu.export_collection(tmp_path, format="PNG", crop=True)


def test_clone_network_retries_rename(p4c):
    import py4cytoscape
    p4c.CyError = py4cytoscape.CyError
    p4c.networks.clone_network.return_value = 50
    p4c.rename_network.side_effect = [py4cytoscape.CyError("unrecognized (table entry)"), None]

    with patch.object(cu, "_RETRY_WAIT", 0):
        assert cu.clone_network(name="Heat", network=7) == 50
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


def test_bundled_style_mappings_use_the_loaders_attributes(pss_interaction_network_edge_path,
                                                          pss_interaction_network_node_path,
                                                          pss_reaction_graph_edge_path,
                                                          pss_reaction_graph_node_path,
                                                          ckn_edge_path, ckn_node_path):
    from skm_tools.ckn import ckn_to_networkx
    from skm_tools.pss import pss_interaction_network_to_networkx, pss_reaction_graph_to_networkx

    cytoscape_types = {"string": str, "boolean": bool, "integer": int}
    networks = {
        "SKM": [pss_interaction_network_to_networkx(pss_interaction_network_edge_path,
                                                    pss_interaction_network_node_path),
                ckn_to_networkx(ckn_edge_path, ckn_node_path)],
        "SKM-reactions": [pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path,
                                                         pss_reaction_graph_node_path)],
    }
    for name, vs in _bundled_styles().items():
        for section in ("node", "edge"):
            for m in vs.find(section).iter():
                if not m.tag.endswith("Mapping"):
                    continue
                attribute, expected = m.get("attributeName"), cytoscape_types[m.get("attributeType")]
                for g in networks[name]:
                    items = g.nodes(data=True) if section == "node" else g.edges(data=True)
                    values = [d[-1].get(attribute) for d in items]
                    assert any(v is not None for v in values), (name, attribute)
                    assert all(type(v) is expected for v in values if v is not None), (name, attribute)


def test_bundled_styles_map_node_type_to_pss_classes_only():
    old = {"gene", "protein_coding", "metabolite", "complex", "biotic", "abiotic", "mirna"}
    for name, vs in _bundled_styles().items():
        for vp in vs.iter("visualProperty"):
            for m in vp.iter("discreteMapping"):
                if m.get("attributeName") == "node_type":
                    values = {e.get("attributeValue") for e in m.iter("discreteMappingEntry")}
                    assert not values & old, (name, vp.get("name"))


@pytest.mark.parametrize("key, expected", [("skm", "SKM"), ("SKM-reactions", "SKM-reactions")])
def test_apply_builtin_style_imports_once_and_applies(p4c, key, expected):
    p4c.styles.get_visual_style_names.return_value = ["default"]
    assert cu.apply_builtin_style(key, network=7) == expected
    p4c.import_visual_styles.assert_called_once()
    p4c.set_visual_style.assert_called_once_with(expected, network=7)


def test_apply_builtin_style_does_not_reimport(p4c):
    p4c.styles.get_visual_style_names.return_value = ["default", "SKM"]
    cu.apply_builtin_style(network=7)
    p4c.import_visual_styles.assert_not_called()


def test_apply_builtin_style_unknown(p4c):
    for style in ("fancy", "pss", "ckn"):
        with pytest.raises(ValueError, match="skm-reactions"):
            cu.apply_builtin_style(style, network=7)


def test_load_network(p4c):
    p4c.networks.create_network_from_networkx.return_value = 7
    p4c.styles.get_visual_style_names.return_value = ["SKM"]
    g = nx.DiGraph([("A", "B")])
    assert cu.load_network(g, "PSS", collection="Tutorial", style="skm") == 7
    p4c.networks.create_network_from_networkx.assert_called_once_with(g, title="PSS", collection="Tutorial")
    p4c.set_visual_style.assert_called_once_with("SKM", network=7)


def test_load_network_without_style(p4c):
    cu.load_network(nx.DiGraph([("A", "B")]), "PSS")
    p4c.set_visual_style.assert_not_called()
