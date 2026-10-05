import networkx as nx

import pytest

from skm_tools.load_networks import pss_gene_network_to_networkx, pss_interaction_network_to_networkx
from skm_tools.pss_utils import filter_pss_nodes, simplify_pss, remove_and_rewire, remove_duplicated_binding_edges


def _base_edge(reaction_id, reaction_effect):
    return {
        "reaction_id": reaction_id,
        "interaction": "positive-influence" if reaction_effect == "activation" else "negative-influence",
        "directed": True,
        "reaction_type": "protein activation",
        "reaction_effect": reaction_effect,
        "source_role": "modifier",
        "source_location": "cytoplasm",
        "source_form": "protein",
        "target_role": "product",
        "target_location": "cytoplasm",
        "target_form": "protein_active",
    }


def _make_multigraph(reaction_effects):
    g = nx.MultiDiGraph()
    g.add_node("A", node_type="PlantCoding")
    g.add_node("B", node_type="PlantCoding")
    for i, effect in enumerate(reaction_effects):
        # keyed by reaction_id, as the loaders do
        g.add_edge("A", "B", key=f"R{i}", **_base_edge(f"R{i}", effect))
    return g


def test_simplify_pss_merges_parallel_edges_by_default():
    g = _make_multigraph(["activation", "inhibition"])

    new_g = simplify_pss(g)

    assert isinstance(new_g, nx.DiGraph)
    assert not isinstance(new_g, nx.MultiDiGraph)
    assert new_g.number_of_edges() == 1
    data = new_g["A"]["B"]
    assert set(data["reaction_id"].split(",")) == {"R0", "R1"}


def test_simplify_pss_keeps_edges_separate_when_split_attr_differs():
    g = _make_multigraph(["activation", "inhibition"])

    new_g = simplify_pss(g, split_on_attrs=["reaction_effect"])

    assert isinstance(new_g, nx.MultiDiGraph)
    assert new_g.number_of_edges("A", "B") == 2
    effects = {d["reaction_effect"] for d in new_g["A"]["B"].values()}
    assert effects == {"activation", "inhibition"}


def test_simplify_pss_keeps_attributes_not_in_the_original_whitelist():
    g = _make_multigraph(["activation", "inhibition"])
    g["A"]["B"]["R0"]["custom_score"] = "high"
    g["A"]["B"]["R1"]["custom_score"] = "high"

    new_g = simplify_pss(g)

    assert new_g["A"]["B"]["custom_score"] == "high"


def test_simplify_pss_merges_within_matching_split_groups():
    g = _make_multigraph(["activation", "activation", "inhibition"])

    new_g = simplify_pss(g, split_on_attrs=["reaction_effect"])

    assert new_g.number_of_edges("A", "B") == 2
    reaction_ids_by_effect = {
        d["reaction_effect"]: set(d["reaction_id"].split(","))
        for d in new_g["A"]["B"].values()
    }
    assert reaction_ids_by_effect["activation"] == {"R0", "R1"}
    assert reaction_ids_by_effect["inhibition"] == {"R2"}


def _binding_edge(reaction_id):
    edge = _base_edge(reaction_id, "inhibition")
    edge["reaction_type"] = "binding/oligomerisation"
    edge["interaction"] = "negative-influence"
    edge["directed"] = False
    edge["source_role"] = edge["target_role"] = "interactor"
    return edge


def test_remove_duplicated_binding_edges_on_digraph():
    g = nx.DiGraph()
    g.add_node("A", node_type="PlantCoding")
    g.add_node("B", node_type="PlantCoding")
    g.add_edge("A", "B", **_binding_edge("R0"))
    g.add_edge("B", "A", **_binding_edge("R0"))

    remove_duplicated_binding_edges(g)

    assert g.number_of_edges() == 1


def test_remove_duplicated_binding_edges_on_multidigraph_keeps_unrelated_parallel_edges():
    g = nx.MultiDiGraph()
    g.add_node("A", node_type="PlantCoding")
    g.add_node("B", node_type="PlantCoding")
    g.add_edge("A", "B", **_binding_edge("R0"))
    g.add_edge("A", "B", **_base_edge("R1", "activation"))  # unrelated parallel edge, not binding
    g.add_edge("B", "A", **_binding_edge("R0"))

    remove_duplicated_binding_edges(g)

    assert g.number_of_edges() == 2
    remaining_reaction_ids = {d["reaction_id"] for _, _, d in g.edges(data=True)}
    assert remaining_reaction_ids == {"R0", "R1"}


def test_simplify_pss_single_edge_with_non_integer_key():
    g = _make_multigraph(["activation"])

    new_g = simplify_pss(g)

    assert new_g["A"]["B"]["reaction_id"] == "R0"


def _chain():
    """A -> B -> C (DiGraph), for rewiring"""
    g = nx.DiGraph()
    for n in "ABC":
        g.add_node(n, node_type="PlantCoding")
    g.add_edge("A", "B", **_base_edge("R0", "activation"))
    g.add_edge("B", "C", **_base_edge("R1", "activation"))
    return g


def test_remove_and_rewire_rejects_multidigraph():
    with pytest.raises(NotImplementedError):
        remove_and_rewire(_make_multigraph(["activation"]), ["A"])


def test_remove_and_rewire_connects_upstream_to_downstream():
    g = _chain()

    remove_and_rewire(g, ["B"])

    assert "B" not in g
    assert set(g["A"]["C"]["reaction_id"].split(",")) == {"R0", "R1"}
    assert g["A"]["C"]["interaction"] == "positive-influence"


def test_remove_and_rewire_does_not_propagate_mutual_binding():
    g = _chain()
    g.add_node("D", node_type="PlantCoding")
    g.add_edge("D", "B", **_binding_edge("R2"))
    g.add_edge("B", "D", **_binding_edge("R2"))

    remove_and_rewire(g, ["B"])

    assert g.has_edge("A", "C")
    assert not g.has_edge("D", "C")


def test_filter_pss_nodes_removes_complexes_by_component():
    g = nx.MultiDiGraph()
    g.add_node("CTR[fc00049]", node_type="PlantCoding")
    g.add_node("ETR[fc00075]", node_type="PlantCoding")
    g.add_node("ET", node_type="Metabolite")
    g.add_node("CTR|ETR", node_type="Complex", components=["CTR[fc00049]", "ETR[fc00075]"])
    g.add_node("ET|ETR", node_type="Complex", components=["ET", "ETR[fc00075]"])

    reasons = filter_pss_nodes(g, node_types=["PlantCoding", "Complex"], remove_isolates=False)

    assert reasons["ET"] == "wrong node type"
    assert reasons["ET|ETR"] == "complex component removed"
    assert "CTR|ETR" in g


def test_filter_pss_nodes_ignores_components_not_in_network():
    # e.g. ETP|SCF lists CUL and RBX, which aren't nodes of the interaction network
    g = nx.MultiDiGraph()
    g.add_node("A", node_type="PlantCoding")
    g.add_node("A|X", node_type="Complex", components=["A", "X"])

    filter_pss_nodes(g, node_types=["PlantCoding", "Complex"], remove_isolates=False)

    assert "A|X" in g


def test_filter_pss_nodes_gene_network_keeps_complex_while_cluster_has_genes():
    g = nx.MultiDiGraph()
    g.add_node("G1", node_type="gene", entity=["ETR[fc00075]"])
    g.add_node("G2", node_type="gene", entity=["ETR[fc00075]", "OTHER[fc1]"])
    g.add_node("ET", node_type="Metabolite", entity=["ET"])
    g.add_node("ETR|X", node_type="Complex", entity=["ETR|X"], components=["ETR[fc00075]"])

    h = g.copy()
    filter_pss_nodes(h, node_types=["gene", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" in h

    h = g.copy()
    h.nodes["G1"]["node_type"] = "drop"
    filter_pss_nodes(h, node_types=["gene", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" in h  # G2 is still in ETR

    h = g.copy()
    h.nodes["G1"]["node_type"] = h.nodes["G2"]["node_type"] = "drop"
    filter_pss_nodes(h, node_types=["gene", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" not in h


def test_filter_pss_nodes_species_on_interaction_network(pss_interaction_network_edge_path, pss_interaction_network_node_path):
    g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)

    filter_pss_nodes(g, species=["ath"], remove_isolates=False)

    for _, d in g.nodes(data=True):
        if d["node_type"] in ("PlantCoding", "PlantNonCoding"):
            assert d["ath_homologues"]
