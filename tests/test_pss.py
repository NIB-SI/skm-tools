"""Tests for skm_tools.pss: loaders (small fixtures cut from the export samples) and utils."""

import networkx as nx
import pytest

from skm_tools.pss import (
    filter_pss_nodes,
    pss_gene_network_to_networkx,
    pss_interaction_network_to_networkx,
    pss_reaction_graph_to_networkx,
    remove_and_rewire,
    remove_duplicated_binding_edges,
    simplify_pss,
)


# ---------------------------------------------------------------------------
# PSS exports (pss-export): shared behaviour of all three loaders
# ---------------------------------------------------------------------------

PSS_LOADERS = {
    "reaction_graph": pss_reaction_graph_to_networkx,
    "interaction_network": pss_interaction_network_to_networkx,
    "gene_network_ath": pss_gene_network_to_networkx,
}


@pytest.fixture(params=list(PSS_LOADERS))
def pss_export(request):
    """(loader, edge path, node path) for each PSS export"""
    name = request.param
    return (
        PSS_LOADERS[name],
        request.getfixturevalue(f"pss_{name}_edge_path"),
        request.getfixturevalue(f"pss_{name}_node_path"),
    )


def _rows(path):
    with open(path) as handle:
        return [line.rstrip("\n").split("\t") for line in handle][1:]


class TestPSSExports:

    def test_returns_multidigraph(self, pss_export):
        loader, edge_path, node_path = pss_export
        assert isinstance(loader(edge_path, node_path), nx.MultiDiGraph)

    def test_all_nodes_and_edges_loaded(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        assert g.number_of_nodes() == len(_rows(node_path))
        assert g.number_of_edges() == len(_rows(edge_path))

    def test_accepts_string_paths(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(str(edge_path), str(node_path))
        assert g.number_of_edges() > 0

    def test_node_attributes_set(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, data in g.nodes(data=True):
            assert data["node_type"]
            assert isinstance(data["display_label"], str) and data["display_label"]
            assert "id" not in data

    def test_edge_endpoints_not_kept_as_attributes(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, _, data in g.edges(data=True):
            assert "source" not in data and "target" not in data

    def test_empty_values_are_none(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        values = [v for _, d in g.nodes(data=True) for v in d.values()]
        values += [v for _, _, d in g.edges(data=True) for v in d.values()]
        assert "" not in values
        assert None in values

    def test_list_columns_are_list_or_none(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, data in g.nodes(data=True):
            for k, v in data.items():
                if k in ("synonyms", "all_pathways", "external_links", "components") or k.endswith("_homologues"):
                    assert v is None or isinstance(v, list)

    def test_lists_split_on_semicolon_not_comma(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        synonyms = [s for _, d in g.nodes(data=True) for s in (d["synonyms"] or [])]
        assert len(synonyms) > g.number_of_nodes() / 2
        assert not any(";" in s for s in synonyms)
        # names with commas stay whole (10,11-EHT, AHK2,3,4, ...)
        assert any("," in s for s in synonyms)

    def test_complex_components_are_lists(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        complexes = [d for _, d in g.nodes(data=True) if d["node_type"] == "Complex"]
        assert complexes
        assert all(isinstance(d["components"], list) for d in complexes)

    def test_booleans_parsed(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        for _, _, data in g.edges(data=True):
            assert isinstance(data["directed"], bool)

    def test_missing_file_raises(self, pss_export, tmp_path):
        loader, edge_path, node_path = pss_export
        with pytest.raises(FileNotFoundError):
            loader(tmp_path / "missing-edges.tsv", node_path)


# ---------------------------------------------------------------------------
# PSS exports: per-format specifics
# ---------------------------------------------------------------------------

class TestPSSReactionGraph:

    def test_reactions_are_nodes(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        reactions = [n for n, d in g.nodes(data=True) if d["node_type"] == "reaction"]
        assert reactions
        for r in reactions:
            assert g.nodes[r]["reaction_type"]

    def test_edges_keyed_by_role(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["role"]

    def test_same_entity_twice_in_one_reaction(self, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
        # WRKY33 is both template and stimulator of its own expression (rx00160)
        g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
        assert set(g["WRKY33[fc00166]"]["rx00160"]) == {"template", "stimulator"}


class TestPSSInteractionNetwork:

    def test_edges_keyed_by_reaction_id(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["reaction_id"]

    def test_mutual_edges_listed_both_ways(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        mutual = [(u, v, k) for u, v, k, d in g.edges(keys=True, data=True) if not d["directed"]]
        assert mutual
        for u, v, k in mutual:
            assert g.has_edge(v, u, key=k)

    def test_homologues_loaded(self, pss_interaction_network_edge_path, pss_interaction_network_node_path):
        g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
        plant = [d for _, d in g.nodes(data=True) if d["node_type"] == "PlantCoding"]
        assert plant
        assert any(d["ath_homologues"] for d in plant)


class TestPSSGeneNetwork:

    def test_genes_have_species_and_entity(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        genes = [d for _, d in g.nodes(data=True) if d["node_type"] == "gene"]
        assert genes
        for d in genes:
            assert d["species"] == "ath"
            assert d["entity"]

    def test_cluster_columns_are_aligned_lists(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        multi = [d for _, d in g.nodes(data=True) if len(d["functional_cluster_id"] or []) > 1]
        assert multi
        for _, d in g.nodes(data=True):
            n = len(d["entity"])
            for c in ("short_name", "pathway", "functional_cluster_id"):
                assert d[c] is None or len(d[c]) == n

    def test_edges_keyed_by_reaction_id(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        for _, _, key, data in g.edges(keys=True, data=True):
            assert key == data["reaction_id"]


# ---------------------------------------------------------------------------
# Filtering, simplifying and rewiring
# ---------------------------------------------------------------------------

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
