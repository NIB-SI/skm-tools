"""Tests for skm_tools.pss: loaders (small fixtures cut from the export samples) and utils."""

import shutil
from unittest.mock import patch

import networkx as nx
import pytest

from skm_tools.pss import (
    filter_pss_nodes,
    pss_gene_network_to_networkx,
    pss_interaction_network_to_networkx,
    pss_reaction_graph_to_networkx,
    remove_and_rewire,
    remove_deadend_complexes,
    remove_duplicated_binding_edges,
    remove_reactions,
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
                if k in ("synonyms", "all_pathways", "mapman", "external_links", "components",
                         "component_cluster_ids") or k.endswith("_homologues"):
                    assert v is None or isinstance(v, list)

    def test_lists_split_on_pipe_only(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        synonyms = [s for _, d in g.nodes(data=True) for s in (d["synonyms"] or [])]
        assert len(synonyms) > g.number_of_nodes() / 2
        assert not any("|" in s for s in synonyms)
        # names with commas stay whole (10,11-EHT, AHK2,3,4, ...)
        assert any("," in s for s in synonyms)

    def test_complex_ids_keep_pipe(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        assert any("|" in n for n, d in g.nodes(data=True) if d["node_type"] == "Complex")

    def test_rank_zero_on_edges(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        assert {d["rank"] for *_, d in g.edges(data=True)} == {0}

    def test_mapman_and_component_cluster_ids(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        assert any(d["mapman"] for _, d in g.nodes(data=True))
        complexes = [d for _, d in g.nodes(data=True) if d["node_type"] == "Complex" and d["component_cluster_ids"]]
        assert complexes
        for d in complexes:
            assert all(c.startswith("fc") for c in d["component_cluster_ids"])

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

    def test_missing_file_downloaded(self, pss_export, tmp_path):
        loader, edge_path, node_path = pss_export
        missing = tmp_path / "missing-edges.tsv"
        with patch("skm_tools.utils.download",
                   side_effect=lambda url, path: shutil.copy(edge_path, path)) as mock_dl:
            g = loader(missing, node_path)
        mock_dl.assert_called_once()
        url = mock_dl.call_args[0][0]
        assert url.startswith("https://skm.nib.si/downloads/pss/public/live/") and url.endswith("/edges")
        assert missing.exists()
        assert sorted(g.edges(keys=True)) == sorted(loader(edge_path, node_path).edges(keys=True))

    def test_network_recorded(self, pss_export):
        loader, edge_path, node_path = pss_export
        g = loader(edge_path, node_path)
        expected = {pss_reaction_graph_to_networkx: "reaction_graph",
                    pss_interaction_network_to_networkx: "interaction_network",
                    pss_gene_network_to_networkx: "gene_network"}[loader]
        assert g.graph["pss_network"] == expected


def test_gene_network_download_urls_and_default_names(
        pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path, tmp_path):
    sources = {"edges": pss_gene_network_ath_edge_path, "nodes": pss_gene_network_ath_node_path}
    with patch("skm_tools.utils.download",
               side_effect=lambda url, path: shutil.copy(sources[url.rsplit("/", 1)[1]], path)) as mock_dl:
        g = pss_gene_network_to_networkx(species="stu", data_dir=tmp_path)
    assert [c[0] for c in mock_dl.call_args_list] == [
        ("https://skm.nib.si/downloads/pss/public/live/gene-network/stu/edges",
         tmp_path / "pss-public-gene-network-stu-edges.tsv"),
        ("https://skm.nib.si/downloads/pss/public/live/gene-network/stu/nodes",
         tmp_path / "pss-public-gene-network-stu-nodes.tsv"),
    ]
    assert g.graph == {"pss_network": "gene_network", "species": "stu"}


def test_gene_network_download_defaults_to_ath(pss_gene_network_ath_node_path, tmp_path):
    with patch("skm_tools.utils.download") as mock_dl:
        mock_dl.side_effect = lambda url, path: shutil.copy(
            pss_gene_network_ath_node_path.with_name("pss_gene_network_ath_edges.tsv"), path)
        pss_gene_network_to_networkx(tmp_path / "e.tsv", pss_gene_network_ath_node_path)
    assert mock_dl.call_args[0][0] == "https://skm.nib.si/downloads/pss/public/live/gene-network/ath/edges"


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

    def test_genes_have_species_and_clusters(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        genes = [d for _, d in g.nodes(data=True) if d.get("species")]
        assert genes
        for d in genes:
            assert d["species"] == "ath"
            # the class of the gene's functional cluster
            assert d["node_type"] in ("PlantCoding", "PlantNonCoding")
            assert d["functional_cluster_id"]
            assert "entity" not in d

    def test_cluster_columns_are_aligned_lists(self, pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
        g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
        multi = [d for _, d in g.nodes(data=True) if len(d["functional_cluster_id"] or []) > 1]
        assert multi
        for _, d in g.nodes(data=True):
            if not d.get("species"):  # not a gene
                continue
            n = len(d["functional_cluster_id"])
            assert isinstance(d["display_label"], str)
            for c in ("short_name", "pathway"):
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
    assert data["reaction_id"] == ["R0", "R1"]


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
        d["reaction_effect"]: d["reaction_id"] for d in new_g["A"]["B"].values()
    }
    assert reaction_ids_by_effect["activation"] == ["R0", "R1"]
    assert reaction_ids_by_effect["inhibition"] == ["R2"]


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

    removed = remove_duplicated_binding_edges(g)

    assert g.number_of_edges() == 1
    assert len(removed) == 1 and not g.has_edge(*removed[0])


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

    assert new_g["A"]["B"]["reaction_id"] == ["R0"]  # a list for every edge
    assert g["A"]["B"]["R0"]["reaction_id"] == "R0"  # g unchanged


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

    reasons = remove_and_rewire(g, ["B"])

    assert reasons == {"B": "rewired"}
    assert "B" not in g
    data = g["A"]["C"]
    assert data["reaction_id"] == ["R0", "R1"]
    assert data["interaction"] == "positive-influence"
    assert data["directed"] is True
    assert data["note"] == "rewired through B"
    # roles of the ends: the source's from A -> B, the target's from B -> C
    assert (data["source_role"], data["target_role"]) == ("modifier", "product")


def test_remove_and_rewire_does_not_propagate_mutual_binding():
    g = _chain()
    g.add_node("D", node_type="PlantCoding")
    g.add_edge("D", "B", **_binding_edge("R2"))
    g.add_edge("B", "D", **_binding_edge("R2"))

    remove_and_rewire(g, ["B"])

    assert g.has_edge("A", "C")
    assert not g.has_edge("D", "C")  # upstream through binding
    assert not g.has_edge("A", "D")  # downstream through binding


def test_remove_and_rewire_through_complex_formation():
    g = _chain()
    g.add_node("B|D", node_type="Complex")
    forming = _binding_edge("R2")
    forming.update(directed=True, interaction="positive-influence", target_role="product")
    g.add_edge("B", "B|D", **forming)

    remove_and_rewire(g, ["B"])

    assert g.has_edge("A", "B|D")


def test_remove_and_rewire_merges_into_existing_edges():
    g = _chain()
    g.add_edge("A", "C", **_base_edge("R5", "activation"), custom="kept")

    remove_and_rewire(g, ["B"])

    assert g["A"]["C"]["reaction_id"] == ["R0", "R1", "R5"]
    assert g["A"]["C"]["custom"] == "kept"


def test_remove_and_rewire_two_removed_nodes_between_the_same_pair():
    g = nx.DiGraph()
    for u, v, r in [("A", "B1", "R0"), ("B1", "C", "R1"), ("A", "B2", "R2"), ("B2", "C", "R3")]:
        g.add_edge(u, v, **_base_edge(r, "activation"))

    remove_and_rewire(g, ["B1", "B2"])

    assert list(g.edges()) == [("A", "C")]
    assert g["A"]["C"]["reaction_id"] == ["R0", "R1", "R2", "R3"]


def test_remove_and_rewire_chain_of_removed_nodes():
    g = _chain()
    g.add_edge("C", "D", **_base_edge("R2", "inhibition"))

    reasons = remove_and_rewire(g, ["B", "C"])

    assert list(g.edges()) == [("A", "D")]
    assert g["A"]["D"]["reaction_id"] == ["R0", "R1", "R2"]
    assert g["A"]["D"]["interaction"] == "negative-influence"
    assert reasons == {"B": "rewired", "C": "rewired"}


@pytest.mark.parametrize("first, second, composed", [
    ("activation", "activation", "positive-influence"),
    ("activation", "inhibition", "negative-influence"),
    ("inhibition", "inhibition", "positive-influence"),
])
def test_remove_and_rewire_sign_of_the_chain(first, second, composed):
    g = nx.DiGraph()
    g.add_edge("A", "B", **_base_edge("R0", first))
    g.add_edge("B", "C", **_base_edge("R1", second))
    remove_and_rewire(g, "B")
    assert g["A"]["C"]["interaction"] == composed


def test_remove_and_rewire_unknown_sign():
    g = _chain()
    g["B"]["C"]["interaction"] = "unknown-influence"
    remove_and_rewire(g, ["B"])
    assert g["A"]["C"]["interaction"] == "unknown-influence"


def test_remove_and_rewire_dry_run():
    g = _chain()
    before = nx.to_dict_of_dicts(g)

    planned = remove_and_rewire(g, ["B"], dry_run=True)

    assert nx.to_dict_of_dicts(g) == before
    assert [(u, v) for u, v, _ in planned] == [("A", "C")]
    assert planned[0][2]["reaction_id"] == ["R0", "R1"]


def test_remove_and_rewire_nothing_to_rewire_and_missing_nodes(caplog):
    g = _chain()
    assert remove_and_rewire(g, ["A", "nope"]) == {"A": "nothing to rewire"}
    assert "nope" in caplog.text


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
    g.add_node("G1", node_type="PlantCoding", species="ath", functional_cluster_id=["fc00075"])
    g.add_node("G2", node_type="PlantCoding", species="ath", functional_cluster_id=["fc00075", "fc1"])
    g.add_node("ET", node_type="Metabolite")
    # components are interaction-network names: not node ids in a gene network
    g.add_node("ETR|X", node_type="Complex", components=["ETR[fc00075]", "X"],
               component_cluster_ids=["fc00075"])

    h = g.copy()
    filter_pss_nodes(h, node_types=["PlantCoding", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" in h

    h = g.copy()
    h.nodes["G1"]["node_type"] = "drop"
    filter_pss_nodes(h, node_types=["PlantCoding", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" in h  # G2 is still in ETR

    h = g.copy()
    h.nodes["G1"]["node_type"] = h.nodes["G2"]["node_type"] = "drop"
    filter_pss_nodes(h, node_types=["PlantCoding", "Complex", "Metabolite"], remove_isolates=False)
    assert "ETR|X" not in h


def test_filter_pss_nodes_species_on_interaction_network(pss_interaction_network_edge_path, pss_interaction_network_node_path):
    g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)

    filter_pss_nodes(g, species=["ath"], remove_isolates=False)

    for _, d in g.nodes(data=True):
        if d.get("functional_cluster_id") and d["node_type"] != "PlantAbstract":
            assert d["ath_homologues"]


def test_filter_pss_nodes_gene_network_non_cluster_component_by_id():
    g = nx.MultiDiGraph()
    g.add_node("G1", node_type="PlantCoding", species="ath", functional_cluster_id=["fc00075"])
    g.add_node("ET", node_type="Metabolite")
    g.add_node("ET|ETR", node_type="Complex", components=["ET", "ETR[fc00075]"],
               component_cluster_ids=["fc00075"])

    reasons = filter_pss_nodes(g, node_types=["PlantCoding", "Complex"], remove_isolates=False)

    assert reasons["ET|ETR"] == "complex component removed"
    assert "G1" in g


def _species_network():
    """rx1: ENZ (no ath genes) catalyses S -> P; rx2: K (ath genes) catalyses S -> X;
    rx3: ENZ and K bind to form ENZ|K, which (rx4) activates X; rx5: ABS (abstract,
    no genes) catalyses X -> Y."""
    g = nx.MultiDiGraph()
    g.add_node("ENZ[fc1]", node_type="PlantCoding", functional_cluster_id="fc1")
    g.add_node("K[fc2]", node_type="PlantCoding", functional_cluster_id="fc2", ath_homologues=["AT1G01010"])
    g.add_node("ABS[fc3]", node_type="PlantAbstract", functional_cluster_id="fc3")
    for n in ("S", "P", "X", "Y"):
        g.add_node(n, node_type="Metabolite")
    g.add_node("ENZ|K", node_type="Complex", components=["ENZ[fc1]", "K[fc2]"],
               component_cluster_ids=["fc1", "fc2"])
    for u, v, r in [("S", "P", "rx1"), ("ENZ[fc1]", "P", "rx1"), ("ENZ[fc1]", "S", "rx1"),
                    ("S", "X", "rx2"), ("K[fc2]", "X", "rx2"),
                    ("ENZ[fc1]", "ENZ|K", "rx3"), ("K[fc2]", "ENZ|K", "rx3"),
                    ("ENZ|K", "X", "rx4"),
                    ("X", "Y", "rx5"), ("ABS[fc3]", "Y", "rx5")]:
        g.add_edge(u, v, key=r, reaction_id=r)
    return g


def test_filter_pss_nodes_species_removes_whole_reactions():
    g = _species_network()

    reasons = filter_pss_nodes(g, species=["ath"])

    # rx1 and rx3 involve ENZ, which has no ath genes: all their edges go, also S -> P;
    # rx4 involves ENZ|K, a complex with ENZ as component
    assert {k for *_, k in g.edges(keys=True)} == {"rx2", "rx5"}
    assert reasons["ENZ[fc1]"] == "species missing"
    assert reasons["ENZ|K"] == "complex component removed"
    assert reasons["P"] == "isolate"
    # abstract clusters have no genes, and don't decide
    assert "ABS[fc3]" in g


def test_filter_pss_nodes_species_on_merged_edges_keeps_other_reactions():
    g = simplify_pss(_species_network())
    g.add_edge("S", "P", reaction_id=["rx1", "rx2"])  # S -> P also through a kept reaction

    filter_pss_nodes(g, species=["ath"])

    assert g["S"]["P"]["reaction_id"] == ["rx2"]


def test_filter_pss_nodes_string_arguments():
    g, h = _species_network(), _species_network()
    assert filter_pss_nodes(g, species="ath", node_types="Metabolite") == \
        filter_pss_nodes(h, species=["ath"], node_types=["Metabolite"])
    assert set(g.edges(keys=True)) == set(h.edges(keys=True))


def test_filter_pss_nodes_species_matches_gene_network(
        pss_interaction_network_edge_path, pss_interaction_network_node_path,
        pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
    g = pss_interaction_network_to_networkx(pss_interaction_network_edge_path, pss_interaction_network_node_path)
    gn = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)

    filter_pss_nodes(g, species=["ath"])

    assert set(g.edges(keys=True)) == {
        (d["source_entity"], d["target_entity"], k) for _, _, k, d in gn.edges(keys=True, data=True)
    }


@pytest.mark.parametrize("func", [
    lambda g: filter_pss_nodes(g, node_types=["Metabolite"]),
    simplify_pss,
    lambda g: remove_and_rewire(simplify_pss(g), []),
    remove_duplicated_binding_edges,
])
def test_functions_reject_reaction_graph(func, pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
    g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
    with pytest.raises(ValueError, match="reaction graph"):
        func(g)


def test_filter_pss_nodes_species_rejected_on_gene_network(
        pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
    g = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)
    n = g.number_of_nodes()
    with pytest.raises(ValueError, match="one species"):
        filter_pss_nodes(g, species=["stu"])
    assert g.number_of_nodes() == n


# ---------------------------------------------------------------------------
# remove_reactions, and the reaction graph
# ---------------------------------------------------------------------------

def _reactions(g):
    return {n for n, d in g.nodes(data=True) if d["node_type"] == "reaction"}


def test_remove_reactions_interaction_network():
    g = _species_network()

    reasons = remove_reactions(g, ["rx1", "rx99"])

    assert {k for *_, k in g.edges(keys=True)} == {"rx2", "rx3", "rx4", "rx5"}
    assert reasons == {"P": "isolate"}


def test_remove_reactions_keeps_isolates_if_asked():
    g = _species_network()
    remove_reactions(g, ["rx1"], remove_isolates=False)
    assert "P" in g


def test_remove_reactions_merged_edges_keep_other_reactions():
    g = nx.DiGraph([("S", "P", {"reaction_id": ["rx1", "rx2"]}), ("P", "X", {"reaction_id": "rx1"})])

    remove_reactions(g, ["rx1"])

    assert list(g.edges(data="reaction_id")) == [("S", "P", ["rx2"])]


def test_remove_reactions_single_string():
    g = _species_network()
    remove_reactions(g, "rx1")
    assert {k for *_, k in g.edges(keys=True)} == {"rx2", "rx3", "rx4", "rx5"}


def test_remove_reactions_reaction_graph(pss_reaction_graph_edge_path, pss_reaction_graph_node_path):
    g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
    participants = set(nx.all_neighbors(g, "rx00001"))

    reasons = remove_reactions(g, ["rx00001"])

    assert "rx00001" not in g and reasons["rx00001"] == "reaction removed"
    assert all(n in g or reasons[n] == "isolate" for n in participants)
    assert all(g.degree(n) > 0 for n in g)


def test_filter_pss_nodes_species_on_reaction_graph_matches_gene_network(
        pss_reaction_graph_edge_path, pss_reaction_graph_node_path,
        pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path):
    g = pss_reaction_graph_to_networkx(pss_reaction_graph_edge_path, pss_reaction_graph_node_path)
    gn = pss_gene_network_to_networkx(pss_gene_network_ath_edge_path, pss_gene_network_ath_node_path)

    filter_pss_nodes(g, species=["ath"])

    # translocations without a transporter give no influence edges, so no gene network edges
    no_edges = {r for r in _reactions(g)
                if g.nodes[r]["reaction_type"] == "translocation"
                and not any(k == "transporter" for _, _, k in g.in_edges(r, keys=True))}
    assert _reactions(g) - no_edges == {k for *_, k in gn.edges(keys=True)}


def test_remove_deadend_complexes_interaction_network_keeps_partner_edges():
    g = nx.MultiDiGraph()
    g.add_nodes_from(["A", "B"], node_type="PlantCoding")
    g.add_node("A|B", node_type="Complex")
    for u, v in [("A", "B"), ("B", "A"), ("A", "A|B"), ("B", "A|B")]:
        g.add_edge(u, v, key="rx1", reaction_id="rx1")

    assert remove_deadend_complexes(g) == {"A|B": "dead-end complex"}
    assert set(g.edges()) == {("A", "B"), ("B", "A")}


def test_remove_deadend_complexes_reaction_graph_removes_forming_reactions():
    g = nx.MultiDiGraph()
    g.add_nodes_from(["A", "B", "C"], node_type="PlantCoding")
    g.add_nodes_from(["A|B", "A|B|C"], node_type="Complex")
    g.add_nodes_from(["rx1", "rx2", "rx3"], node_type="reaction")
    # rx1: A + B -> A|B; rx2: A|B + C -> A|B|C (a dead end, so A|B becomes one too);
    # rx3: C -> B
    for u, v, role in [("A", "rx1", "interactor"), ("B", "rx1", "interactor"), ("rx1", "A|B", "product"),
                       ("A|B", "rx2", "interactor"), ("C", "rx2", "interactor"), ("rx2", "A|B|C", "product"),
                       ("C", "rx3", "modifier"), ("rx3", "B", "product")]:
        g.add_edge(u, v, key=role, role=role)

    reasons = remove_deadend_complexes(g)

    assert sorted(n for n, r in reasons.items() if r == "dead-end complex") == ["A|B", "A|B|C"]
    assert reasons["rx1"] == reasons["rx2"] == "reaction removed"
    assert set(g.nodes()) == {"B", "C", "rx3"}


def test_remove_deadend_complexes_node_without_node_type():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1", reaction_id="rx1")
    assert remove_deadend_complexes(g) == {}


def _conflicting_edges():
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1", reaction_id="rx1", interaction="positive-influence", reaction_type="catalysis")
    g.add_edge("A", "B", key="rx2", reaction_id="rx2", interaction="negative-influence", reaction_type="catalysis")
    g.add_edge("B", "C", key="rx3", reaction_id="rx3", interaction="positive-influence", reaction_type="binding")
    g.add_edge("B", "C", key="rx4", reaction_id="rx4", interaction="negative-influence", reaction_type="catalysis")
    return g


def test_simplify_pss_reports_a_summary(caplog):
    with caplog.at_level("INFO", logger="skm_tools"):
        simplify_pss(_conflicting_edges())
    assert len(caplog.records) == 1
    assert "interaction (2)" in caplog.text and "reaction_type (1)" in caplog.text


def test_simplify_pss_verbose_reports_every_edge(caplog):
    with caplog.at_level("INFO", logger="skm_tools"):
        simplify_pss(_conflicting_edges(), verbose=True)
    assert "['rx1', 'rx2'] --> interaction" in caplog.text
    assert "['rx3', 'rx4'] --> reaction_type" in caplog.text


def test_simplify_pss_no_report_without_differences(caplog):
    g = nx.MultiDiGraph()
    g.add_edge("A", "B", key="rx1", reaction_id="rx1", interaction="positive-influence")
    g.add_edge("A", "B", key="rx2", reaction_id="rx2", interaction="positive-influence")
    with caplog.at_level("INFO", logger="skm_tools"):
        simplify_pss(g)
    assert caplog.text == ""


def test_simplify_pss_merge_rule():
    g = nx.MultiDiGraph(pss_network="interaction_network")
    # added out of reaction id order: the kept values don't depend on the order
    g.add_edge("A", "B", key="rx2", reaction_id="rx2", interaction="negative-influence", directed=True,
               rank=2, reaction_type="binding", forms=["x"])
    g.add_edge("A", "B", key="rx1", reaction_id="rx1", interaction="positive-influence", directed=False,
               rank=0, reaction_type="catalysis", forms=["x"])
    h = simplify_pss(g)
    assert h["A"]["B"] == {"reaction_id": ["rx1", "rx2"], "interaction": "unknown-influence",
                           "directed": True, "rank": 0, "reaction_type": "catalysis", "forms": ["x"]}
    assert h.graph == {"pss_network": "interaction_network"}


def test_simplify_pss_rejects_simple_graph():
    with pytest.raises(ValueError, match="already simplified"):
        simplify_pss(simplify_pss(_conflicting_edges()))
