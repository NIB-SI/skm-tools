from pathlib import Path
import pytest

FIXTURES = Path(__file__).parent / "fixtures"

@pytest.fixture
def ckn_edge_path():
    return FIXTURES / "ckn_v2.0.1_edges.tsv"

@pytest.fixture
def ckn_node_path():
    return FIXTURES / "ckn_v2.0.1_nodes.tsv.gz"

@pytest.fixture
def pss_reaction_graph_edge_path():
    return FIXTURES / "pss_reaction_graph_edges.tsv"

@pytest.fixture
def pss_reaction_graph_node_path():
    return FIXTURES / "pss_reaction_graph_nodes.tsv"

@pytest.fixture
def pss_interaction_network_edge_path():
    return FIXTURES / "pss_interaction_network_edges.tsv"

@pytest.fixture
def pss_interaction_network_node_path():
    return FIXTURES / "pss_interaction_network_nodes.tsv"

@pytest.fixture
def pss_gene_network_ath_edge_path():
    return FIXTURES / "pss_gene_network_ath_edges.tsv"

@pytest.fixture
def pss_gene_network_ath_node_path():
    return FIXTURES / "pss_gene_network_ath_nodes.tsv"
