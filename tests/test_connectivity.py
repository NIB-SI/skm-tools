import networkx as nx
import numpy as np
import pytest

from skm_tools.connectivity import _components_largest, induced_connectivity


def chain_with_clique(n_chain=60, clique=6):
    '''A clique of `clique` nodes (0..clique-1) attached to a long path.'''
    g = nx.complete_graph(clique)
    g.add_edges_from((i, i + 1) for i in range(clique - 1, clique + n_chain))
    return g


class TestComponents:
    def test_largest_component(self):
        u, v = np.array([0, 1, 3]), np.array([1, 2, 4])
        assert _components_largest(6, u, v) == 3

    def test_no_edges(self):
        assert _components_largest(4, np.array([], int), np.array([], int)) == 1

    def test_long_path_converges(self):
        n = 200
        assert _components_largest(n, np.arange(n - 1), np.arange(1, n)) == n


class TestInducedConnectivity:
    def test_clique_is_more_connected(self):
        g = chain_with_clique()
        r = induced_connectivity(g, range(6), n_random=199, seed=1)
        assert r["n_nodes"] == 6
        assert r["edges"] == 15 and r["lcc"] == 6
        assert r["edges_p"] == pytest.approx(1 / 200)
        assert r["edges_fold"] > 3

    def test_scattered_nodes_are_not(self):
        g = chain_with_clique()
        r = induced_connectivity(g, [10, 20, 30, 40, 50, 60], n_random=199, seed=1)
        assert r["edges"] == 0 and r["lcc"] == 1
        assert r["edges_p"] > 0.5

    def test_undirected_and_parallel_edges_count_once(self):
        g = nx.MultiDiGraph()
        g.add_edges_from([(0, 1), (1, 0), (0, 1), (1, 2), (2, 2)])
        g.add_nodes_from(range(3, 20))
        r = induced_connectivity(g, [0, 1, 2], n_random=10, seed=0)
        assert r["edges"] == 2 and r["lcc"] == 3

    def test_universe_restricts_nodes_and_draws(self):
        g = chain_with_clique()
        universe = set(range(30))
        r = induced_connectivity(g, [0, 1, 2, 50], universe=universe, n_random=20, seed=0)
        assert r["n_nodes"] == 3 and r["n_universe"] == 30

    def test_reproducible_with_seed(self):
        g = nx.gnm_random_graph(80, 160, seed=3)
        a = induced_connectivity(g, range(10), n_random=50, seed=7)
        b = induced_connectivity(g, range(10), n_random=50, seed=7)
        assert a == b

    def test_random_sets_match_degree_bins(self, monkeypatch):
        g = chain_with_clique()
        seen = []
        import skm_tools.connectivity as c
        orig = c._Pairs.stats

        def spy(self, idx):
            seen.append(np.sort(self.degree[idx]))
            return orig(self, idx)
        monkeypatch.setattr(c._Pairs, "stats", spy)
        induced_connectivity(g, range(6), n_random=5, seed=2, n_bins=2)
        # the clique nodes have the highest degrees: every random set draws from the same (top) bin
        assert all(len(s) == 6 for s in seen)
        assert all(s.min() >= 2 for s in seen[1:])

    def test_too_few_nodes(self):
        g = chain_with_clique()
        assert induced_connectivity(g, [0, "missing"]) == {"n_nodes": 1, "n_universe": g.number_of_nodes()}
