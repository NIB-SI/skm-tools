import numpy as np
import pytest

from skm_tools.stats import bh_fdr, hypergeom_sf


@pytest.mark.parametrize("k,M,n,N", [(3, 20, 5, 4), (0, 50, 10, 5), (5, 100, 20, 10), (10, 18000, 200, 300),
                                     (6, 10, 5, 5), (1, 10, 0, 3)])
def test_hypergeom_sf_matches_scipy(k, M, n, N):
    stats = pytest.importorskip("scipy.stats")
    assert hypergeom_sf(k, M, n, N) == pytest.approx(stats.hypergeom.sf(k - 1, M, n, N), rel=1e-9, abs=1e-300)


def test_hypergeom_sf_bounds():
    assert hypergeom_sf(0, 30, 10, 5) == pytest.approx(1.0)
    assert hypergeom_sf(6, 30, 10, 5) == 0.0          # more hits than draws


def test_bh_fdr_matches_statsmodels():
    multitest = pytest.importorskip("statsmodels.stats.multitest")
    p = np.random.default_rng(1).uniform(size=50) ** 3
    assert np.allclose(bh_fdr(p), multitest.multipletests(p, method="fdr_bh")[1])


def test_bh_fdr_empty():
    assert len(bh_fdr([])) == 0
