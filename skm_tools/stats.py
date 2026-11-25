'''Small statistics helpers (no scipy): hypergeometric tail and Benjamini-Hochberg FDR.'''

import math

import numpy as np


def hypergeom_sf(k, population, successes, draws):
    '''P(X >= k) for X ~ Hypergeometric(population, successes, draws): the one-sided over-representation p-value
    of k hits (e.g. k of a TF's `successes` targets among `draws` DEGs, out of `population` genes).

    Parameters
    ----------
    k, population, successes, draws : int

    Returns
    -------
    float

    Examples
    --------
    >>> round(hypergeom_sf(3, 20, 5, 4), 5)
    0.03199
    '''
    lo, hi = max(k, 0, draws - (population - successes)), min(successes, draws)
    if lo > hi:
        return 0.0 if k > hi else 1.0
    log_total = _log_comb(population, draws)
    terms = [_log_comb(successes, i) + _log_comb(population - successes, draws - i) - log_total
             for i in range(lo, hi + 1)]
    m = max(terms)
    return float(min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms)))


def _log_comb(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def bh_fdr(p):
    '''Benjamini-Hochberg adjusted p-values (FDR), in the order of `p`.

    Examples
    --------
    >>> bh_fdr([0.01, 0.04, 0.03, 0.2]).round(3).tolist()
    [0.04, 0.053, 0.053, 0.2]
    '''
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out
