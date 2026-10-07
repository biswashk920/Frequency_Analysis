import numpy as np
from scipy import stats
from floodfreq import FITTERS, fit_gumbel, fit_gev, fit_lp3, lmoments, quantiles, bootstrap_ci

RNG = np.random.default_rng(1)


def test_gumbel_recovery():
    f = fit_gumbel(stats.gumbel_r(loc=100, scale=30).rvs(20000, random_state=1))
    assert abs(f.params["location"] - 100) < 1.5 and abs(f.params["scale"] - 30) < 1


def test_gev_recovery():
    f = fit_gev(stats.genextreme(c=-0.15, loc=500, scale=120).rvs(20000, random_state=2))
    assert abs(f.params["shape_k"] + 0.15) < 0.03
    assert abs(f.params["location"] - 500) < 5 and abs(f.params["scale"] - 120) < 4


def test_lp3_recovery():
    y = stats.pearson3(0.4, loc=3.0, scale=0.25).rvs(20000, random_state=3)
    f = fit_lp3(10 ** y)
    assert abs(f.params["skew_log10"] - 0.4) < 0.08 and abs(f.params["mean_log10"] - 3.0) < 0.01


def test_roundtrip_and_monotonic():
    x = stats.genextreme(c=-0.1, loc=800, scale=250).rvs(50, random_state=4)
    F = np.array([0.1, 0.5, 0.9, 0.99])
    for name, fitter in FITTERS.items():
        f = fitter(x)
        assert np.allclose(f.cdf(f.ppf(F)), F, atol=1e-6), name
        assert np.all(np.diff(quantiles(f, [2, 5, 10, 100, 500])) > 0), name


def test_lmoments_known():
    l1, l2, t3 = lmoments(np.arange(1, 101))  # uniform-like data: t3 ~ 0
    assert abs(l1 - 50.5) < 1e-9 and abs(t3) < 1e-9


def test_ci_brackets_estimate():
    x = stats.gumbel_r(loc=100, scale=30).rvs(40, random_state=5)
    f = fit_gumbel(x)
    lo, hi = bootstrap_ci(f, fit_gumbel, 40, [10, 100], nboot=300)
    q = quantiles(f, [10, 100])
    assert np.all(lo < q) and np.all(q < hi)


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
