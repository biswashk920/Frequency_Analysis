"""Distribution fitting. Each fit_* returns a `Fitted` object with ppf/cdf/logpdf/rvs.

Gumbel and GEV: L-moments (Hosking 1990).  Lognormal: moments of ln(x).
LP3: moments of log10(x) with sample skew (Bulletin 17B style, no regional skew).
"""
import numpy as np
from scipy import stats
from scipy.special import gamma as G

EULER = 0.5772156649015329


def lmoments(x):
    """Return (l1, l2, t3): first two sample L-moments and L-skewness."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    i = np.arange(1, n + 1)
    b0 = x.mean()
    b1 = np.sum((i - 1) / (n - 1) * x) / n
    b2 = np.sum((i - 1) * (i - 2) / ((n - 1) * (n - 2)) * x) / n
    l1, l2, l3 = b0, 2 * b1 - b0, 6 * b2 - 6 * b1 + b0
    return l1, l2, l3 / l2


class Fitted:
    """A fitted distribution in the original (flow / rainfall) units."""

    def __init__(self, name, rv, params, n_par, log10=False):
        self.name, self.rv, self.params, self.n_par, self.log10 = name, rv, params, n_par, log10

    def ppf(self, F):
        v = self.rv.ppf(F)
        return 10 ** v if self.log10 else v

    def cdf(self, x):
        x = np.asarray(x, float)
        return self.rv.cdf(np.log10(x) if self.log10 else x)

    def logpdf(self, x):
        x = np.asarray(x, float)
        if not self.log10:
            return self.rv.logpdf(x)
        return self.rv.logpdf(np.log10(x)) - np.log(x * np.log(10))

    def rvs(self, n, rng):
        """Inverse-transform sampling (works for every distribution here)."""
        return self.ppf(rng.uniform(1e-9, 1 - 1e-9, size=n))


def fit_gumbel(x):
    l1, l2, _ = lmoments(x)
    a = l2 / np.log(2)
    xi = l1 - EULER * a
    return Fitted("Gumbel", stats.gumbel_r(loc=xi, scale=a), {"location": xi, "scale": a}, 2)


def fit_gev(x):
    l1, l2, t3 = lmoments(x)
    c = 2 / (3 + t3) - np.log(2) / np.log(3)
    k = float(np.clip(7.8590 * c + 2.9554 * c ** 2, -0.5, 0.5))  # Hosking's approximation range
    if abs(k) < 1e-6:
        g = fit_gumbel(x)
        return Fitted("GEV", stats.genextreme(c=0.0, loc=g.params["location"], scale=g.params["scale"]),
                      {"location": g.params["location"], "scale": g.params["scale"], "shape_k": 0.0}, 3)
    a = l2 * k / ((1 - 2 ** -k) * G(1 + k))
    xi = l1 - a * (1 - G(1 + k)) / k
    # scipy's genextreme `c` equals Hosking's k (k > 0: bounded upper tail)
    return Fitted("GEV", stats.genextreme(c=k, loc=xi, scale=a),
                  {"location": xi, "scale": a, "shape_k": k}, 3)


def fit_lognormal(x):
    ly = np.log(np.asarray(x, float))
    mu, s = ly.mean(), ly.std(ddof=1)
    return Fitted("Lognormal", stats.lognorm(s=s, scale=np.exp(mu)), {"mu_ln": mu, "sigma_ln": s}, 2)


def fit_lp3(x):
    y = np.log10(np.asarray(x, float))
    m, s, g = y.mean(), y.std(ddof=1), float(stats.skew(y, bias=False))
    return Fitted("LP3", stats.pearson3(g, loc=m, scale=s),
                  {"mean_log10": m, "sd_log10": s, "skew_log10": g}, 3, log10=True)
