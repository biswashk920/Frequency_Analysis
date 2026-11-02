"""Fit all distributions, compute quantiles, plotting positions and goodness-of-fit."""
import warnings
import numpy as np
from .distributions import fit_gumbel, fit_gev, fit_lognormal, fit_lp3

FITTERS = {"Gumbel": fit_gumbel, "GEV": fit_gev, "Lognormal": fit_lognormal, "LP3": fit_lp3}
LABEL = {"Gumbel": "Gumbel", "GEV": "GEV", "Lognormal": "Log Normal", "LP3": "Log Pearson III"}
RETURN_PERIODS = (2, 5, 10, 25, 50, 100, 200, 500)
_PP_A = {"weibull": 0.0, "cunnane": 0.4, "gringorten": 0.44}


def plotting_positions(n, method="weibull"):
    """Non-exceedance probabilities of the sorted sample (ascending)."""
    a = _PP_A[method]
    return (np.arange(1, n + 1) - a) / (n + 1 - 2 * a)


def quantiles(fit, T):
    """Value for return period(s) T (years): x_T = F^-1(1 - 1/T)."""
    return np.asarray(fit.ppf(1 - 1 / np.asarray(T, float)))


def fit_all(x):
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, fitter in FITTERS.items():
            out[name] = fitter(x)
    return out


def gof(x, fit, pp):
    """KS statistic, RMSE (sorted data vs fitted quantiles at the plotting positions), AIC and BIC.

    AIC/BIC are NaN (and support_ok False) when any observation lies outside the range of the
    fitted distribution, because the likelihood is then zero. KS p-values are not reported:
    with parameters estimated from the same data they are optimistic.
    """
    xs = np.sort(np.asarray(x, float))
    n = len(xs)
    F, i = fit.cdf(xs), np.arange(1, n + 1)
    ks = float(max(np.max(i / n - F), np.max(F - (i - 1) / n)))
    lp = fit.logpdf(xs)
    ok = bool(np.all(np.isfinite(lp)))
    ll = float(np.sum(lp)) if ok else np.nan
    return {
        "KS": ks,
        "RMSE": float(np.sqrt(np.mean((xs - fit.ppf(pp)) ** 2))),
        "AIC": 2 * fit.n_par - 2 * ll if ok else np.nan,
        "BIC": fit.n_par * np.log(n) - 2 * ll if ok else np.nan,
        "support_ok": ok,
    }


def akaike_weights(aic):
    """Return (delta_AIC, Akaike weights). Models with NaN AIC get NaN delta and weight 0."""
    a = np.asarray(aic, float)
    ok = np.isfinite(a)
    delta, w = np.full(a.shape, np.nan), np.zeros(a.shape)
    if ok.any():
        delta[ok] = a[ok] - a[ok].min()
        raw = np.exp(-0.5 * delta[ok])
        w[ok] = raw / raw.sum()
    return delta, w
