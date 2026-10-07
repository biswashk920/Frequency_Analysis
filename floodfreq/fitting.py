"""Fit all distributions, compute quantiles, plotting positions and goodness-of-fit."""
import warnings
import numpy as np
from scipy import stats
from .distributions import fit_gumbel, fit_gev, fit_lognormal, fit_lp3

FITTERS = {"Gumbel": fit_gumbel, "GEV": fit_gev, "Lognormal": fit_lognormal, "LP3": fit_lp3}
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
    """KS, RMSE (sorted data vs fitted quantiles), AIC, BIC.
    KS p-values are optimistic because the parameters were estimated from the same data."""
    x = np.asarray(x, float)
    n = len(x)
    ks = stats.kstest(x, fit.cdf)
    ll = float(np.sum(fit.logpdf(x)))
    return {
        "KS": ks.statistic, "KS_p": ks.pvalue,
        "RMSE": float(np.sqrt(np.mean((np.sort(x) - fit.ppf(pp)) ** 2))),
        "AIC": 2 * fit.n_par - 2 * ll, "BIC": fit.n_par * np.log(n) - 2 * ll,
    }
