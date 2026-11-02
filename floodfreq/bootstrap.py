"""Parametric bootstrap confidence intervals for quantiles."""
import warnings
import numpy as np
from .fitting import quantiles


def bootstrap_ci(fit, fitter, n, T, nboot=2000, conf=0.95, seed=42):
    """Resample n values from the fitted distribution, refit, collect quantiles.
    Returns (lower, upper) percentile bounds at each return period in T."""
    rng = np.random.default_rng(seed)
    T = np.asarray(T, float)
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for _ in range(nboot):
            try:
                q = quantiles(fitter(fit.rvs(n, rng)), T)
            except Exception:
                continue
            if np.all(np.isfinite(q)):
                out.append(q)
    if len(out) < 0.5 * nboot:
        raise RuntimeError("Too many failed bootstrap fits; data may be unsuitable for this distribution.")
    a = (1 - conf) / 2
    out = np.array(out)
    return np.quantile(out, a, axis=0), np.quantile(out, 1 - a, axis=0)
