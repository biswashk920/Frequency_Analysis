"""One call that runs the whole analysis and returns everything the CLI, Excel export and tests need."""
import numpy as np
import pandas as pd
from .fitting import FITTERS, RETURN_PERIODS, fit_all, quantiles, plotting_positions, gof, akaike_weights
from .bootstrap import bootstrap_ci
from .screening import mann_kendall
from .report import best_model_summary


def validate_values(x):
    x = np.asarray(x, float)
    if x.size < 5:
        raise ValueError("Need at least 5 values.")
    if not np.all(np.isfinite(x)):
        raise ValueError("Values must be finite numbers (found NaN or infinity).")
    if np.any(x <= 0):
        raise ValueError("All values must be positive (Lognormal and LP3 use logarithms).")
    if np.ptp(x) == 0:
        raise ValueError("All values are identical; a distribution cannot be fitted.")
    return x


def analyze(x, conf=0.95, nboot=2000, pp_method="weibull", trend=False, seed=42, years=None):
    x = validate_values(x)
    if not 0 < conf < 1:
        raise ValueError("Confidence level must be between 0 and 1.")
    if nboot < 50:
        raise ValueError("Use at least 50 bootstrap simulations.")
    n = len(x)
    RP = np.array(RETURN_PERIODS, float)
    grid = np.unique(np.concatenate([np.geomspace(1.02, 500, 45), RP]))
    sel = np.isin(grid, RP)
    pp = plotting_positions(n, pp_method)
    fits = fit_all(x)
    bands, rows, srows = {}, [], []
    for name, fit in fits.items():
        lo, hi = bootstrap_ci(fit, FITTERS[name], n, grid, nboot, conf, seed)
        bands[name] = (lo, hi)
        for T, q, l, h in zip(RP, quantiles(fit, RP), lo[sel], hi[sel]):
            rows.append({"distribution": name, "return_period_yr": int(T), "estimate": q, "lower": l, "upper": h})
        srows.append({"distribution": name, **gof(x, fit, pp), **{f"param_{k}": v for k, v in fit.params.items()}})
    st = pd.DataFrame(srows)
    delta, w = akaike_weights(st["AIC"])
    st.insert(5, "delta_AIC", delta)
    st.insert(6, "akaike_weight", w)
    st["rank_AIC"] = st["AIC"].rank(method="min")
    mk = mann_kendall(x) if trend else None
    res = pd.DataFrame(rows)
    q100 = {d: float(res[(res.distribution == d) & (res.return_period_yr == 100)].estimate.iloc[0]) for d in fits}
    summary = best_model_summary(st.to_dict("records"), q100, n, fits["GEV"].params["shape_k"], mk)
    return {"x": x, "years": years, "n": n, "conf": conf, "nboot": nboot, "pp_method": pp_method,
            "grid": grid, "pp": pp, "fits": fits, "bands": bands, "return_periods": res,
            "fit_statistics": st, "q100": q100, "mk": mk, "summary": summary}
