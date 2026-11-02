"""Plain-language summary of which model fits best."""
import numpy as np
from .fitting import LABEL


def fmt(v):
    v = float(v)
    return f"{v:.0f}" if v >= 1000 else f"{v:.1f}" if v >= 100 else f"{v:.2f}"


def _pfmt(p):
    return "< 0.001" if p < 0.001 else f"{p:.3f}"


def _names(rs):
    n = [LABEL[r["distribution"]] for r in rs]
    return n[0] if len(n) == 1 else ", ".join(n[:-1]) + " and " + n[-1]


def best_model_summary(rows, q100, n, gev_k=None, mk=None):
    """rows: dicts with distribution, KS, RMSE, AIC, delta_AIC, akaike_weight, support_ok.
    q100: {distribution: 100-year estimate}. Returns a list of sentences."""
    lines = []
    ok = sorted([r for r in rows if np.isfinite(r["AIC"])], key=lambda r: r["AIC"])
    if ok:
        b = ok[0]
        lines.append(f"Best fit by AIC: {LABEL[b['distribution']]} (Akaike weight {b['akaike_weight']:.2f}).")
        close = [r for r in ok[1:] if r["delta_AIC"] < 2]
        if close:
            lines.append(f"{_names(close)} {'is' if len(close) == 1 else 'are'} within \u0394AIC < 2 of it, "
                         f"so the data cannot separate {'the two' if len(close) == 1 else 'them'}; "
                         "compare the upper tails and use hydrological judgement.")
        else:
            lines.append(f"No other model is within \u0394AIC = 2, so AIC favours {LABEL[b['distribution']]} over the rest.")
        ks = min(rows, key=lambda r: r["KS"])["distribution"]
        rm = min(rows, key=lambda r: r["RMSE"])["distribution"]
        if ks == rm == b["distribution"]:
            lines.append(f"KS, RMSE and AIC all point to {LABEL[ks]}.")
        else:
            lines.append(f"Other criteria: lowest KS = {LABEL[ks]}, lowest RMSE = {LABEL[rm]}. These can disagree "
                         "with AIC, and RMSE depends on the plotting position.")
    else:
        lines.append("AIC could not be computed for any model (observations fall outside every fitted range).")
    bad = [r for r in rows if not r["support_ok"]]
    if bad:
        lines.append(f"AIC/BIC not available for {_names(bad)}: some observations fall outside the range of the fitted distribution.")
    qs = {d: v for d, v in q100.items() if np.isfinite(v)}
    if len(qs) > 1:
        lo, hi = min(qs, key=qs.get), max(qs, key=qs.get)
        lines.append(f"The 100-year estimate ranges from {fmt(qs[lo])} ({LABEL[lo]}) to {fmt(qs[hi])} ({LABEL[hi]}) "
                     f"across the models, a spread of {(qs[hi] / qs[lo] - 1) * 100:.0f}%. This model-to-model difference "
                     "is separate from the bootstrap bands, which only reflect sampling noise within each fitted model.")
    if n < 30:
        lines.append(f"n = {n} is below 30, so the 100\u2013500 year estimates are heavy extrapolation.")
    if gev_k is not None and abs(gev_k) > 0.4:
        lines.append(f"GEV shape k = {gev_k:.2f} is near the limit of the L-moment approximation.")
    if mk is not None:
        if mk["trend"] != "none":
            lines.append(f"Mann-Kendall found a significant {mk['trend']} trend (p = {_pfmt(mk['p'])}); "
                         "the stationarity assumption behind all four fits is questionable.")
        else:
            lines.append(f"Mann-Kendall found no significant trend (p = {_pfmt(mk['p'])}); serial correlation was not tested.")
    return lines
