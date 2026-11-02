"""Optional trend screening: Mann-Kendall test with Sen's slope."""
import numpy as np
from scipy import stats


def mann_kendall(x, alpha=0.05):
    """Mann-Kendall trend test (tie-corrected, normal approximation with continuity correction).

    `x` must be in chronological order. Serial correlation is NOT tested or corrected.
    Returns a dict with S, var_S, Z, p (two-sided), tau, sen_slope (per time step) and
    trend in {"increasing", "decreasing", "none"} at significance level `alpha`.
    """
    x = np.asarray(x, float)
    n = len(x)
    if n < 4:
        raise ValueError("The Mann-Kendall test needs at least 4 values.")
    i, j = np.triu_indices(n, 1)
    d = x[j] - x[i]
    S = int(np.sum(np.sign(d)))
    _, counts = np.unique(x, return_counts=True)
    t = counts[counts > 1].astype(float)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(t * (t - 1) * (2 * t + 5))) / 18.0
    if var <= 0:
        raise ValueError("All values are tied; the Mann-Kendall test is undefined.")
    z = 0.0 if S == 0 else (S - np.sign(S)) / np.sqrt(var)
    p = float(2 * stats.norm.sf(abs(z)))
    trend = "none" if p >= alpha else ("increasing" if z > 0 else "decreasing")
    return {"n": n, "S": S, "var_S": float(var), "Z": float(z), "p": p, "alpha": alpha,
            "tau": S / (n * (n - 1) / 2), "sen_slope": float(np.median(d / (j - i))), "trend": trend}
