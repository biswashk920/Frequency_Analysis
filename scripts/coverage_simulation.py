"""How often do the nominal 95% bootstrap bands contain the true 100- and 500-year flood?

    python scripts/coverage_simulation.py [replications] [n]

Truth is GEV(k = -0.1). Takes about a minute for the defaults. The standard error of each percentage is
roughly 3-5 points at 100 replications, so treat the output as indicative.
"""
import sys
import numpy as np
from scipy import stats
from floodfreq import FITTERS, bootstrap_ci

R = int(sys.argv[1]) if len(sys.argv) > 1 else 100
n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
true = stats.genextreme(c=-0.1, loc=800, scale=250)
Q100, Q500 = true.ppf(0.99), true.ppf(1 - 1 / 500)
hit = {d: [0, 0] for d in FITTERS}
for s in range(R):
    x = true.rvs(n, random_state=1000 + s)
    for d, fit in FITTERS.items():
        lo, hi = bootstrap_ci(fit(x), fit, n, [100, 500], nboot=150, conf=0.95, seed=s)
        hit[d][0] += lo[0] <= Q100 <= hi[0]
        hit[d][1] += lo[1] <= Q500 <= hi[1]
print(f"truth GEV(k=-0.1), n={n}, {R} replications, nominal 95% bands")
for d in FITTERS:
    print(f"  {d:10s} coverage of 100-yr: {hit[d][0] / R:4.0%}   of 500-yr: {hit[d][1] / R:4.0%}")
