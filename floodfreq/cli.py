"""Command line interface.

  python -m floodfreq.cli data.csv --type discharge --station "Gauge 1"
  python -m floodfreq.cli --values "520 610 480 905 730 ..." --type rainfall
"""
import argparse, os, re, sys
import numpy as np
import pandas as pd
from .fitting import FITTERS, RETURN_PERIODS, fit_all, quantiles, plotting_positions, gof
from .bootstrap import bootstrap_ci
from . import plotting as P

UNITS = {"discharge": "Discharge (m³/s)", "rainfall": "Rainfall (mm)"}


def read_values(src, values, column):
    if values:
        raw = re.split(r"[,\s;]+", values.strip())
    elif src and os.path.isfile(src):
        df = pd.read_csv(src)
        numeric = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if not numeric:
            sys.exit("No numeric column found in the CSV.")
        col = column or numeric[-1]  # last numeric column (a year column usually comes first)
        print(f"Using column '{col}'")
        raw = df[col].dropna().tolist()
    else:
        sys.exit("Provide a CSV file path or --values '1 2 3 ...'")
    try:
        return np.array([float(v) for v in raw if str(v) != ""])
    except ValueError as e:
        sys.exit(f"Could not parse numbers: {e}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Flood frequency analysis (Gumbel, GEV, Lognormal, LP3)")
    ap.add_argument("source", nargs="?", help="CSV file with annual maxima")
    ap.add_argument("--values", help="Numbers pasted as a string")
    ap.add_argument("--column", help="CSV column name (default: last numeric)")
    ap.add_argument("--type", choices=UNITS, default="discharge")
    ap.add_argument("--station", default="")
    ap.add_argument("--conf", type=float, default=0.95)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--pp", choices=["weibull", "cunnane", "gringorten"], default="weibull")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args(argv)

    x = read_values(a.source, a.values, a.column)
    n = len(x)
    if n < 5:
        sys.exit("Need at least 5 values.")
    if np.any(x <= 0):
        sys.exit("All values must be positive (Lognormal and LP3 use logarithms).")
    if n < 30:
        print(f"WARNING: n = {n} < 30. Long return periods (100-500 yr) are heavy extrapolation.\n")

    os.makedirs(a.outdir, exist_ok=True)
    unit, RP = UNITS[a.type], np.array(RETURN_PERIODS, float)
    grid = np.unique(np.concatenate([np.geomspace(1.02, 500, 45), RP]))
    pp = plotting_positions(n, a.pp)
    fits = fit_all(x)
    rows, stats_rows, bands = [], [], {}
    for name, fit in fits.items():
        lo, hi = bootstrap_ci(fit, FITTERS[name], n, grid, a.nboot, a.conf)
        bands[name] = (lo, hi)
        sel = np.isin(grid, RP)
        for T, q, l, h in zip(RP, quantiles(fit, RP), lo[sel], hi[sel]):
            rows.append({"distribution": name, "return_period_yr": int(T), "estimate": q, "lower": l, "upper": h})
        stats_rows.append({"distribution": name, **gof(x, fit, pp), **{f"param_{k}": v for k, v in fit.params.items()}})
        P.plot_one(os.path.join(a.outdir, f"fit_{name}.png"), x, name, fit, grid, lo, hi, pp, unit, a.conf, a.station)
    res, st = pd.DataFrame(rows), pd.DataFrame(stats_rows)
    st["best_AIC"] = st["AIC"] == st["AIC"].min()
    res.to_csv(os.path.join(a.outdir, "return_periods.csv"), index=False)
    st.to_csv(os.path.join(a.outdir, "fit_statistics.csv"), index=False)
    P.plot_all(os.path.join(a.outdir, "fit_all.png"), x, fits, bands, grid, pp, unit, a.conf, a.station)
    P.plot_compare(os.path.join(a.outdir, "compare.png"), x, fits, grid, pp, unit, a.station)

    k = fits["GEV"].params["shape_k"]
    if abs(k) > 0.4:
        print(f"WARNING: GEV shape k = {k:.2f} is near the limit of the L-moment approximation.\n")
    wide = res.pivot(index="return_period_yr", columns="distribution", values="estimate")[list(FITTERS)]
    print(f"{a.station or 'Station'}  n={n}  ({unit})\n")
    print(wide.round(1).to_string()); print()
    print(st[["distribution", "KS", "RMSE", "AIC", "BIC", "best_AIC"]].round(3).to_string(index=False))
    print(f"\nResults written to {a.outdir}/")


if __name__ == "__main__":
    main()
