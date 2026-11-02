"""Command line interface.

  python -m floodfreq.cli data.csv --type discharge --station "Gauge 1"
  python -m floodfreq.cli --values "520 610 480 905 730 ..." --type rainfall --trend
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from .fitting import LABEL, FITTERS, RETURN_PERIODS
from .inputs import parse_numbers, read_csv_values
from .analysis import analyze
from . import plotting as P

UNITS = {"discharge": "Discharge (m\u00b3/s)", "rainfall": "Rainfall (mm)"}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Flood frequency analysis (Gumbel, GEV, Log Normal, Log Pearson III)")
    ap.add_argument("source", nargs="?", help="CSV file with annual maxima (a year column, if present, is detected)")
    ap.add_argument("--values", help="Numbers pasted as a string (in chronological order if you use --trend)")
    ap.add_argument("--column", help="CSV column holding the values (default: last numeric non-year column)")
    ap.add_argument("--type", choices=UNITS, default="discharge")
    ap.add_argument("--station", default="")
    ap.add_argument("--conf", type=float, default=0.95)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--pp", choices=["weibull", "cunnane", "gringorten"], default="weibull")
    ap.add_argument("--trend", action="store_true", help="also run the (optional) Mann-Kendall trend test")
    ap.add_argument("--no-xlsx", action="store_true", help="skip the Excel workbook")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args(argv)

    years = None
    try:
        if a.values:
            x, notes = parse_numbers(a.values)
        elif a.source:
            if not os.path.isfile(a.source):
                sys.exit(f"File not found: {a.source}")
            x, years, notes, col = read_csv_values(a.source, a.column)
            print(f"Using column '{col}'")
        else:
            sys.exit("Provide a CSV file path or --values '1 2 3 ...'")
        for note in notes:
            print("Note:", note)
        res = analyze(x, a.conf, a.nboot, a.pp, a.trend, years=years)
    except ValueError as e:
        sys.exit(f"Error: {e}")

    os.makedirs(a.outdir, exist_ok=True)
    unit, x, n, grid, pp = UNITS[a.type], res["x"], res["n"], res["grid"], res["pp"]
    fits, bands = res["fits"], res["bands"]
    for name, fit in fits.items():
        P.plot_one(os.path.join(a.outdir, f"fit_{name}.png"), x, name, fit, grid, *bands[name], pp, unit, a.conf, a.station)
    P.plot_all(os.path.join(a.outdir, "fit_all.png"), x, fits, bands, grid, pp, unit, a.conf, a.station)
    P.plot_compare(os.path.join(a.outdir, "compare.png"), x, fits, grid, pp, unit, a.station)
    res["return_periods"].to_csv(os.path.join(a.outdir, "return_periods.csv"), index=False)
    res["fit_statistics"].to_csv(os.path.join(a.outdir, "fit_statistics.csv"), index=False)
    if not a.no_xlsx:
        from .excel import write_workbook
        write_workbook(os.path.join(a.outdir, "floodfreq_results.xlsx"), res, unit, a.station, a.type)

    wide = res["return_periods"].pivot(index="return_period_yr", columns="distribution", values="estimate")[list(FITTERS)]
    st = res["fit_statistics"]
    print(f"\n{a.station or 'Station'}  n={n}  ({unit})\n")
    print(wide.round(1).to_string()); print()
    print(st[["distribution", "KS", "RMSE", "AIC", "delta_AIC", "akaike_weight"]].round(3).to_string(index=False))
    print("\nWhich model fits best?")
    for line in res["summary"]:
        print(" -", line)
    print(f"\nResults written to {a.outdir}/")


if __name__ == "__main__":
    main()
