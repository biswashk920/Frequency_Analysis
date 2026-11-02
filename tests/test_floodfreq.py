import os, tempfile
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from floodfreq import (FITTERS, fit_all, fit_gumbel, fit_gev, fit_lp3, lmoments, quantiles, bootstrap_ci, gof,
                       plotting_positions, akaike_weights, mann_kendall, parse_numbers, read_csv_values,
                       best_model_summary, analyze, validate_values)

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "examples", "synthetic_peaks.csv")


# ---------- fitting (unchanged behaviour) ----------
def test_gumbel_recovery():
    f = fit_gumbel(stats.gumbel_r(loc=100, scale=30).rvs(20000, random_state=1))
    assert abs(f.params["location"] - 100) < 1.5 and abs(f.params["scale"] - 30) < 1


def test_gev_recovery():
    f = fit_gev(stats.genextreme(c=-0.15, loc=500, scale=120).rvs(20000, random_state=2))
    assert abs(f.params["shape_k"] + 0.15) < 0.03
    assert abs(f.params["location"] - 500) < 5 and abs(f.params["scale"] - 120) < 4


def test_lp3_recovery():
    y = stats.pearson3(0.4, loc=3.0, scale=0.25).rvs(20000, random_state=3)
    f = fit_lp3(10 ** y)
    assert abs(f.params["skew_log10"] - 0.4) < 0.08 and abs(f.params["mean_log10"] - 3.0) < 0.01


def test_roundtrip_and_monotonic():
    x = stats.genextreme(c=-0.1, loc=800, scale=250).rvs(50, random_state=4)
    F = np.array([0.1, 0.5, 0.9, 0.99])
    for name, fitter in FITTERS.items():
        f = fitter(x)
        assert np.allclose(f.cdf(f.ppf(F)), F, atol=1e-6), name
        assert np.all(np.diff(quantiles(f, [2, 5, 10, 100, 500])) > 0), name


def test_lmoments_known():
    l1, l2, t3 = lmoments(np.arange(1, 101))
    assert abs(l1 - 50.5) < 1e-9 and abs(t3) < 1e-9


def test_plotting_positions():
    assert np.allclose(plotting_positions(4, "weibull"), [0.2, 0.4, 0.6, 0.8])
    assert np.allclose(plotting_positions(4, "cunnane"), (np.arange(1, 5) - 0.4) / 4.2)
    assert np.allclose(plotting_positions(4, "gringorten"), (np.arange(1, 5) - 0.44) / 4.12)


def test_ci_brackets_estimate():
    x = stats.gumbel_r(loc=100, scale=30).rvs(40, random_state=5)
    f = fit_gumbel(x)
    lo, hi = bootstrap_ci(f, fit_gumbel, 40, [10, 100], nboot=300)
    q = quantiles(f, [10, 100])
    assert np.all(lo < q) and np.all(q < hi)


def test_bootstrap_coverage_sanity():
    """With the correct model the 95% band should contain the true 100-year flood most of the time.
    (Measured coverage is a little below nominal for small n; this guards against gross regressions.)"""
    true = stats.genextreme(c=-0.1, loc=800, scale=250); q100 = true.ppf(0.99); hit = 0; R = 40
    for s in range(R):
        x = true.rvs(40, random_state=500 + s); f = fit_gev(x)
        lo, hi = bootstrap_ci(f, fit_gev, 40, [100], nboot=100, seed=s)
        hit += lo[0] <= q100 <= hi[0]
    assert hit / R >= 0.75


# ---------- goodness of fit, AIC, support ----------
def test_gof_matches_scipy_ks():
    x = stats.gumbel_r(loc=100, scale=30).rvs(40, random_state=7)
    f = fit_gumbel(x)
    assert abs(gof(x, f, plotting_positions(40))["KS"] - stats.kstest(x, f.cdf).statistic) < 1e-12


def _violating_sample():
    for s in range(400):
        x = stats.genextreme(c=0.2, loc=500, scale=150).rvs(30, random_state=s)
        if not gof(x, fit_gev(x), plotting_positions(30))["support_ok"]:
            return x
    raise AssertionError("no support violation found")


def test_support_violation_is_flagged_not_inf():
    x = _violating_sample()
    g = gof(x, fit_gev(x), plotting_positions(30))
    assert g["support_ok"] is False and np.isnan(g["AIC"]) and np.isnan(g["BIC"])
    res = analyze(x, nboot=60)
    st = res["fit_statistics"].set_index("distribution")
    assert st.loc["GEV", "akaike_weight"] == 0 and np.isclose(st["akaike_weight"].sum(), 1)
    assert any("not available" in s for s in res["summary"])


def test_akaike_weights():
    d, w = akaike_weights([100.0, 101.0, np.nan, 110.0])
    assert np.isclose(w.sum(), 1) and w[2] == 0 and np.isnan(d[2]) and d[0] == 0 and w[0] > w[1] > w[3]


def test_summary_wording_cases():
    row = lambda d, aic, ks=0.1, rm=1.0: {"distribution": d, "AIC": aic, "KS": ks, "RMSE": rm, "support_ok": True}
    def rows(aics):
        out = [row(d, a) for d, a in zip(["Gumbel", "GEV", "Lognormal", "LP3"], aics)]
        dl, w = akaike_weights([r["AIC"] for r in out])
        for r, a, b in zip(out, dl, w): r["delta_AIC"], r["akaike_weight"] = a, b
        return out
    q = {"Gumbel": 100, "GEV": 110, "Lognormal": 105, "LP3": 104}
    tie = best_model_summary(rows([10, 11, 30, 30]), q, 60)
    assert "cannot separate" in tie[1] and "GEV" in tie[1]
    clear = best_model_summary(rows([10, 20, 30, 30]), q, 60)
    assert "No other model is within" in clear[1]
    assert any("below 30" in s for s in best_model_summary(rows([10, 20, 30, 30]), q, 20))


# ---------- inputs ----------
def test_parse_numbers():
    v, notes = parse_numbers("1,200 1,350 980")
    assert list(v) == [1200, 1350, 980] and len(notes) == 2
    assert list(parse_numbers("520,610,480")[0]) == [520, 610, 480]
    assert list(parse_numbers("520, 610, 480; 700\n800")[0]) == [520, 610, 480, 700, 800]
    assert list(parse_numbers("12,500.5")[0]) == [12500.5]
    with pytest.raises(ValueError):
        parse_numbers("1 2 abc")


def test_csv_year_column_never_analysed(tmp_path):
    a = tmp_path / "a.csv"; pd.DataFrame({"peak": [520, 610, 480, 905, 730, 810], "year": [2001, 2002, 2003, 2004, 2005, 2006]}).to_csv(a, index=False)
    v, years, notes, col = read_csv_values(str(a))
    assert col == "peak" and list(v) == [520, 610, 480, 905, 730, 810] and years[0] == 2001
    b = tmp_path / "b.csv"; pd.DataFrame({"Year": [2004, 2002, 2003, 2001, 2005, 2006], "Q": [4, 2, 3, 1, 5, 6]}).to_csv(b, index=False)
    v, years, _, col = read_csv_values(str(b))
    assert col == "Q" and list(v) == [1, 2, 3, 4, 5, 6] and years == [2001, 2002, 2003, 2004, 2005, 2006]
    c = tmp_path / "c.csv"; pd.DataFrame({"flow": [5, 4, 9, 7, 3, 8], "wy": [1990, 1991, 1992, 1993, 1994, 1995]}).to_csv(c, index=False)
    assert read_csv_values(str(c))[3] == "flow"


@pytest.mark.parametrize("bad,msg", [([1, 2, 3], "at least 5"), ([1, 2, 3, np.nan, 5], "finite"), ([1, 2, 3, np.inf, 5], "finite"),
                                      ([1, 2, 0, 4, 5], "positive"), ([7] * 8, "identical")])
def test_validation(bad, msg):
    with pytest.raises(ValueError, match=msg):
        validate_values(bad)


# ---------- Mann-Kendall ----------
def test_mann_kendall_known_values():
    r = mann_kendall(np.arange(1.0, 11.0))
    assert r["S"] == 45 and r["trend"] == "increasing" and abs(r["tau"] - 1) < 1e-12 and abs(r["sen_slope"] - 1) < 1e-12
    assert abs(r["Z"] - 44 / np.sqrt(125)) < 1e-12 and abs(r["p"] - 2 * stats.norm.sf(44 / np.sqrt(125))) < 1e-12
    d = mann_kendall(np.arange(10.0, 0.0, -1))
    assert d["S"] == -45 and d["trend"] == "decreasing" and abs(d["sen_slope"] + 1) < 1e-12


def test_mann_kendall_ties_and_tau():
    x = np.array([1, 1, 2, 2, 3.0])
    r = mann_kendall(x)
    var = (5 * 4 * 15 - (2 * 1 * 9 + 2 * 1 * 9)) / 18
    assert abs(r["var_S"] - var) < 1e-12
    y = stats.norm.rvs(size=30, random_state=3)
    assert abs(mann_kendall(y)["tau"] - stats.kendalltau(np.arange(30), y).statistic) < 1e-12
    assert mann_kendall(stats.norm.rvs(size=40, random_state=11) + 10)["trend"] == "none"


# ---------- whole pipeline, Excel, CLI ----------
def test_analyze_pipeline():
    x = pd.read_csv(EXAMPLE).iloc[:, -1].to_numpy()
    r = analyze(x, nboot=80, trend=True)
    st = r["fit_statistics"]
    assert np.isclose(st["akaike_weight"].sum(), 1) and r["mk"] is not None and r["summary"][0].startswith("Best fit by AIC")
    rp = r["return_periods"]
    assert (rp["lower"] <= rp["estimate"]).all() and (rp["estimate"] <= rp["upper"]).all()


def test_excel_workbook(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    from floodfreq.excel import write_workbook
    x = pd.read_csv(EXAMPLE).iloc[:, -1].to_numpy()
    r = analyze(x, nboot=80, trend=True, years=list(range(1961, 1961 + len(x))))
    p = tmp_path / "out.xlsx"; write_workbook(str(p), r, "Discharge (m\u00b3/s)", "Test")
    wb = openpyxl.load_workbook(p)
    assert wb.sheetnames == ["Summary", "Return Periods", "Fit Statistics", "Parameters", "Data", "Curves", "Charts"]
    assert len(wb["Charts"]._charts) == 5
    text = " ".join(str(c.value) for row in wb["Summary"].iter_rows() for c in row if c.value)
    assert r["summary"][0] in text and "Mann-Kendall" in text
    rp = wb["Return Periods"]
    for k, T in enumerate([2, 5, 10, 25, 50, 100, 200, 500]):
        assert rp.cell(5 + k, 1).value == T
        assert abs(rp.cell(5 + k, 2).value - quantiles(r["fits"]["Gumbel"], [T])[0]) < 1e-9
    assert wb["Data"].cell(2, 1).value == 1961


def test_cli_end_to_end(tmp_path, capsys):
    from floodfreq.cli import main
    main([EXAMPLE, "--nboot", "60", "--trend", "--outdir", str(tmp_path), "--station", "T"])
    for f in ["return_periods.csv", "fit_statistics.csv", "fit_all.png", "compare.png", "floodfreq_results.xlsx", "fit_GEV.png"]:
        assert (tmp_path / f).exists(), f
    out = capsys.readouterr().out
    assert "Which model fits best?" in out and "Using column" in out
    st = pd.read_csv(tmp_path / "fit_statistics.csv")
    assert {"delta_AIC", "akaike_weight", "support_ok"} <= set(st.columns) and "KS_p" not in st.columns


def test_cli_rejects_bad_input():
    from floodfreq.cli import main
    for argv in (["--values", "1 2 3 nan 5 6"], ["--values", "4 4 4 4 4 4"], ["nonexistent.csv"]):
        with pytest.raises(SystemExit):
            main(argv + ["--outdir", tempfile.mkdtemp()])


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
