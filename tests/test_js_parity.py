"""The browser implementation (docs/floodfreq.js) must agree with the Python package."""
import base64, io, json, os, shutil, subprocess, zlib, struct
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from floodfreq import FITTERS, quantiles, gof, plotting_positions, mann_kendall, fit_gev

HERE = os.path.dirname(__file__)
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="Node.js not installed")
RP = [2, 5, 10, 25, 50, 100, 200, 500]


def _png():
    raw = b"\x00\xff\x00\x00" * 4
    def chunk(t, d): return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xFFFFFFFF)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 4, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw[:4 * 4])) + chunk(b"IEND", b"")


def _datasets():
    ex = pd.read_csv(os.path.join(HERE, "..", "examples", "synthetic_peaks.csv")).iloc[:, -1].to_numpy()
    d = {"example": ex,
         "lp3_g1.5_n30": 10 ** stats.pearson3(1.5, loc=3, scale=0.3).rvs(30, random_state=11),
         "lp3_g-1.0_n40": 10 ** stats.pearson3(-1.0, loc=3, scale=0.3).rvs(40, random_state=12),
         "lp3_g2.5_n25": 10 ** stats.pearson3(2.5, loc=2.5, scale=0.4).rvs(25, random_state=13),
         "gev_k-0.3_n35": stats.genextreme(c=-0.3, loc=500, scale=150).rvs(35, random_state=14),
         "ties": np.round(stats.gumbel_r(loc=100, scale=30).rvs(30, random_state=15))}
    for s in range(400):  # a sample whose GEV fit excludes some observations
        x = stats.genextreme(c=0.2, loc=500, scale=150).rvs(30, random_state=s)
        if not gof(x, fit_gev(x), plotting_positions(30))["support_ok"]:
            d["support_violation"] = x; break
    return d


@pytest.fixture(scope="module")
def js(tmp_path_factory):
    data = _datasets(); xl = str(tmp_path_factory.mktemp("x") / "js.xlsx")
    payload = {"datasets": {k: list(map(float, v)) for k, v in data.items()}, "xlsx_dataset": "example", "xlsx_path": xl,
               "png_b64": base64.b64encode(_png()).decode(),
               "parse": ["1,200 1,350 980", "520,610,480", "520, 610, 480; 700", "12,500.5", "1 2 abc"],
               "invalid": {"short": [1, 2, 3], "nan": [1, 2, 3, 4, None], "zero": [1, 2, 0, 4, 5], "const": [5] * 6}}
    out = subprocess.run(["node", os.path.join(HERE, "js_harness.js")], input=json.dumps(payload), capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr
    return {"data": data, "out": json.loads(out.stdout), "xlsx": xl}


def test_quantiles_match(js):
    for name, x in js["data"].items():
        for d, fitter in FITTERS.items():
            py = quantiles(fitter(x), RP); jq = np.array(js["out"]["datasets"][name]["q"][d])
            assert np.allclose(jq, py, rtol=2e-6), (name, d, np.max(np.abs(jq / py - 1)))


def test_statistics_match(js):
    for name, x in js["data"].items():
        pp = plotting_positions(len(x))
        for row in js["out"]["datasets"][name]["stats"]:
            g = gof(x, FITTERS[row["distribution"]](x), pp)
            assert row["support_ok"] == g["support_ok"], (name, row["distribution"])
            for k in ("KS", "RMSE", "AIC", "BIC"):
                if np.isnan(g[k]): assert row[k] is None
                else: assert abs(row[k] - g[k]) <= 1e-6 * max(1, abs(g[k])), (name, row["distribution"], k, row[k], g[k])


def test_support_violation_present_and_agreed(js):
    assert "support_violation" in js["data"]
    rows = js["out"]["datasets"]["support_violation"]["stats"]
    assert any(not r["support_ok"] for r in rows) and abs(sum(r["akaike_weight"] for r in rows) - 1) < 1e-9
    assert any("not available" in s for s in js["out"]["datasets"]["support_violation"]["summary"])


def test_mann_kendall_matches(js):
    for name, x in js["data"].items():
        py, j = mann_kendall(x), js["out"]["datasets"][name]["mk"]
        assert j["S"] == py["S"] and j["trend"] == py["trend"]
        for k in ("Z", "tau", "sen_slope", "var_S"): assert abs(j[k] - py[k]) < 1e-9, (name, k)
        assert abs(j["p"] - py["p"]) < 1e-9


def test_parsing_and_validation_match(js):
    p = js["out"]["parse"]
    assert p["1,200 1,350 980"] == [1200, 1350, 980] and p["520,610,480"] == [520, 610, 480]
    assert p["520, 610, 480; 700"] == [520, 610, 480, 700] and p["12,500.5"] == [12500.5] and p["1 2 abc"] == "error"
    e = js["out"]["errors"]
    assert "at least 5" in e["short"] and "finite" in e["nan"] and "positive" in e["zero"] and "identical" in e["const"]


def test_js_workbook_opens_and_matches_python(js):
    openpyxl = pytest.importorskip("openpyxl")
    wb = openpyxl.load_workbook(js["xlsx"])
    assert wb.sheetnames == ["Summary", "Return Periods", "Fit Statistics", "Parameters", "Data", "Curves", "Figures"]
    assert len(wb["Figures"]._images) == 2
    x = js["data"]["example"]; rp = wb["Return Periods"]
    for k, T in enumerate(RP):
        for i, d in enumerate(FITTERS):
            assert abs(rp.cell(5 + k, 2 + 3 * i).value - quantiles(FITTERS[d](x), [T])[0]) < 1e-3
    assert wb["Summary"]["A11"].value == js["out"]["datasets"]["example"]["summary"][0]
    assert wb["Data"].max_row == len(x) + 1
