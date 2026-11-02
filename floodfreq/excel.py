"""Excel workbook with tables for every method and native (editable) Excel charts."""
import math
import numpy as np
import xlsxwriter
from ._version import __version__
from .fitting import FITTERS, LABEL, RETURN_PERIODS

COLORS = {"Gumbel": "#1f77b4", "GEV": "#d62728", "Lognormal": "#2ca02c", "LP3": "#9467bd"}
PARAM_INFO = {
    "location": "Location \u03be (same units as the data)", "scale": "Scale \u03b1 (same units as the data)",
    "shape_k": "Shape k (k < 0: heavy upper tail; k > 0: bounded upper tail)",
    "mu_ln": "Mean of ln(x)", "sigma_ln": "Standard deviation of ln(x)",
    "mean_log10": "Mean of log10(x)", "sd_log10": "Standard deviation of log10(x)",
    "skew_log10": "Sample skew of log10(x) (bias-corrected, no regional skew)",
}
UNITS = {"discharge": "Discharge (m\u00b3/s)", "rainfall": "Rainfall (mm)"}


def _num(ws, r, c, v, fmt=None):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        ws.write_string(r, c, "n/a", fmt)
    else:
        ws.write_number(r, c, float(v), fmt)


def write_workbook(path, res, unit_label, station="", data_type="discharge"):
    """res: the dict returned by floodfreq.analyze()."""
    x, n, grid, pp = res["x"], res["n"], res["grid"], res["pp"]
    fits, bands, st, rp = res["fits"], res["bands"], res["fit_statistics"], res["return_periods"]
    names = list(FITTERS)
    wb = xlsxwriter.Workbook(path)
    f = {
        "title": wb.add_format({"bold": True, "font_size": 14}),
        "b": wb.add_format({"bold": True}),
        "h": wb.add_format({"bold": True, "font_color": "white", "bg_color": "#0E6B83", "align": "center",
                            "valign": "vcenter", "text_wrap": True, "border": 1}),
        "n1": wb.add_format({"num_format": "#,##0.0"}), "n2": wb.add_format({"num_format": "0.00"}),
        "n3": wb.add_format({"num_format": "0.000"}), "n4": wb.add_format({"num_format": "0.0000"}),
        "w": wb.add_format({"text_wrap": True, "valign": "top"}),
    }
    q100 = res["q100"]

    # ---- Summary
    ws = wb.add_worksheet("Summary")
    ws.set_column(0, 0, 24); ws.set_column(1, 8, 14)
    ws.write(0, 0, "Flood frequency analysis", f["title"])
    meta = [("Station", station or "\u2014"), ("Data type", unit_label), ("Records (n)", n),
            ("Confidence level", f"{res['conf']:.0%}"), ("Bootstrap simulations", res["nboot"]),
            ("Plotting position", res["pp_method"]), ("Software", f"floodfreq {__version__}")]
    for i, (k, v) in enumerate(meta, 1):
        ws.write(i, 0, k, f["b"]); ws.write(i, 1, v)
    r = len(meta) + 2
    ws.write(r, 0, "Which model fits best?", f["b"]); r += 1
    for line in res["summary"]:
        ws.merge_range(r, 0, r, 8, line, f["w"])
        ws.set_row(r, 15 * max(1, math.ceil(len(line) / 105))); r += 1
    r += 1
    ws.write(r, 0, "Model ranking (by AIC; lowest is best)", f["b"]); r += 1
    for c, h in enumerate(["Rank", "Model", "AIC", "\u0394AIC", "Akaike weight", "KS", "RMSE", "100-yr estimate"]):
        ws.write(r, c, h, f["h"])
    r += 1
    order = st.assign(_a=st["AIC"].fillna(np.inf)).sort_values("_a", kind="stable")
    for rank, (_, s) in enumerate(order.iterrows(), 1):
        d = s["distribution"]
        ws.write(r, 0, rank if s["support_ok"] else "n/a"); ws.write(r, 1, LABEL[d])
        _num(ws, r, 2, s["AIC"], f["n2"]); _num(ws, r, 3, s["delta_AIC"], f["n2"]); _num(ws, r, 4, s["akaike_weight"], f["n2"])
        _num(ws, r, 5, s["KS"], f["n3"]); _num(ws, r, 6, s["RMSE"], f["n1"]); _num(ws, r, 7, q100[d], f["n1"])
        r += 1
    if res["mk"]:
        mk = res["mk"]; r += 1
        ws.write(r, 0, "Trend screening (Mann-Kendall)", f["b"]); r += 1
        for k, v, fm in [("S", mk["S"], None), ("Z", mk["Z"], f["n3"]), ("p-value (two-sided)", mk["p"], f["n4"]),
                         ("Kendall tau", mk["tau"], f["n3"]), ("Sen slope (per record step)", mk["sen_slope"], f["n3"]),
                         ("Result at 5%", mk["trend"], None)]:
            ws.write(r, 0, k); (_num(ws, r, 1, v, fm) if not isinstance(v, str) else ws.write(r, 1, v)); r += 1
        ws.write(r, 0, "Assumes the values are in chronological order; serial correlation is not tested.")

    # ---- Return Periods
    ws = wb.add_worksheet("Return Periods")
    ws.set_column(0, 0, 18); ws.set_column(1, 12, 13)
    ws.write(0, 0, f"Return period estimates ({unit_label}), {res['conf']:.0%} bootstrap intervals", f["title"])
    ws.merge_range(2, 0, 3, 0, "Return period (yr)", f["h"])
    for i, d in enumerate(names):
        ws.merge_range(2, 1 + 3 * i, 2, 3 + 3 * i, LABEL[d], f["h"])
        for j, h in enumerate(["Estimate", "Lower", "Upper"]):
            ws.write(3, 1 + 3 * i + j, h, f["h"])
    for k, T in enumerate(RETURN_PERIODS):
        ws.write_number(4 + k, 0, T)
        for i, d in enumerate(names):
            row = rp[(rp.distribution == d) & (rp.return_period_yr == T)].iloc[0]
            for j, key in enumerate(["estimate", "lower", "upper"]):
                _num(ws, 4 + k, 1 + 3 * i + j, row[key], f["n1"])
    ws.freeze_panes(4, 1)

    # ---- Fit Statistics
    ws = wb.add_worksheet("Fit Statistics")
    ws.set_column(0, 0, 18); ws.set_column(1, 9, 14); ws.set_column(10, 10, 44)
    heads = ["Model", "Parameters", "KS", "RMSE", "AIC", "\u0394AIC", "Akaike weight", "BIC", "AIC rank", "Support OK", "Note"]
    for c, h in enumerate(heads):
        ws.write(0, c, h, f["h"])
    for i, (_, s) in enumerate(st.iterrows(), 1):
        d = s["distribution"]
        ws.write(i, 0, LABEL[d]); ws.write(i, 1, fits[d].n_par)
        _num(ws, i, 2, s["KS"], f["n3"]); _num(ws, i, 3, s["RMSE"], f["n1"]); _num(ws, i, 4, s["AIC"], f["n2"])
        _num(ws, i, 5, s["delta_AIC"], f["n2"]); _num(ws, i, 6, s["akaike_weight"], f["n3"]); _num(ws, i, 7, s["BIC"], f["n2"])
        _num(ws, i, 8, s["rank_AIC"]); ws.write(i, 9, "yes" if s["support_ok"] else "no")
        ws.write(i, 10, "" if s["support_ok"] else "Observations outside fitted range; AIC/BIC undefined")
    ws.write(6, 0, "KS p-values are not reported: parameters are estimated from the same data, which makes them optimistic.")
    ws.write(7, 0, "RMSE compares sorted data with fitted quantiles at the plotting positions, so it depends on the plotting position used.")

    # ---- Parameters
    ws = wb.add_worksheet("Parameters")
    ws.set_column(0, 0, 18); ws.set_column(1, 1, 16); ws.set_column(2, 2, 14); ws.set_column(3, 3, 64)
    for c, h in enumerate(["Model", "Parameter", "Value", "Meaning"]):
        ws.write(0, c, h, f["h"])
    r = 1
    for d in names:
        for k, v in fits[d].params.items():
            ws.write(r, 0, LABEL[d]); ws.write(r, 1, k); _num(ws, r, 2, v, f["n4"]); ws.write(r, 3, PARAM_INFO.get(k, "")); r += 1

    # ---- Data
    ws = wb.add_worksheet("Data")
    ws.set_column(0, 1, 16); ws.set_column(2, 2, 3); ws.set_column(3, 13, 16)
    years = res.get("years")
    for c, h in [(0, "Year" if years else "Order"), (1, f"Value as entered ({unit_label})"), (3, "Rank (ascending)"),
                 (4, "Sorted value"), (5, "Non-exceedance probability"), (6, "Empirical return period (yr)")] + \
                [(7 + i, f"Fitted: {LABEL[d]}") for i, d in enumerate(names)]:
        ws.write(0, c, h, f["h"])
    ws.set_row(0, 45)
    xs = np.sort(x)
    for i in range(n):
        ws.write(i + 1, 0, years[i] if years else i + 1); ws.write_number(i + 1, 1, x[i])
        ws.write_number(i + 1, 3, i + 1); ws.write_number(i + 1, 4, xs[i], f["n1"])
        ws.write_number(i + 1, 5, pp[i], f["n4"]); ws.write_number(i + 1, 6, 1 / (1 - pp[i]), f["n2"])
        for j, d in enumerate(names):
            ws.write_number(i + 1, 7 + j, float(fits[d].ppf(pp[i])), f["n1"])
    ws.freeze_panes(1, 0)

    # ---- Curves (chart source)
    ws = wb.add_worksheet("Curves")
    ws.set_column(0, 12, 14)
    ws.write(0, 0, "Return period (yr)", f["h"])
    for i, d in enumerate(names):
        for j, h in enumerate(["", " lower", " upper"]):
            ws.write(0, 1 + 3 * i + j, LABEL[d] + h, f["h"])
    ng = len(grid)
    for k, T in enumerate(grid):
        ws.write_number(k + 1, 0, float(T), f["n2"])
        for i, d in enumerate(names):
            for j, v in enumerate([float(fits[d].ppf(1 - 1 / T)), bands[d][0][k], bands[d][1][k]]):
                ws.write_number(k + 1, 1 + 3 * i + j, float(v), f["n1"])

    # ---- Charts (native, editable)
    ws = wb.add_worksheet("Charts")
    ws.write(0, 0, "Charts are native Excel charts: click one to restyle it. Source data: 'Curves' and 'Data' sheets.", f["b"])

    def chart(title, dists, band):
        ch = wb.add_chart({"type": "scatter", "subtype": "straight"})
        for d in dists:
            c0 = 1 + 3 * names.index(d)
            ch.add_series({"name": LABEL[d], "categories": ["Curves", 1, 0, ng, 0], "values": ["Curves", 1, c0, ng, c0],
                           "line": {"color": COLORS[d], "width": 2.25}, "marker": {"type": "none"}})
            if band:
                for off, nm in ((1, f"{res['conf']:.0%} lower"), (2, f"{res['conf']:.0%} upper")):
                    ch.add_series({"name": nm, "categories": ["Curves", 1, 0, ng, 0], "values": ["Curves", 1, c0 + off, ng, c0 + off],
                                   "line": {"color": COLORS[d], "width": 1, "dash_type": "dash"}, "marker": {"type": "none"}})
        ch.add_series({"name": "Observed", "categories": ["Data", 1, 6, n, 6], "values": ["Data", 1, 4, n, 4],
                       "line": {"none": True}, "marker": {"type": "circle", "size": 5, "fill": {"color": "black"}, "border": {"color": "black"}}})
        ch.set_title({"name": f"{station} {title}".strip(), "name_font": {"size": 12}})
        ch.set_x_axis({"name": "Return period (years)", "log_base": 10, "min": 1, "max": 1000,
                       "major_gridlines": {"visible": True, "line": {"color": "#d9d9d9"}}})
        ch.set_y_axis({"name": unit_label, "major_gridlines": {"visible": True, "line": {"color": "#d9d9d9"}}})
        ch.set_legend({"position": "bottom"})
        ch.set_size({"width": 640, "height": 400})
        return ch

    spots = ["A3", "L3", "A24", "L24", "A45"]
    for spot, (title, dists, band) in zip(spots, [(f"{LABEL[d]}", [d], True) for d in names] + [("Distribution comparison", names, False)]):
        ws.insert_chart(spot, chart(title, dists, band))
    wb.close()
