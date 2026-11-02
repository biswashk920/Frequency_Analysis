/* floodfreq core (browser + Node). Mirrors the Python package: same fits, same statistics.
   Gumbel/GEV: L-moments. Log Normal: moments of ln x. LP3: moments of log10 x with sample skew and the
   exact Pearson III distribution (incomplete gamma), as in SciPy. */
(function (root) {
"use strict";
const RP = [2, 5, 10, 25, 50, 100, 200, 500], DIST = ["Gumbel", "GEV", "Lognormal", "LP3"];
const LABEL = { Gumbel: "Gumbel", GEV: "GEV", Lognormal: "Log Normal", LP3: "Log Pearson III" };
const VERSION = "0.2.0", EULER = 0.5772156649015329, SQ2PI = Math.sqrt(2 * Math.PI);
const mean = a => a.reduce((s, v) => s + v, 0) / a.length;
const sd = a => { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); };
const sortNum = a => [...a].sort((p, q) => p - q);

/* ---------- special functions ---------- */
function nCdf(x) { // Marsaglia (2004), ~1e-15 absolute
  if (x < -8) return 0; if (x > 8) return 1;
  let s = x, t = 0, b = x, i = 1; const q = x * x;
  while (s !== t) { t = s; i += 2; b *= q / i; s = t + b; }
  return 0.5 + s * Math.exp(-0.5 * q - 0.91893853320467274178);
}
function nInv(p) { // Acklam + one Halley step
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239],
    b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572],
    c = [-.007784894002430293, -.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783],
    d = [.007784695709041462, .3224671290700398, 2.445134137286659, 3.754408661907416];
  let q, r, x;
  if (p < .02425) { q = Math.sqrt(-2 * Math.log(p)); x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else if (p > .97575) { q = Math.sqrt(-2 * Math.log(1 - p)); x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else { q = p - .5; r = q * q; x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); }
  if (p > 1e-8 && p < 1 - 1e-8) { const e = nCdf(x) - p, u = e * SQ2PI * Math.exp(x * x / 2); x -= u / (1 + x * u / 2); }
  return x;
}
const LC = [.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
function lgam(z) {
  if (z < .5) return Math.log(Math.PI / Math.abs(Math.sin(Math.PI * z))) - lgam(1 - z);
  z -= 1; let x = LC[0]; for (let i = 1; i < 9; i++) x += LC[i] / (z + i);
  const t = z + 7.5; return .5 * Math.log(2 * Math.PI) + (z + .5) * Math.log(t) - t + Math.log(x);
}
const gam = z => z < .5 ? Math.PI / (Math.sin(Math.PI * z) * gam(1 - z)) : Math.exp(lgam(z));
// Regularised lower incomplete gamma P(a, x). For a > 1e4 (|skew| < 0.02) Wilson-Hilferty is used (error << 1e-4).
function gammaP(a, x) {
  if (x <= 0) return 0;
  if (a > 1e4) return nCdf((Math.pow(x / a, 1 / 3) - (1 - 1 / (9 * a))) * 3 * Math.sqrt(a));
  const gln = lgam(a);
  if (x < a + 1) {
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < 20000; n++) { ap += 1; del *= x / ap; sum += del; if (Math.abs(del) < Math.abs(sum) * 1e-15) break; }
    return sum * Math.exp(-x + a * Math.log(x) - gln);
  }
  let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
  for (let i = 1; i < 20000; i++) {
    const an = -i * (i - a); b += 2; d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
    c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300; d = 1 / d; const del = d * c; h *= del; if (Math.abs(del - 1) < 1e-15) break;
  }
  return 1 - Math.exp(-x + a * Math.log(x) - gln) * h;
}
function gammaPinv(p, a) {
  if (p <= 0) return 0; if (p >= 1) return Infinity;
  const z = nInv(p); let x = a * Math.pow(1 - 1 / (9 * a) + z / (3 * Math.sqrt(a)), 3);
  if (a > 1e4) return Math.max(x, 0);
  if (!(x > 0) || a < 1) x = Math.exp((Math.log(p) + lgam(a + 1)) / a);
  const lg = lgam(a); let lo = 0, hi = Infinity;
  for (let it = 0; it < 200; it++) {
    const f = gammaP(a, x) - p; if (f > 0) hi = Math.min(hi, x); else lo = Math.max(lo, x);
    let xn = x - f / Math.exp(-x + (a - 1) * Math.log(x) - lg);
    if (!(xn > lo && xn < hi)) xn = isFinite(hi) ? (lo + hi) / 2 : x * 2;
    if (Math.abs(xn - x) < 1e-14 * Math.max(1, x)) { x = xn; break; } x = xn;
  }
  return x;
}

/* ---------- fits: {name, par, nPar, ppf, cdf, logpdf} ---------- */
function lm(x) {
  const s = sortNum(x), n = s.length; let b1 = 0, b2 = 0;
  for (let i = 1; i <= n; i++) { b1 += (i - 1) / (n - 1) * s[i - 1]; b2 += (i - 1) * (i - 2) / ((n - 1) * (n - 2)) * s[i - 1]; }
  const b0 = mean(s); b1 /= n; b2 /= n; const l2 = 2 * b1 - b0; return [b0, l2, (6 * b2 - 6 * b1 + b0) / l2];
}
const gum = (xi, a) => ({
  ppf: F => xi - a * Math.log(-Math.log(F)), cdf: v => Math.exp(-Math.exp(-(v - xi) / a)),
  logpdf: v => { const z = (v - xi) / a; return -Math.log(a) - z - Math.exp(-z); }
});
const gevf = (xi, a, k) => ({
  ppf: F => xi + a / k * (1 - Math.pow(-Math.log(F), k)),
  cdf: v => { const t = 1 - k * (v - xi) / a; return t <= 0 ? (k > 0 ? 1 : 0) : Math.exp(-Math.pow(t, 1 / k)); },
  logpdf: v => { const t = 1 - k * (v - xi) / a; return t <= 0 ? -Infinity : -Math.log(a) + (1 / k - 1) * Math.log(t) - Math.pow(t, 1 / k); }
});
function fitGumbel(x) {
  const [l1, l2] = lm(x), a = l2 / Math.LN2, xi = l1 - EULER * a;
  return { name: "Gumbel", par: { location: xi, scale: a }, nPar: 2, ...gum(xi, a) };
}
function fitGEV(x) {
  const [l1, l2, t3] = lm(x), c = 2 / (3 + t3) - Math.LN2 / Math.log(3);
  const k = Math.max(-.5, Math.min(.5, 7.859 * c + 2.9554 * c * c));
  if (Math.abs(k) < 1e-6) { const g = fitGumbel(x); return { ...g, name: "GEV", par: { ...g.par, shape_k: 0 }, nPar: 3 }; }
  const a = l2 * k / ((1 - Math.pow(2, -k)) * gam(1 + k)), xi = l1 - a * (1 - gam(1 + k)) / k;
  return { name: "GEV", par: { location: xi, scale: a, shape_k: k }, nPar: 3, ...gevf(xi, a, k) };
}
function fitLN(x) {
  const y = x.map(v => Math.log(v)), m = mean(y), s = sd(y);
  return {
    name: "Lognormal", par: { mu_ln: m, sigma_ln: s }, nPar: 2, ppf: F => Math.exp(m + s * nInv(F)),
    cdf: v => v <= 0 ? 0 : nCdf((Math.log(v) - m) / s),
    logpdf: v => { const z = (Math.log(v) - m) / s; return -Math.log(v) - Math.log(s) - Math.log(SQ2PI) - 0.5 * z * z; }
  };
}
function fitLP3(x) {
  const y = x.map(v => Math.log10(v)), n = y.length, m = mean(y), s = sd(y);
  const g = n / ((n - 1) * (n - 2)) * y.reduce((t, v) => t + ((v - m) / s) ** 3, 0);
  const norm = Math.abs(g) < 1e-6, a = 4 / (g * g), sa = Math.sqrt(a), pos = g > 0;
  const zq = F => norm ? nInv(F) : pos ? (gammaPinv(F, a) - a) / sa : -(gammaPinv(1 - F, a) - a) / sa;
  const zc = z => { if (norm) return nCdf(z); if (pos) { const t = a + sa * z; return t <= 0 ? 0 : gammaP(a, t); } const t = a - sa * z; return t <= 0 ? 1 : 1 - gammaP(a, t); };
  const zl = z => { if (norm) return -z * z / 2 - Math.log(SQ2PI); const t = pos ? a + sa * z : a - sa * z; return t <= 0 ? -Infinity : Math.log(sa) + (a - 1) * Math.log(t) - t - lgam(a); };
  return {
    name: "LP3", par: { mean_log10: m, sd_log10: s, skew_log10: g }, nPar: 3,
    ppf: F => Math.pow(10, m + s * zq(F)), cdf: v => zc((Math.log10(v) - m) / s),
    logpdf: v => zl((Math.log10(v) - m) / s) - Math.log(s) - Math.log(v * Math.LN10)
  };
}
const FIT = { Gumbel: fitGumbel, GEV: fitGEV, Lognormal: fitLN, LP3: fitLP3 };

/* ---------- inputs ---------- */
const THOUS = /^\d{1,2}(,\d{3})+(\.\d+)?$/;
function parseNumbers(text) {
  const values = [], notes = [];
  for (const tok of text.trim().split(/[\s;]+/)) {
    if (!tok) continue; let parts;
    if (THOUS.test(tok)) { notes.push(`read '${tok}' as ${tok.replace(/,/g, "")}`); parts = [tok.replace(/,/g, "")]; }
    else parts = tok.split(",").filter(Boolean);
    for (const p of parts) { const v = Number(p); if (Number.isNaN(v)) throw new Error(`Could not read '${p}' as a number.`); values.push(v); }
  }
  return { values, notes };
}
function validate(x) {
  if (x.length < 5) throw new Error("Need at least 5 values.");
  if (!x.every(v => typeof v === "number" && isFinite(v))) throw new Error("Values must be finite numbers (found NaN or infinity).");
  if (x.some(v => v <= 0)) throw new Error("All values must be positive (Log Normal and LP3 use logarithms).");
  if (Math.max(...x) === Math.min(...x)) throw new Error("All values are identical; a distribution cannot be fitted.");
  return x;
}

/* ---------- statistics ---------- */
const PPA = { weibull: 0, cunnane: .4, gringorten: .44 };
const plottingPositions = (n, m = "weibull") => Array.from({ length: n }, (_, i) => (i + 1 - PPA[m]) / (n + 1 - 2 * PPA[m]));
const quantiles = (f, T) => T.map(t => f.ppf(1 - 1 / t));
function gof(x, f, pp) {
  const s = sortNum(x), n = s.length; let D = 0;
  s.forEach((v, i) => { const F = f.cdf(v); D = Math.max(D, (i + 1) / n - F, F - i / n); });
  const rmse = Math.sqrt(mean(s.map((v, i) => (v - f.ppf(pp[i])) ** 2)));
  let ll = 0, ok = true;
  for (const v of s) { const l = f.logpdf(v); if (!isFinite(l)) { ok = false; break; } ll += l; }
  return { KS: D, RMSE: rmse, AIC: ok ? 2 * f.nPar - 2 * ll : NaN, BIC: ok ? f.nPar * Math.log(n) - 2 * ll : NaN, support_ok: ok };
}
function akaike(aic) {
  const fin = aic.filter(isFinite), mn = Math.min(...fin), delta = aic.map(a => isFinite(a) ? a - mn : NaN);
  const raw = delta.map(d => isFinite(d) ? Math.exp(-d / 2) : 0), tot = raw.reduce((s, v) => s + v, 0);
  return { delta, w: raw.map(v => tot ? v / tot : 0) };
}
function mannKendall(x, alpha = 0.05) {
  const n = x.length; if (n < 4) throw new Error("The Mann-Kendall test needs at least 4 values.");
  let S = 0; const slopes = [];
  for (let i = 0; i < n - 1; i++) for (let j = i + 1; j < n; j++) { const d = x[j] - x[i]; S += Math.sign(d); slopes.push(d / (j - i)); }
  const cnt = new Map(); x.forEach(v => cnt.set(v, (cnt.get(v) || 0) + 1));
  let tie = 0; cnt.forEach(t => { if (t > 1) tie += t * (t - 1) * (2 * t + 5); });
  const v = (n * (n - 1) * (2 * n + 5) - tie) / 18; if (v <= 0) throw new Error("All values are tied; the Mann-Kendall test is undefined.");
  const z = S === 0 ? 0 : (S - Math.sign(S)) / Math.sqrt(v), p = 2 * nCdf(-Math.abs(z));
  const ss = sortNum(slopes), h = ss.length, sen = h % 2 ? ss[(h - 1) / 2] : (ss[h / 2 - 1] + ss[h / 2]) / 2;
  return { n, S, var_S: v, Z: z, p, alpha, tau: S / (n * (n - 1) / 2), sen_slope: sen, trend: p >= alpha ? "none" : z > 0 ? "increasing" : "decreasing" };
}

/* ---------- summary (same wording as floodfreq/report.py) ---------- */
const fmt = v => v >= 1000 ? v.toFixed(0) : v >= 100 ? v.toFixed(1) : v.toFixed(2);
const pfmt = p => p < 0.001 ? "< 0.001" : p.toFixed(3);
const names = rs => { const n = rs.map(r => LABEL[r.distribution]); return n.length === 1 ? n[0] : n.slice(0, -1).join(", ") + " and " + n[n.length - 1]; };
function bestModelSummary(rows, q100, n, gevK, mk) {
  const L = [], ok = rows.filter(r => isFinite(r.AIC)).sort((a, b) => a.AIC - b.AIC);
  if (ok.length) {
    const b = ok[0], close = ok.slice(1).filter(r => r.delta_AIC < 2);
    L.push(`Best fit by AIC: ${LABEL[b.distribution]} (Akaike weight ${b.akaike_weight.toFixed(2)}).`);
    if (close.length) L.push(`${names(close)} ${close.length === 1 ? "is" : "are"} within \u0394AIC < 2 of it, so the data cannot separate ${close.length === 1 ? "the two" : "them"}; compare the upper tails and use hydrological judgement.`);
    else L.push(`No other model is within \u0394AIC = 2, so AIC favours ${LABEL[b.distribution]} over the rest.`);
    const ks = rows.reduce((a, c) => c.KS < a.KS ? c : a).distribution, rm = rows.reduce((a, c) => c.RMSE < a.RMSE ? c : a).distribution;
    if (ks === rm && rm === b.distribution) L.push(`KS, RMSE and AIC all point to ${LABEL[ks]}.`);
    else L.push(`Other criteria: lowest KS = ${LABEL[ks]}, lowest RMSE = ${LABEL[rm]}. These can disagree with AIC, and RMSE depends on the plotting position.`);
  } else L.push("AIC could not be computed for any model (observations fall outside every fitted range).");
  const bad = rows.filter(r => !r.support_ok);
  if (bad.length) L.push(`AIC/BIC not available for ${names(bad)}: some observations fall outside the range of the fitted distribution.`);
  const ds = Object.keys(q100).filter(d => isFinite(q100[d]));
  if (ds.length > 1) {
    const lo = ds.reduce((a, c) => q100[c] < q100[a] ? c : a), hi = ds.reduce((a, c) => q100[c] > q100[a] ? c : a);
    L.push(`The 100-year estimate ranges from ${fmt(q100[lo])} (${LABEL[lo]}) to ${fmt(q100[hi])} (${LABEL[hi]}) across the models, a spread of ${((q100[hi] / q100[lo] - 1) * 100).toFixed(0)}%. This model-to-model difference is separate from the bootstrap bands, which only reflect sampling noise within each fitted model.`);
  }
  if (n < 30) L.push(`n = ${n} is below 30, so the 100\u2013500 year estimates are heavy extrapolation.`);
  if (gevK !== undefined && Math.abs(gevK) > 0.4) L.push(`GEV shape k = ${gevK.toFixed(2)} is near the limit of the L-moment approximation.`);
  if (mk) L.push(mk.trend !== "none"
    ? `Mann-Kendall found a significant ${mk.trend} trend (p = ${pfmt(mk.p)}); the stationarity assumption behind all four fits is questionable.`
    : `Mann-Kendall found no significant trend (p = ${pfmt(mk.p)}); serial correlation was not tested.`);
  return L;
}

/* ---------- bootstrap + analysis ---------- */
function rng(seed) { return () => { seed |= 0; seed = seed + 0x6D2B79F5 | 0; let t = Math.imul(seed ^ seed >>> 15, 1 | seed); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
const qtile = (s, p) => { const h = p * (s.length - 1), i = Math.floor(h), f = h - i; return i + 1 < s.length ? s[i] + f * (s[i + 1] - s[i]) : s[i]; };
function bootstrap(f, fit, n, T, B, conf, seed) {
  const r = rng(seed), out = [];
  for (let b = 0; b < B; b++) {
    const s = Array.from({ length: n }, () => f.ppf(1e-9 + r() * (1 - 2e-9)));
    if (s.some(v => !isFinite(v))) continue;
    try { const q = quantiles(fit(s), T); if (q.every(isFinite)) out.push(q); } catch (e) { /* skip failed refit */ }
  }
  if (out.length < B / 2) throw new Error("Too many failed bootstrap fits; the data may be unsuitable for this distribution.");
  const a = (1 - conf) / 2, cols = T.map((_, j) => sortNum(out.map(o => o[j])));
  return { lo: cols.map(c => qtile(c, a)), hi: cols.map(c => qtile(c, 1 - a)) };
}
const tick = () => new Promise(r => setTimeout(r, 0));
async function analyze(xIn, opts = {}) {
  const o = Object.assign({ conf: .95, nboot: 1000, pp: "weibull", trend: false, seed: 42, onProgress: null }, opts);
  const x = validate(Array.from(xIn)), n = x.length;
  if (!(o.conf > 0 && o.conf < 1)) throw new Error("Confidence level must be between 0 and 1.");
  const gm = Array.from({ length: 45 }, (_, i) => i === 44 ? 500 : 1.02 * Math.pow(500 / 1.02, i / 44));
  const grid = [...new Set([...gm, ...RP])].sort((a, b) => a - b), idx = RP.map(t => grid.indexOf(t));
  const pp = plottingPositions(n, o.pp), srt = sortNum(x), res = {}, rows = [];
  for (const d of DIST) {
    if (o.onProgress) o.onProgress(d); await tick();
    const f = FIT[d](x), bd = bootstrap(f, FIT[d], n, grid, o.nboot, o.conf, o.seed);
    res[d] = { f, q: quantiles(f, grid), lo: bd.lo, hi: bd.hi, stat: gof(x, f, pp) };
    rows.push({ distribution: d, ...res[d].stat });
  }
  const ak = akaike(rows.map(r => r.AIC));
  rows.forEach((r, i) => { r.delta_AIC = ak.delta[i]; r.akaike_weight = ak.w[i]; res[r.distribution].stat = r; });
  const q100 = {}; DIST.forEach(d => q100[d] = res[d].q[idx[5]]);
  const mk = o.trend ? mannKendall(x) : null;
  return { x, n, conf: o.conf, nboot: o.nboot, pp_method: o.pp, grid, idx, pp, srt, res, rows, q100, mk,
    summary: bestModelSummary(rows, q100, n, res.GEV.f.par.shape_k, mk) };
}

/* ---------- workbook layout (same sheets as the Python/XlsxWriter export, figures added by the page) ---------- */
function workbookSheets(r, meta) {
  const unit = meta.unitLabel, S = [], conf = `${Math.round(r.conf * 100)}%`;
  const num = (v, s) => isFinite(v) ? { v, s } : "n/a";
  const order = [...r.rows].sort((a, b) => (isFinite(a.AIC) ? a.AIC : Infinity) - (isFinite(b.AIC) ? b.AIC : Infinity));
  // Summary
  const sm = [[{ v: "Flood frequency analysis", s: "title" }], ["Station", meta.station || "\u2014"], ["Data type", unit], ["Records (n)", r.n],
    ["Confidence level", conf], ["Bootstrap simulations", r.nboot], ["Plotting position", r.pp_method], ["Software", "floodfreq " + VERSION], [],
    [{ v: "Which model fits best?", s: "b" }]];
  const merges = [], heights = {};
  r.summary.forEach(line => { sm.push([{ v: line, s: "w" }]); merges.push(`A${sm.length}:I${sm.length}`); heights[sm.length - 1] = 15 * Math.max(1, Math.ceil(line.length / 105)); });
  sm.push([], [{ v: "Model ranking (by AIC; lowest is best)", s: "b" }],
    ["Rank", "Model", "AIC", "\u0394AIC", "Akaike weight", "KS", "RMSE", "100-yr estimate"].map(h => ({ v: h, s: "h" })));
  order.forEach((s, i) => sm.push([s.support_ok ? i + 1 : "n/a", LABEL[s.distribution], num(s.AIC, "n2"), num(s.delta_AIC, "n2"), num(s.akaike_weight, "n2"),
    num(s.KS, "n3"), num(s.RMSE, "n1"), num(r.q100[s.distribution], "n1")]));
  if (r.mk) { const m = r.mk; sm.push([], [{ v: "Trend screening (Mann-Kendall)", s: "b" }], ["S", m.S], ["Z", num(m.Z, "n3")], ["p-value (two-sided)", num(m.p, "n3")],
    ["Kendall tau", num(m.tau, "n3")], ["Sen slope (per record step)", num(m.sen_slope, "n3")], ["Result at 5%", m.trend],
    ["Assumes the values are in chronological order; serial correlation is not tested."]); }
  S.push({ name: "Summary", rows: sm, colWidths: [24, 14, 14, 14, 14, 14, 14, 14, 14], merges, rowHeights: heights });
  // Return Periods
  const rp = [[{ v: `Return period estimates (${unit}), ${conf} bootstrap intervals`, s: "title" }], [],
    [{ v: "Return period (yr)", s: "h" }, ...DIST.flatMap(d => [{ v: LABEL[d], s: "h" }, { v: "", s: "h" }, { v: "", s: "h" }])],
    [{ v: "", s: "h" }, ...DIST.flatMap(() => ["Estimate", "Lower", "Upper"].map(h => ({ v: h, s: "h" })))]];
  const rpm = ["A3:A4"]; DIST.forEach((_, i) => { const c = 1 + 3 * i, L = k => String.fromCharCode(65 + k); rpm.push(`${L(c)}3:${L(c + 2)}3`); });
  RP.forEach((T, k) => rp.push([T, ...DIST.flatMap(d => [num(r.res[d].q[r.idx[k]], "n1"), num(r.res[d].lo[r.idx[k]], "n1"), num(r.res[d].hi[r.idx[k]], "n1")])]));
  S.push({ name: "Return Periods", rows: rp, colWidths: [18, ...Array(12).fill(13)], merges: rpm });
  // Fit statistics
  const fs = [["Model", "Parameters", "KS", "RMSE", "AIC", "\u0394AIC", "Akaike weight", "BIC", "AIC rank", "Support OK", "Note"].map(h => ({ v: h, s: "h" }))];
  r.rows.forEach(s => fs.push([LABEL[s.distribution], r.res[s.distribution].f.nPar, num(s.KS, "n3"), num(s.RMSE, "n1"), num(s.AIC, "n2"), num(s.delta_AIC, "n2"),
    num(s.akaike_weight, "n3"), num(s.BIC, "n2"), s.support_ok ? order.indexOf(s) + 1 : "n/a", s.support_ok ? "yes" : "no", s.support_ok ? "" : "Observations outside fitted range; AIC/BIC undefined"]));
  fs.push([], ["KS p-values are not reported: parameters are estimated from the same data, which makes them optimistic."],
    ["RMSE compares sorted data with fitted quantiles at the plotting positions, so it depends on the plotting position used."]);
  S.push({ name: "Fit Statistics", rows: fs, colWidths: [18, ...Array(9).fill(14), 44] });
  // Parameters
  const INFO = { location: "Location \u03be (same units as the data)", scale: "Scale \u03b1 (same units as the data)", shape_k: "Shape k (k < 0: heavy upper tail; k > 0: bounded upper tail)",
    mu_ln: "Mean of ln(x)", sigma_ln: "Standard deviation of ln(x)", mean_log10: "Mean of log10(x)", sd_log10: "Standard deviation of log10(x)",
    skew_log10: "Sample skew of log10(x) (bias-corrected, no regional skew)" };
  const pr = [["Model", "Parameter", "Value", "Meaning"].map(h => ({ v: h, s: "h" }))];
  DIST.forEach(d => Object.entries(r.res[d].f.par).forEach(([k, v]) => pr.push([LABEL[d], k, num(v, "n4"), INFO[k] || ""])));
  S.push({ name: "Parameters", rows: pr, colWidths: [18, 16, 14, 64] });
  // Data
  const da = [[meta.years ? "Year" : "Order", `Value as entered (${unit})`, "", "Rank (ascending)", "Sorted value", "Non-exceedance probability", "Empirical return period (yr)",
    ...DIST.map(d => "Fitted: " + LABEL[d])].map(h => h === "" ? "" : { v: h, s: "h" })];
  for (let i = 0; i < r.n; i++) da.push([meta.years ? meta.years[i] : i + 1, r.x[i], "", i + 1, { v: r.srt[i], s: "n1" }, { v: r.pp[i], s: "n4" }, { v: 1 / (1 - r.pp[i]), s: "n2" },
    ...DIST.map(d => ({ v: r.res[d].f.ppf(r.pp[i]), s: "n1" }))]);
  S.push({ name: "Data", rows: da, colWidths: [16, 16, 3, ...Array(11).fill(16)] });
  // Curves
  const cu = [[{ v: "Return period (yr)", s: "h" }, ...DIST.flatMap(d => [LABEL[d], LABEL[d] + " lower", LABEL[d] + " upper"].map(h => ({ v: h, s: "h" })))]];
  r.grid.forEach((T, k) => cu.push([{ v: T, s: "n2" }, ...DIST.flatMap(d => [{ v: r.res[d].q[k], s: "n1" }, { v: r.res[d].lo[k], s: "n1" }, { v: r.res[d].hi[k], s: "n1" }])]));
  S.push({ name: "Curves", rows: cu, colWidths: Array(13).fill(14) });
  return S;
}

const FF = { RP, DIST, LABEL, VERSION, FIT, nInv, nCdf, gammaP, gammaPinv, parseNumbers, validate, plottingPositions, quantiles, gof, akaike,
  mannKendall, bestModelSummary, analyze, workbookSheets, fmt };
if (typeof module !== "undefined" && module.exports) module.exports = FF; else root.FF = FF;
})(typeof self !== "undefined" ? self : this);
