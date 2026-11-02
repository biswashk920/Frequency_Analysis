// Runs the browser implementation (docs/floodfreq.js, docs/xlsxmini.js) on datasets given as JSON on stdin.
const fs = require("fs"), path = require("path");
const FF = require(path.join(__dirname, "..", "docs", "floodfreq.js"));
const XL = require(path.join(__dirname, "..", "docs", "xlsxmini.js"));
const input = JSON.parse(fs.readFileSync(0, "utf8"));
(async () => {
  const out = { datasets: {}, parse: {} };
  for (const [name, x] of Object.entries(input.datasets)) {
    const r = await FF.analyze(x, { nboot: 60, trend: true });
    out.datasets[name] = {
      q: Object.fromEntries(FF.DIST.map(d => [d, FF.RP.map((_, k) => r.res[d].q[r.idx[k]])])),
      stats: r.rows.map(s => ({ ...s })), mk: r.mk, summary: r.summary, q100: r.q100,
    };
    if (name === input.xlsx_dataset) {
      const sheets = FF.workbookSheets(r, { station: "JS test", unitLabel: "Discharge (m\u00b3/s)" });
      const png = Uint8Array.from(Buffer.from(input.png_b64, "base64"));
      sheets.push({ name: "Figures", rows: [["Figures are embedded as images."]], images: [{ data: png, col: 0, row: 2, w: 40, h: 40, alt: "test" }, { data: png, col: 8, row: 2, w: 40, h: 40, alt: "test2" }] });
      fs.writeFileSync(input.xlsx_path, XL.build(sheets));
    }
  }
  for (const s of input.parse) { try { out.parse[s] = FF.parseNumbers(s).values; } catch (e) { out.parse[s] = "error"; } }
  out.errors = {};
  for (const [k, v] of Object.entries(input.invalid)) { try { await FF.analyze(v.map(t => t === null ? NaN : t), { nboot: 50 }); out.errors[k] = null; } catch (e) { out.errors[k] = e.message; } }
  console.log(JSON.stringify(out, (k, v) => (typeof v === "number" && !isFinite(v) ? null : v)));
})();
