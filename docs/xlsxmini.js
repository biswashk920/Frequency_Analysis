/* Minimal dependency-free .xlsx writer (tables + embedded PNG images). Works in browsers and Node. */
(function (root) {
"use strict";
const enc = new TextEncoder();
const esc = s => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
const clean = s => String(s).replace(/[\u0000-\u0008\u000B\u000C\u000E-\u001F]/g, "");
const colName = n => { let s = ""; n++; while (n > 0) { const r = (n - 1) % 26; s = String.fromCharCode(65 + r) + s; n = Math.floor((n - 1) / 26); } return s; };
const STYLE = { h: 1, b: 2, title: 3, n1: 4, n2: 5, n3: 6, w: 7, n4: 8 };
const NS = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"';
const HEAD = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
const xf = (n, f, fill, extra) => `<xf numFmtId="${n}" fontId="${f}" fillId="${fill}" borderId="0" xfId="0" applyNumberFormat="1" applyFont="1" applyFill="1"${extra ? ' applyAlignment="1">' + extra + "</xf>" : "/>"}`;
const STYLES = HEAD + `<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<numFmts count="4"><numFmt numFmtId="164" formatCode="#,##0.0"/><numFmt numFmtId="165" formatCode="0.00"/><numFmt numFmtId="166" formatCode="0.000"/><numFmt numFmtId="167" formatCode="0.0000"/></numFmts>
<fonts count="4"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="14"/><name val="Calibri"/></font></fonts>
<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF0E6B83"/><bgColor indexed="64"/></patternFill></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="9">${xf(0, 0, 0)}${xf(0, 1, 2, '<alignment horizontal="center" vertical="center" wrapText="1"/>')}${xf(0, 2, 0)}${xf(0, 3, 0)}${xf(164, 0, 0)}${xf(165, 0, 0)}${xf(166, 0, 0)}${xf(0, 0, 0, '<alignment vertical="top" wrapText="1"/>')}${xf(167, 0, 0)}</cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>`;

function sheetXml(sh, drawing) {
  let x = HEAD + `<worksheet ${NS}>`;
  if (sh.colWidths) x += "<cols>" + sh.colWidths.map((w, i) => `<col min="${i + 1}" max="${i + 1}" width="${w}" customWidth="1"/>`).join("") + "</cols>";
  x += "<sheetData>";
  sh.rows.forEach((row, r) => {
    if (!row || !row.length) return;
    const ht = sh.rowHeights && sh.rowHeights[r];
    x += `<row r="${r + 1}"${ht ? ` ht="${ht}" customHeight="1"` : ""}>`;
    row.forEach((c, ci) => {
      if (c === null || c === undefined || c === "") return;
      let v = c, s = 0; if (typeof c === "object") { v = c.v; s = STYLE[c.s] || 0; }
      if (v === "" || v === null || v === undefined) { if (s) x += `<c r="${colName(ci)}${r + 1}" s="${s}"/>`; return; }
      const ref = `${colName(ci)}${r + 1}`;
      if (typeof v === "number" && isFinite(v)) x += `<c r="${ref}" s="${s}"><v>${v}</v></c>`;
      else x += `<c r="${ref}" s="${s}" t="inlineStr"><is><t xml:space="preserve">${typeof v === "number" ? "n/a" : esc(clean(v))}</t></is></c>`;
    });
    x += "</row>";
  });
  x += "</sheetData>";
  if (sh.merges && sh.merges.length) x += `<mergeCells count="${sh.merges.length}">` + sh.merges.map(m => `<mergeCell ref="${m}"/>`).join("") + "</mergeCells>";
  if (drawing) x += '<drawing r:id="rId1"/>';
  return x + "</worksheet>";
}
function drawingXml(imgs) {
  let x = HEAD + '<xdr:wsDr xmlns:xdr="http://schemas.openxmlformats.org/drawingml/2006/spreadsheetDrawing" xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">';
  imgs.forEach((im, i) => {
    const cx = Math.round(im.w * 9525), cy = Math.round(im.h * 9525);
    x += `<xdr:oneCellAnchor><xdr:from><xdr:col>${im.col}</xdr:col><xdr:colOff>0</xdr:colOff><xdr:row>${im.row}</xdr:row><xdr:rowOff>0</xdr:rowOff></xdr:from><xdr:ext cx="${cx}" cy="${cy}"/>` +
      `<xdr:pic><xdr:nvPicPr><xdr:cNvPr id="${i + 2}" name="Picture ${i + 1}" descr="${esc(im.alt || "")}"/><xdr:cNvPicPr><a:picLocks noChangeAspect="1"/></xdr:cNvPicPr></xdr:nvPicPr>` +
      `<xdr:blipFill><a:blip r:embed="rId${i + 1}"/><a:stretch><a:fillRect/></a:stretch></xdr:blipFill><xdr:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="${cx}" cy="${cy}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></xdr:spPr></xdr:pic><xdr:clientData/></xdr:oneCellAnchor>`;
  });
  return x + "</xdr:wsDr>";
}
const CRC = (() => { const t = new Uint32Array(256); for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xEDB88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; } return t; })();
const crc32 = b => { let c = 0xFFFFFFFF; for (let i = 0; i < b.length; i++) c = CRC[(c ^ b[i]) & 255] ^ (c >>> 8); return (c ^ 0xFFFFFFFF) >>> 0; };
const u16 = v => [v & 255, (v >>> 8) & 255], u32 = v => [v & 255, (v >>> 8) & 255, (v >>> 16) & 255, (v >>> 24) & 255];
function zip(files) { // stored (uncompressed) zip
  const parts = [], cen = []; let off = 0;
  for (const f of files) {
    const nm = enc.encode(f.name), crc = crc32(f.data), sz = f.data.length;
    const lh = new Uint8Array([0x50, 0x4b, 3, 4, ...u16(20), ...u16(0), ...u16(0), ...u16(0), ...u16(0x21), ...u32(crc), ...u32(sz), ...u32(sz), ...u16(nm.length), ...u16(0)]);
    parts.push(lh, nm, f.data);
    cen.push(new Uint8Array([0x50, 0x4b, 1, 2, ...u16(20), ...u16(20), ...u16(0), ...u16(0), ...u16(0), ...u16(0x21), ...u32(crc), ...u32(sz), ...u32(sz), ...u16(nm.length), ...u16(0), ...u16(0), ...u16(0), ...u16(0), ...u32(0), ...u32(off)]), nm);
    off += lh.length + nm.length + sz;
  }
  const cdSize = cen.reduce((s, a) => s + a.length, 0);
  const end = new Uint8Array([0x50, 0x4b, 5, 6, ...u16(0), ...u16(0), ...u16(files.length), ...u16(files.length), ...u32(cdSize), ...u32(off), ...u16(0)]);
  const all = [...parts, ...cen, end], out = new Uint8Array(all.reduce((s, a) => s + a.length, 0)); let p = 0;
  for (const a of all) { out.set(a, p); p += a.length; }
  return out;
}
/* sheets: [{name, rows:[[cell|{v,s}]], colWidths, merges:["A1:C1"], rowHeights:{rowIndex:pts}, images:[{data:Uint8Array,col,row,w,h,alt}]}] */
function build(sheets) {
  const F = [], T = s => enc.encode(s), REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";
  let media = 0, drawings = 0, ct = "", wbRels = "", wbSheets = "";
  sheets.forEach((sh, i) => {
    const n = i + 1, imgs = sh.images || [];
    ct += `<Override PartName="/xl/worksheets/sheet${n}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>`;
    wbRels += `<Relationship Id="rId${n}" Type="${REL}/worksheet" Target="worksheets/sheet${n}.xml"/>`;
    wbSheets += `<sheet name="${esc(sh.name.slice(0, 31))}" sheetId="${n}" r:id="rId${n}"/>`;
    F.push({ name: `xl/worksheets/sheet${n}.xml`, data: T(sheetXml(sh, imgs.length > 0)) });
    if (imgs.length) {
      const d = ++drawings; let rels = "";
      imgs.forEach((im, k) => { const m = ++media; F.push({ name: `xl/media/image${m}.png`, data: im.data }); rels += `<Relationship Id="rId${k + 1}" Type="${REL}/image" Target="../media/image${m}.png"/>`; });
      F.push({ name: `xl/worksheets/_rels/sheet${n}.xml.rels`, data: T(HEAD + `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/drawing" Target="../drawings/drawing${d}.xml"/></Relationships>`) });
      F.push({ name: `xl/drawings/drawing${d}.xml`, data: T(drawingXml(imgs)) });
      F.push({ name: `xl/drawings/_rels/drawing${d}.xml.rels`, data: T(HEAD + `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">${rels}</Relationships>`) });
      ct += `<Override PartName="/xl/drawings/drawing${d}.xml" ContentType="application/vnd.openxmlformats-officedocument.drawing+xml"/>`;
    }
  });
  const ns = sheets.length;
  F.unshift({ name: "[Content_Types].xml", data: T(HEAD + `<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>${ct}</Types>`) });
  F.push({ name: "_rels/.rels", data: T(HEAD + `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/officeDocument" Target="xl/workbook.xml"/></Relationships>`) });
  F.push({ name: "xl/workbook.xml", data: T(HEAD + `<workbook ${NS}><sheets>${wbSheets}</sheets></workbook>`) });
  F.push({ name: "xl/_rels/workbook.xml.rels", data: T(HEAD + `<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">${wbRels}<Relationship Id="rId${ns + 1}" Type="${REL}/styles" Target="styles.xml"/></Relationships>`) });
  F.push({ name: "xl/styles.xml", data: T(STYLES) });
  return zip(F);
}
const dataUrlToBytes = u => { const b = atob(u.split(",")[1]), a = new Uint8Array(b.length); for (let i = 0; i < b.length; i++) a[i] = b.charCodeAt(i); return a; };
const X = { build, dataUrlToBytes };
if (typeof module !== "undefined" && module.exports) module.exports = X; else root.XLSXMINI = X;
})(typeof self !== "undefined" ? self : this);
