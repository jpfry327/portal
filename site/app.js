/* Chordoma Atlas — static SPA on DuckDB-WASM.
 * Data URLs come from config.js (local ./data in dev; object storage in prod).
 */
import * as duckdb from "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.29.0/+esm";
import { DATA_BASE, ZENODO_MODE } from "./config.js";

// DuckDB-WASM's HTTP reader only engages on absolute http(s) URLs; a relative
// "./data/x.parquet" is treated as a (empty) virtual-filesystem glob. Resolve
// the configured base against the page URL at runtime.
const ABS_BASE = /^https?:\/\//.test(DATA_BASE)
  ? DATA_BASE.replace(/\/$/, "")
  : new URL(DATA_BASE.replace(/\/$/, "") + "/", location.href).href.replace(/\/$/, "");
// Published builds rename the Parquet files to opaque tokens (see
// portal/obfuscate_names.py) so a stray URL does not advertise its contents.
// datamap.json carries logical name -> token; when it is absent (local dev)
// the readable names are used unchanged.
let DATAMAP = {};
async function loadDataMap() {
  try {
    const r = await fetch(new URL("datamap.json", location.href), { cache: "no-store" });
    if (r.ok) DATAMAP = await r.json();
  } catch { /* local dev without a map: fall through to readable names */ }
}
const pqUrl = (name) =>
  `${ABS_BASE}/${DATAMAP[name] ?? name}.parquet${ZENODO_MODE ? "/content" : ""}`;

const S = globalThis.ChordStats;
const COLORS = { C1: "#440154", C2: "#31688E", C3: "#35B779", C4: "#FDE725" };
const ORDER = ["C1", "C2", "C3", "C4"];
const SUBTYPE_INFO = {
  C1: ["C1 hypermetabolic", "Loss of notochordal differentiation; nucleotide biosynthesis. Poor prognosis."],
  C2: ["C2 neural proliferative", "MYCN-driven, highly proliferative; neural/GABAergic programs. Poor prognosis."],
  C3: ["C3 immunogenic", "Immune and inflammatory infiltration. Intermediate outcomes."],
  C4: ["C4 LYST-activated", "Hypomethylated enhancer landscape; LYST super-enhancer. Best prognosis."],
};
const COHORTS = { discovery: "Discovery", zhang: "Zhang (validation)", inhouse: "In-house (validation)", bai: "Bai (validation)" };
const EXPR_UNITS = { discovery: "log2 RPKM", zhang: "normalized log expr", inhouse: "normalized log expr", bai: "normalized log expr" };
const SURV = { // first = default
  discovery: [["Recurrence-free survival", "time_rfs", "status_rfs"], ["Overall survival", "time_os", "status_os"]],
  zhang: [["Progression-free survival", "time_pfs", "status_pfs"], ["Overall survival", "time_os", "status_os"]],
};
const PLOTCFG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["select2d", "lasso2d"] };
const BASELAYOUT = {
  font: { family: "-apple-system, Helvetica, Arial, sans-serif", size: 13 },
  paper_bgcolor: "#fff", plot_bgcolor: "#fff",
  margin: { t: 50, r: 20, b: 55, l: 60 },
};

let db, conn, samples = [], geneSet = null;

// ---------------------------------------------------------------- boot
// Injected post-unlock so the served HTML carries no identifying text
// (see the note in index.html).
function fillHeader() {
  const h = document.getElementById("hdr");
  if (!h) return;
  h.innerHTML = `<h1>Chordoma Multi-omic Atlas</h1>
    <p>Companion resource to the integrative multi-omic classification of chordoma — four molecular
       subgroups (C1–C4) defined in 135 discovery tumors and reproduced across three validation cohorts.
       Query expression, proteomics, methylation, and outcome entirely in your browser.</p>`;
  document.title = "Chordoma Multi-omic Atlas";
  const cite = document.getElementById("citation");
  if (cite) {
    cite.textContent = "[Authors]. An integrative multi-omic classification of chordoma. " +
      "[Journal, Year]. Data: [Zenodo DOI].";
    document.getElementById("citebox").style.display = "";
  }
}

async function boot() {
  fillHeader();
  await loadDataMap();   // must precede any pqUrl() call
  const bundles = duckdb.getJsDelivrBundles();
  const bundle = await duckdb.selectBundle(bundles);
  const workerUrl = URL.createObjectURL(
    new Blob([`importScripts("${bundle.mainWorker}");`], { type: "text/javascript" }));
  db = new duckdb.AsyncDuckDB(new duckdb.ConsoleLogger(duckdb.LogLevel.WARNING),
    new Worker(workerUrl));
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
  conn = await db.connect();
  samples = await q(`SELECT sample_id, cohort, subtype FROM read_parquet('${pqUrl("samples")}')`);
  setStatus("");
  render();
}

async function q(sql) {
  const res = await conn.query(sql);
  return res.toArray().map((r) => Object.fromEntries(r));
}

function setStatus(msg, err = false) {
  const el = document.getElementById("status");
  el.textContent = msg;
  el.className = err ? "status err" : "status";
}

const sub = {}; // sample_id -> subtype per cohort
function subtypeOf(cohort) {
  if (!sub[cohort]) {
    sub[cohort] = new Map(samples.filter((s) => s.cohort === cohort)
      .map((s) => [s.sample_id, s.subtype]));
  }
  return sub[cohort];
}

// wide row -> [{sample, value, subtype}]
function meltRow(row, idCol, cohort) {
  const m = subtypeOf(cohort);
  const out = [];
  for (const [k, v] of Object.entries(row)) {
    if (k === idCol || v == null || Number.isNaN(v)) continue;
    out.push({ sample: k, value: Number(v), subtype: m.get(k) ?? null });
  }
  return out;
}

async function fetchGene(gene, cohort) {
  const table = cohort === "protein" ? "protein_zhang" : `expr_${cohort}`;
  const idCol = cohort === "protein" ? "gene_symbol" : "gene";
  const rows = await q(`SELECT * FROM read_parquet('${pqUrl(table)}')
                        WHERE "${idCol}" = '${gene.replace(/'/g, "''")}'`);
  return rows.length ? meltRow(rows[0], idCol, cohort === "protein" ? "zhang" : cohort) : [];
}

async function geneExists(gene) {
  if (!geneSet) {
    const g = await q(`SELECT DISTINCT gene FROM read_parquet('${pqUrl("gene_index")}')`);
    geneSet = new Set(g.map((r) => r.gene));
  }
  return geneSet.has(gene);
}

// ------------------------------------------------------- gene-set utilities
// Free-text -> canonical symbol list. Accepts any mix of commas, semicolons,
// tabs, newlines, pipes and spaces as separators, tolerates quotes/brackets and
// trailing punctuation, dedupes case-insensitively while keeping input order.
function parseGeneList(raw) {
  const seen = new Set();
  const out = [];
  for (let tok of String(raw ?? "").split(/[\s,;|]+/)) {
    tok = tok.replace(/^["'`([{<]+|["'`)\]}>.,;:]+$/g, "").trim();
    if (!tok) continue;
    const key = tok.toUpperCase();
    if (seen.has(key)) continue;
    seen.add(key);
    out.push(tok);
  }
  return out;
}

let geneUpper = null; // UPPERCASE -> canonical symbol, across all cohorts
async function resolveSymbols(tokens) {
  if (!geneSet) await geneExists("TBXT");
  if (!geneUpper) {
    geneUpper = new Map();
    for (const g of geneSet) if (!geneUpper.has(g.toUpperCase())) geneUpper.set(g.toUpperCase(), g);
  }
  const resolved = [], unknown = [];
  for (const t of tokens) {
    const hit = geneUpper.get(t.toUpperCase());
    (hit ? resolved : unknown).push(hit ?? t);
  }
  return { resolved: [...new Set(resolved)], unknown };
}

// One range-read per cohort for the whole set: [{gene, data:[{sample,value,subtype}]}]
async function fetchGeneSetRows(genes, cohort) {
  if (!genes.length) return [];
  const table = cohort === "protein" ? "protein_zhang" : `expr_${cohort}`;
  const idCol = cohort === "protein" ? "gene_symbol" : "gene";
  const list = genes.map((g) => `'${g.replace(/'/g, "''")}'`).join(",");
  const rows = await q(`SELECT * FROM read_parquet('${pqUrl(table)}')
                        WHERE "${idCol}" IN (${list})`);
  const cohortKey = cohort === "protein" ? "zhang" : cohort;
  return rows.map((r) => ({ gene: String(r[idCol]), data: meltRow(r, idCol, cohortKey) }));
}

// Sum-of-z-scores signature score.
//   1. per gene, z = (x - mean) / sd across the SUBTYPED samples of this cohort
//   2. per sample, score = sum of those z (mode "mean" divides by gene count)
// Genes with zero variance or <3 observations contribute nothing and are
// dropped from `used`; samples covered by <50% of the usable genes are dropped.
// Sparse samples are rescaled (mean z x gene count) so the sum stays on the
// same scale as a fully covered sample.
function scoreGeneSet(geneRows, mode = "sum") {
  const acc = new Map();
  const used = [], dropped = [];
  for (const { gene, data } of geneRows) {
    const d = data.filter((x) => x.subtype && Number.isFinite(x.value));
    const v = d.map((x) => x.value);
    if (d.length < 3 || !(S.sd(v) > 0)) { dropped.push(gene); continue; }
    const m = S.mean(v), s = S.sd(v);
    used.push(gene);
    for (const x of d) {
      const e = acc.get(x.sample) ?? { sum: 0, n: 0, subtype: x.subtype };
      e.sum += (x.value - m) / s;
      e.n++;
      acc.set(x.sample, e);
    }
  }
  const k = used.length;
  const minCov = Math.max(2, Math.ceil(0.5 * k));
  const data = [];
  for (const [sample, e] of acc) {
    if (e.n < minCov) continue;
    const meanZ = e.sum / e.n;
    data.push({ sample, subtype: e.subtype, value: mode === "mean" ? meanZ : meanZ * k, nGenes: e.n });
  }
  return { data, used, dropped };
}

// ---------------------------------------------------------------- plotting
function pText(p) { return p < 1e-4 ? "p < 0.0001" : "p = " + p.toPrecision(3); }

// Survival-tab p-value display: plain decimal down to 0.001, scientific (1.0e-4 style) below.
function pFmt(p) { return p < 1e-3 ? p.toExponential(1) : p.toFixed(3); }
function pFmtEq(p) { return "p = " + pFmt(p); }

function violinTraces(data) {
  return ORDER.filter((s) => data.some((d) => d.subtype === s)).map((s) => {
    const vals = data.filter((d) => d.subtype === s).map((d) => d.value);
    return {
      type: "violin", x: vals.map(() => s), y: vals, name: s,
      fillcolor: COLORS[s], opacity: 0.95, spanmode: "hard", points: false,
      line: { color: "black", width: 1.5 },
      box: { visible: true, width: 0.22, fillcolor: "#DCDCDC", line: { color: "black", width: 1.2 } },
      showlegend: false,
    };
  });
}

function violinLayout(title, ylab, data) {
  const groups = ORDER.map((s) => data.filter((d) => d.subtype === s).map((d) => d.value))
    .filter((g) => g.length >= 3);
  const kw = groups.length >= 2 ? S.kruskal(groups) : null;
  const n = data.filter((d) => ORDER.includes(d.subtype)).length; // total plotted across subtypes
  return {
    ...BASELAYOUT, title: { text: `${title} (N = ${n})`, font: { size: 15 } },
    yaxis: { title: ylab, zeroline: false }, xaxis: { categoryorder: "array", categoryarray: ORDER },
    annotations: kw ? [{ x: 0.02, y: 0.98, xref: "paper", yref: "paper", showarrow: false,
      text: pText(kw.p), font: { size: 14 } }] : [],
  };
}

function kmTraces(groups) { // groups: [{label, color, time, event}]
  const tr = [];
  for (const g of groups) {
    const km = S.kmCurve(g.time, g.event);
    tr.push({
      x: km.points.map((p) => p.t / 365.25), y: km.points.map((p) => p.s),
      mode: "lines", line: { shape: "hv", width: 2.5, color: g.color },
      name: `${g.label} (n=${g.time.length})`, type: "scatter",
      legendgroup: g.label,
    });
    if (km.censors.length) {
      tr.push({ // ggsurvplot-style "+" censor ticks
        x: km.censors.map((p) => p.t / 365.25), y: km.censors.map((p) => p.s),
        mode: "markers", type: "scatter", showlegend: false, legendgroup: g.label,
        marker: { symbol: "cross-thin-open", size: 7, color: g.color,
                  line: { width: 1.5, color: g.color } },
        hoverinfo: "skip",
      });
    }
  }
  return tr;
}

// axis limits matching the manuscript figures (days -> years via 365.25)
const KM_XMAX_YEARS = { discovery: 10, zhang: 12 };

function downloadCSV(rows, filename) {
  if (!rows.length) return;
  const cols = Object.keys(rows[0]);
  const csv = [cols.join(",")].concat(
    rows.map((r) => cols.map((c) => r[c]).join(","))).join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  a.download = filename;
  a.click();
}

// ---------------------------------------------------------------- UI shell
const PAGES = {
  "Gene expression": pageExpression,
  "Gene set score": pageGeneSet,
  "Survival": pageSurvival,
  "Correlation": pageCorrelation,
  "Methylation": pageMethylation,
  "Cell lines": pageCellLines,
};
let activePage = "Gene expression";

function render() {
  const nav = document.getElementById("nav");
  nav.innerHTML = "";
  for (const name of Object.keys(PAGES)) {
    const b = document.createElement("button");
    b.textContent = name;
    b.className = name === activePage ? "active" : "";
    b.onclick = () => { activePage = name; render(); };
    nav.appendChild(b);
  }
  const chips = document.getElementById("chips");
  chips.innerHTML = ORDER.map((s) =>
    `<div class="chip"><div class="swatch" style="background:${COLORS[s]}"></div>
     <div><b>${SUBTYPE_INFO[s][0]}</b><small>${SUBTYPE_INFO[s][1]}</small></div></div>`).join("");
  document.getElementById("pages").innerHTML = "";
  PAGES[activePage](document.getElementById("pages"));
}

function el(html) {
  const d = document.createElement("div");
  d.innerHTML = html;
  return d.firstElementChild;
}

function geneInput(id, value) {
  return `<div class="autocomplete">
    <label for="${id}">Gene</label>
    <input type="text" id="${id}" value="${value}" autocomplete="off" spellcheck="false">
    <div class="ac-list" id="${id}-ac"></div>
  </div>`;
}

async function attachAutocomplete(id, onPick) {
  const input = document.getElementById(id);
  const list = document.getElementById(id + "-ac");
  if (!geneSet) await geneExists("TBXT"); // warm the set
  const genes = [...geneSet].sort();
  input.addEventListener("input", () => {
    const v = input.value.trim().toUpperCase();
    if (v.length < 2) { list.style.display = "none"; return; }
    const hits = genes.filter((g) => g.toUpperCase().startsWith(v)).slice(0, 12);
    list.innerHTML = hits.map((h) => `<div>${h}</div>`).join("");
    list.style.display = hits.length ? "block" : "none";
    [...list.children].forEach((c) => c.onclick = () => {
      input.value = c.textContent; list.style.display = "none"; onPick(c.textContent);
    });
  });
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter") { list.style.display = "none"; onPick(input.value.trim()); }
  });
  input.addEventListener("blur", () => setTimeout(() => list.style.display = "none", 200));
}

// ---------------------------------------------------------------- pages
async function pageExpression(root) {
  root.appendChild(el(`<div class="card">
    <div class="row">${geneInput("g-expr", "LYST")}</div>
    <div id="expr-plots"></div><div class="note" id="expr-note"></div>
    <button class="dl" id="expr-dl" style="display:none">Download CSV</button>
  </div>`));
  const run = async (gene) => {
    gene = gene.toUpperCase() === gene ? gene : gene; // symbols case-sensitive; try as typed
    const holder = document.getElementById("expr-plots");
    holder.innerHTML = "";
    setStatus(`Querying ${gene}…`);
    const allRows = [];
    const panels = [];
    for (const c of Object.keys(COHORTS)) {
      const data = await fetchGene(gene, c);
      if (data.length) panels.push({ kind: "rna", cohort: c, data });
      data.forEach((d) => allRows.push({ gene, cohort: c, modality: "rna", ...d }));
    }
    const prot = await fetchGene(gene, "protein");
    if (prot.length) panels.push({ kind: "protein", cohort: "zhang", data: prot });
    prot.forEach((d) => allRows.push({ gene, cohort: "zhang", modality: "protein", ...d }));
    setStatus("");
    if (!panels.length) { setStatus(`${gene}: not found in any expression matrix.`, true); return; }
    const grid = el(`<div class="grid2"></div>`);
    holder.appendChild(grid);
    for (const p of panels) {
      const dv = el(`<div class="plot"></div>`);
      grid.appendChild(dv);
      const title = p.kind === "protein" ? `${gene} protein — Zhang` : `${gene} RNA — ${COHORTS[p.cohort]}`;
      const unit = p.kind === "protein" ? "corrected log abundance" : EXPR_UNITS[p.cohort];
      Plotly.newPlot(dv, violinTraces(p.data), violinLayout(title, unit, p.data), PLOTCFG);
    }
    document.getElementById("expr-note").textContent =
      "Violin: subtype distribution with inner quartile box; p: Kruskal–Wallis across subtypes. " +
      "Discovery subtypes assigned by multi-omic clustering; validation cohorts classifier-predicted.";
    const dl = document.getElementById("expr-dl");
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(allRows, `${gene}_expression.csv`);
  };
  await attachAutocomplete("g-expr", run);
  run("LYST");
}

// Seed set: the C3 immunogenic T-cell / cytotoxic program.
const DEFAULT_GENE_SET = "CD2, CD3D, CD3E, CD8A, GZMA, GZMK, PRF1, CXCL9, CXCL10, CXCL13, IL7R";

async function pageGeneSet(root) {
  root.appendChild(el(`<div class="card">
    <div class="block">
      <label for="gs-input">Gene set</label>
      <div class="hint">One score per sample from the summed per-gene z-scores.
        Paste symbols separated by commas, spaces, tabs or newlines — all three work.</div>
      <textarea id="gs-input" spellcheck="false" rows="7"
        placeholder="CD2, CD3D, CD3E&#10;CD8A GZMK PRF1&#10;CXCL9">${DEFAULT_GENE_SET}</textarea>
    </div>
    <div class="row">
      <div><label>Score</label><select id="gs-mode">
        <option value="sum">Sum of z-scores</option>
        <option value="mean">Mean of z-scores</option></select></div>
      <div><label>&nbsp;</label><button class="dl" id="gs-go" style="margin:0">Score gene set</button></div>
    </div>
    <div id="gs-summary" class="note"></div>
    <div id="gs-plots"></div>
    <div class="note" id="gs-note"></div>
    <button class="dl" id="gs-dl" style="display:none">Download CSV</button>
  </div>`));

  const input = document.getElementById("gs-input");
  const modeSel = document.getElementById("gs-mode");

  async function run() {
    const mode = modeSel.value;
    const tokens = parseGeneList(input.value);
    const holder = document.getElementById("gs-plots");
    const summary = document.getElementById("gs-summary");
    const note = document.getElementById("gs-note");
    const dl = document.getElementById("gs-dl");
    holder.innerHTML = ""; note.textContent = ""; dl.style.display = "none";
    if (!tokens.length) { summary.innerHTML = ""; setStatus("Enter at least one gene symbol.", true); return; }

    setStatus("Resolving symbols…");
    const { resolved, unknown } = await resolveSymbols(tokens);
    if (!resolved.length) {
      summary.innerHTML = "";
      setStatus(`None of the ${tokens.length} symbol(s) are present in any matrix: ${unknown.slice(0, 12).join(", ")}`, true);
      return;
    }

    setStatus(`Scoring ${resolved.length} gene(s) across cohorts…`);
    const allRows = [];
    const panels = [];
    for (const c of Object.keys(COHORTS)) {
      const geneRows = await fetchGeneSetRows(resolved, c);
      const sc = scoreGeneSet(geneRows, mode);
      if (sc.data.length) panels.push({ kind: "rna", cohort: c, ...sc });
      sc.data.forEach((d) => allRows.push({ cohort: c, modality: "rna", n_genes_used: sc.used.length, ...d }));
    }
    const protRows = await fetchGeneSetRows(resolved, "protein");
    const protSc = scoreGeneSet(protRows, mode);
    if (protSc.data.length) panels.push({ kind: "protein", cohort: "zhang", ...protSc });
    protSc.data.forEach((d) => allRows.push({ cohort: "zhang", modality: "protein", n_genes_used: protSc.used.length, ...d }));
    setStatus("");

    if (!panels.length) { summary.innerHTML = ""; setStatus("No cohort had usable data for this set.", true); return; }

    const unit = mode === "sum" ? "sum of z-scores" : "mean z-score";
    const grid = el(`<div class="grid2"></div>`);
    holder.appendChild(grid);
    for (const p of panels) {
      const dv = el(`<div class="plot"></div>`);
      grid.appendChild(dv);
      const label = p.kind === "protein" ? "protein — Zhang" : `RNA — ${COHORTS[p.cohort]}`;
      const title = `Gene set (${p.used.length} genes) ${label}`;
      Plotly.newPlot(dv, violinTraces(p.data), violinLayout(title, unit, p.data), PLOTCFG);
    }

    // per-cohort mapping table: how many of the requested genes each matrix carries
    const rowsHtml = panels.map((p) => {
      const name = p.kind === "protein" ? "Zhang protein" : COHORTS[p.cohort] + " RNA";
      const miss = resolved.filter((g) => !p.used.includes(g));
      return `<tr><td style="padding:3px 12px 3px 0"><b>${name}</b></td>
        <td style="padding:3px 12px 3px 0">${p.used.length}/${resolved.length} genes</td>
        <td style="padding:3px 12px 3px 0">${p.data.length} samples</td>
        <td style="padding:3px 0;color:var(--muted)">${miss.length ? "not in matrix: " + miss.join(", ") : "complete"}</td></tr>`;
    }).join("");
    summary.innerHTML =
      `<div style="margin-bottom:6px">${resolved.length} unique symbol(s) resolved: ` +
      resolved.map((g) => `<span class="pill">${g}</span>`).join("") + `</div>` +
      (unknown.length ? `<div style="margin-bottom:6px" class="err">Not found in any matrix (${unknown.length}): ` +
        unknown.join(", ") + `</div>` : "") +
      `<table style="font-size:12.5px;border-collapse:collapse">${rowsHtml}</table>`;

    note.textContent =
      `Each gene is z-scored across the subtyped samples of that cohort, then ${mode === "sum" ? "summed" : "averaged"} ` +
      `per sample. Scores are cohort-internal — comparable within a panel, never across panels. ` +
      `Genes absent from a matrix or with zero variance are dropped (see the table); samples covered by ` +
      `fewer than half the usable genes are excluded, and partially covered samples are rescaled to the full-set scale. ` +
      `p: Kruskal–Wallis across subtypes.`;
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(
      allRows.map((r) => ({ cohort: r.cohort, modality: r.modality, sample: r.sample,
        subtype: r.subtype, score: r.value, n_genes_scored: r.nGenes, n_genes_used: r.n_genes_used })),
      `geneset_${mode}_zscore_${resolved.slice(0, 3).join("-")}${resolved.length > 3 ? `_plus${resolved.length - 3}` : ""}.csv`);
  }

  document.getElementById("gs-go").onclick = run;
  modeSel.onchange = run;
  // Ctrl/Cmd+Enter scores without leaving the box; plain Enter stays a newline.
  input.addEventListener("keydown", (e) => {
    if (e.key === "Enter" && (e.metaKey || e.ctrlKey)) { e.preventDefault(); run(); }
  });
  run();
}

async function pageSurvival(root) {
  root.appendChild(el(`<div class="card">
    <div class="row">
      <div><label>Cohort</label><select id="sv-cohort">
        <option value="discovery">Discovery</option><option value="zhang">Zhang (validation)</option></select></div>
      <div><label>Endpoint</label><select id="sv-endpoint"></select></div>
      <div><label>Stratify by</label><select id="sv-mode"></select></div>
      <div id="sv-genewrap" style="display:none">${geneInput("g-surv", "LYST")}</div>
    </div>
    <div class="metrics" id="sv-metrics"></div>
    <div class="plot plot-km" id="sv-plot"></div>
    <div class="note" id="sv-note"></div>
    <button class="dl" id="sv-dl" style="display:none">Download CSV</button>
  </div>`));

  const cohortSel = document.getElementById("sv-cohort");
  const endSel = document.getElementById("sv-endpoint");
  const modeSel = document.getElementById("sv-mode");
  const geneWrap = document.getElementById("sv-genewrap");
  let gene = "LYST";

  function fillModes() {
    // Protein expression stratification only exists for the Zhang cohort.
    const prevMode = modeSel.value;
    const modes = ["Molecular subtype", "Gene expression"];
    if (cohortSel.value === "zhang") modes.push("Protein expression");
    modeSel.innerHTML = modes.map((m) => `<option>${m}</option>`).join("");
    if (modes.includes(prevMode)) modeSel.value = prevMode;
  }

  function syncControls() {
    const perGene = modeSel.value === "Gene expression" || modeSel.value === "Protein expression";
    geneWrap.style.display = perGene ? "" : "none";
  }

  function kmAxes(cohort, tcol) { // ggsurvplot-matched: 2-year ticks, fixed xlim
    // Zhang OS xlim 14 yr (5113.5 d in the manuscript code); else cohort default
    const xmax = cohort === "zhang" && tcol === "time_os" ? 14 : KM_XMAX_YEARS[cohort];
    return {
      xaxis: { title: "Years", range: [0, xmax], dtick: 2, tick0: 0, zeroline: false },
      yaxis: { title: "Survival probability", range: [0, 1.02], zeroline: false },
    };
  }

  function fillEndpoints() {
    endSel.innerHTML = SURV[cohortSel.value]
      .map((e, i) => `<option value="${i}">${e[0]}</option>`).join("");
  }

  async function run() {
    const cohort = cohortSel.value;
    const [label, tcol, scol] = SURV[cohort][Number(endSel.value)];
    setStatus("Querying clinical data…");
    // Zhang OS: recurrent samples inherit the matched primary's time-to-death
    // (time zero = diagnosis, not this resection) -> restrict to primaries,
    // as in the paper's methods.
    const primaryOnly = cohort === "zhang" && tcol === "time_os";
    const clin = (await q(`SELECT * FROM read_parquet('${pqUrl("clinical_" + cohort)}')
                           WHERE ${tcol} IS NOT NULL AND ${scol} IS NOT NULL
                           ${primaryOnly ? "AND tumor_sample_type = 1" : ""}`))
      .map((r) => ({ ...r, [tcol]: Number(r[tcol]), [scol]: Number(r[scol]),
        subtype: subtypeOf(cohort).get(r.sample_id) }))
      .filter((r) => r.subtype);
    setStatus("");
    const primaryNote = primaryOnly
      ? " OS restricted to primary tumors (recurrent samples share the matched primary's death time)."
      : "";
    const metrics = document.getElementById("sv-metrics");
    const note = document.getElementById("sv-note");
    const dl = document.getElementById("sv-dl");
    metrics.innerHTML = ""; note.textContent = ""; dl.style.display = "none";

    if (modeSel.value === "Molecular subtype") {
      const groups = ORDER.filter((s) => clin.some((r) => r.subtype === s)).map((s) => {
        const g = clin.filter((r) => r.subtype === s);
        return { label: s, color: COLORS[s], time: g.map((r) => r[tcol]), event: g.map((r) => r[scol]) };
      });
      // pooled logrank p via pairwise chi2 is wrong; do k-group logrank
      const kg = kGroupLogrank(groups);
      Plotly.newPlot(document.getElementById("sv-plot"), kmTraces(groups), {
        ...BASELAYOUT, ...kmAxes(cohort, tcol),
        title: { text: `${label} by subtype — ${COHORTS[cohort]} (log-rank ${pFmtEq(kg.p)})`, font: { size: 15 } },
        legend: { x: 0.99, xanchor: "right", y: 0.99, yanchor: "top" },
      }, PLOTCFG);
      note.textContent = `${clin.length} samples, ${clin.reduce((a, r) => a + r[scol], 0)} events. ` +
        `Censored observations marked with +.` + primaryNote;
      return;
    }

    // per-gene mode (RNA, or protein when "Protein expression" is selected)
    const molecule = modeSel.value === "Protein expression" ? "protein" : "rna";
    const molLabel = molecule === "protein" ? "protein" : "RNA";
    setStatus(`Querying ${gene} ${molLabel}…`);
    const expr = await fetchGene(gene, molecule === "protein" ? "protein" : cohort);
    setStatus("");
    if (!expr.length) { setStatus(`${gene}: not found in the Zhang ${molLabel} matrix.`, true); return; }
    const em = new Map(expr.map((d) => [d.sample, d.value]));
    const df = clin.filter((r) => em.has(r.sample_id))
      .map((r) => ({ ...r, value: em.get(r.sample_id) }));
    if (df.length < 20) { setStatus(`${gene}: <20 samples with ${molLabel} + follow-up.`, true); return; }
    const time = df.map((r) => r[tcol]), event = df.map((r) => r[scol]), x = df.map((r) => r.value);
    const cox = S.coxEfron(time, event, x);
    const cut = S.optimalCutpoint(time, event, x);
    if (!cox || !cut) { setStatus("Not enough events for a stable model.", true); return; }
    metrics.innerHTML = `
      <div class="metric"><div class="v">${cox.hr.toFixed(2)}</div><div class="k">Cox HR per SD</div></div>
      <div class="metric"><div class="v">${cox.hrLo.toFixed(2)}–${cox.hrHi.toFixed(2)}</div><div class="k">95% CI</div></div>
      <div class="metric"><div class="v">${pFmt(cox.lrP)}</div><div class="k">Likelihood-ratio p</div></div>
      <div class="metric"><div class="v">${cox.n} (${cox.events})</div><div class="k">n (events)</div></div>`;
    const hiIdx = df.filter((r) => r.value > cut.cutpoint);
    const loIdx = df.filter((r) => r.value <= cut.cutpoint);
    const groups = [
      { label: `${gene} ${molLabel} high`, color: "#B2182B", time: hiIdx.map((r) => r[tcol]), event: hiIdx.map((r) => r[scol]) },
      { label: `${gene} ${molLabel} low`, color: "#2166AC", time: loIdx.map((r) => r[tcol]), event: loIdx.map((r) => r[scol]) },
    ];
    Plotly.newPlot(document.getElementById("sv-plot"), kmTraces(groups), {
      ...BASELAYOUT, ...kmAxes(cohort, tcol),
      title: { text: `${gene} ${molLabel} — ${label}, ${COHORTS[cohort]} (optimal cutpoint ${cut.cutpoint.toFixed(2)})`, font: { size: 15 } },
      legend: { x: 0.99, xanchor: "right", y: 0.99, yanchor: "top" },
    }, PLOTCFG);
    note.textContent =
      `KM at the optimal cutpoint (max-selected log-rank, 10th–90th pct scan, ≥10/arm) is for display; ` +
      `report the Cox statistics on continuous ${molLabel} (LR p = ${pFmt(cox.lrP)}, ` +
      `HR/SD = ${cox.hr.toFixed(2)}). Cutpoint log-rank p (${pFmt(cut.logrankP)}) is ` +
      `optimistically biased by cutpoint selection. Censored observations marked with +.` + primaryNote;
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(df.map((r) => ({ sample_id: r.sample_id, subtype: r.subtype,
      [tcol]: r[tcol], [scol]: r[scol], [molecule]: r.value })), `${gene}_${molLabel}_${cohort}_survival.csv`);
  }

  cohortSel.onchange = () => { fillEndpoints(); fillModes(); syncControls(); run(); };
  endSel.onchange = run;
  modeSel.onchange = () => { syncControls(); run(); };
  await attachAutocomplete("g-surv", (g) => { gene = g; run(); });
  fillEndpoints();
  fillModes();
  syncControls();
  run();
}

function kGroupLogrank(groups) { // exact multivariate log-rank (matches lifelines/R survdiff)
  const rows = [];
  groups.forEach((g, gi) => g.time.forEach((t, i) => rows.push([t, g.event[i], gi])));
  rows.sort((a, b) => a[0] - b[0]);
  const k = groups.length;
  const n = new Array(k).fill(0);
  rows.forEach((r) => n[r[2]]++);
  let N = rows.length;
  const Z = new Array(k).fill(0);                       // O - E per group
  const V = Array.from({ length: k }, () => new Array(k).fill(0)); // covariance
  let i = 0;
  while (i < rows.length) {
    const t = rows[i][0];
    let d = 0;
    const dg = new Array(k).fill(0), rem = new Array(k).fill(0);
    let j = i;
    while (j < rows.length && rows[j][0] === t) {
      if (rows[j][1]) { d++; dg[rows[j][2]]++; }
      rem[rows[j][2]]++;
      j++;
    }
    if (d > 0 && N > 1) {
      const f = d * (N - d) / (N - 1);
      for (let a = 0; a < k; a++) {
        Z[a] += dg[a] - d * n[a] / N;
        for (let b = 0; b < k; b++) {
          V[a][b] += f * ((a === b ? n[a] / N : 0) - n[a] * n[b] / (N * N));
        }
      }
    }
    for (let g = 0; g < k; g++) n[g] -= rem[g];
    N -= rem.reduce((a, v) => a + v, 0);
    i = j;
  }
  // stat = Z' V^{-1} Z on the first k-1 components (V is singular on all k)
  const m = k - 1;
  const A = Array.from({ length: m }, (_, a) => V[a].slice(0, m).concat(Z[a])); // augmented
  for (let col = 0; col < m; col++) { // gaussian elimination w/ partial pivot
    let piv = col;
    for (let r = col + 1; r < m; r++) if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
    [A[col], A[piv]] = [A[piv], A[col]];
    if (Math.abs(A[col][col]) < 1e-12) return { stat: 0, p: 1 };
    for (let r = 0; r < m; r++) {
      if (r === col) continue;
      const fct = A[r][col] / A[col][col];
      for (let c = col; c <= m; c++) A[r][c] -= fct * A[col][c];
    }
  }
  let stat = 0;
  for (let a = 0; a < m; a++) stat += Z[a] * (A[a][m] / A[a][a]);
  return { stat, p: S.chi2Sf(stat, m) };
}

async function pageCorrelation(root) {
  root.appendChild(el(`<div class="card">
    <div class="row">
      ${geneInput("g-corr-a", "LYST").replace('for="g-corr-a">Gene', 'for="g-corr-a">Gene A')}
      ${geneInput("g-corr-b", "TBXT").replace('for="g-corr-b">Gene', 'for="g-corr-b">Gene B')}
      <div><label>Cohort</label><select id="corr-cohort">
        ${Object.entries(COHORTS).map(([v, l]) => `<option value="${v}">${l}</option>`).join("")}
      </select></div>
      <div><label>Method</label><select id="corr-method"><option>Spearman</option><option>Pearson</option></select></div>
    </div>
    <div class="plot" id="corr-plot"></div>
    <div class="note" id="corr-note"></div>
    <button class="dl" id="corr-dl" style="display:none">Download CSV</button>
  </div>`));
  let a = "LYST", b = "TBXT";
  async function run() {
    const cohort = document.getElementById("corr-cohort").value;
    const method = document.getElementById("corr-method").value;
    setStatus(`Querying ${a} and ${b}…`);
    const [da, dbb] = await Promise.all([fetchGene(a, cohort), fetchGene(b, cohort)]);
    setStatus("");
    const mb = new Map(dbb.map((d) => [d.sample, d.value]));
    const pairs = da.filter((d) => mb.has(d.sample))
      .map((d) => ({ sample: d.sample, subtype: d.subtype, x: d.value, y: mb.get(d.sample) }));
    if (pairs.length < 5) { setStatus("Fewer than 5 samples with both genes.", true); return; }
    const fn = method === "Spearman" ? S.spearman : S.pearson;
    const res = fn(pairs.map((p) => p.x), pairs.map((p) => p.y));
    const sym = method === "Spearman" ? "ρ" : "r";
    const traces = ORDER.filter((s) => pairs.some((p) => p.subtype === s)).map((s) => ({
      x: pairs.filter((p) => p.subtype === s).map((p) => p.x),
      y: pairs.filter((p) => p.subtype === s).map((p) => p.y),
      text: pairs.filter((p) => p.subtype === s).map((p) => p.sample),
      mode: "markers", type: "scatter", name: s,
      marker: { color: COLORS[s], size: 8, opacity: 0.8, line: { color: "rgba(0,0,0,.4)", width: 0.5 } },
    }));
    Plotly.newPlot(document.getElementById("corr-plot"), traces, {
      ...BASELAYOUT,
      title: { text: `${a} vs ${b} — ${COHORTS[cohort]}: ${method} ${sym} = ${res.r.toFixed(2)}, ${pText(res.p)} (n=${res.n})`, font: { size: 15 } },
      xaxis: { title: `${a} (${EXPR_UNITS[cohort]})` }, yaxis: { title: `${b} (${EXPR_UNITS[cohort]})` },
    }, PLOTCFG);
    // within-subtype table
    const per = ORDER.map((s) => {
      const g = pairs.filter((p) => p.subtype === s);
      if (g.length < 5) return null;
      const r = fn(g.map((p) => p.x), g.map((p) => p.y));
      return `${s}: ${sym} = ${r.r.toFixed(2)} (${pText(r.p)}, n=${g.length})`;
    }).filter(Boolean).join(" · ");
    document.getElementById("corr-note").textContent = "Within-subtype: " + per;
    const dl = document.getElementById("corr-dl");
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(pairs, `${a}_${b}_${cohort}.csv`);
  }
  await attachAutocomplete("g-corr-a", (g) => { a = g; run(); });
  await attachAutocomplete("g-corr-b", (g) => { b = g; run(); });
  document.getElementById("corr-cohort").onchange = run;
  document.getElementById("corr-method").onchange = run;
  run();
}

// Parse "chr1:1,000,000-2,000,000" (also "1:1000000-2000000", whitespace ok).
function parseLocus(s) {
  const m = String(s).trim().replace(/,/g, "")
    .match(/^(?:chr)?([0-9]{1,2}|X|Y)\s*[:\s]\s*([0-9]+)\s*[-–]\s*([0-9]+)$/i);
  if (!m) return null;
  const start = Number(m[2]), end = Number(m[3]);
  if (!(end > start)) return null;
  return { chr: "chr" + m[1].toUpperCase().replace("CHR", ""), start, end };
}

async function pageMethylation(root) {
  root.appendChild(el(`<div class="card">
    <div class="row">
      <div><label>Query by</label><select id="meth-mode">
        <option value="probe">Probe ID</option><option value="range">Genomic range</option></select></div>
      <div id="meth-probewrap"><label>Probe ID</label>
        <input type="text" id="meth-probe" value="cg26928153" spellcheck="false"></div>
      <div id="meth-rangewrap" style="display:none"><label>Region (hg19)</label>
        <input type="text" id="meth-range" value="chr7:5,527,000-5,532,000" spellcheck="false"
               style="min-width:250px"></div>
      <div><button class="dl" id="meth-go" style="margin:0">Query</button></div>
    </div>
    <div class="metrics" id="meth-metrics"></div>
    <div class="plot" id="meth-plot"></div>
    <div class="note" id="meth-note"></div>
    <button class="dl" id="meth-dl" style="display:none">Download CSV</button>
  </div>`));

  const modeSel = document.getElementById("meth-mode");
  const probeWrap = document.getElementById("meth-probewrap");
  const rangeWrap = document.getElementById("meth-rangewrap");

  function syncControls() {
    const byRange = modeSel.value === "range";
    probeWrap.style.display = byRange ? "none" : "";
    rangeWrap.style.display = byRange ? "" : "none";
  }

  // probe -> chromosome, so we know which per-chromosome file to range-read
  async function chrOfProbe(probe) {
    const r = await q(`SELECT chr FROM read_parquet('${pqUrl("meth_probe_index")}')
                       WHERE meth_site = '${probe.replace(/'/g, "''")}'`);
    return r.length ? r[0].chr : null;
  }

  async function runProbe(probe) {
    setStatus(`Locating ${probe}…`);
    const chr = await chrOfProbe(probe);
    if (!chr) { setStatus(`Probe '${probe}' not found (expects cg… identifiers).`, true); return; }
    setStatus(`Reading ${chr}…`);
    const rows = await q(`SELECT * FROM read_parquet('${pqUrl("meth_" + chr)}')
                          WHERE meth_site = '${probe.replace(/'/g, "''")}'`);
    setStatus("");
    if (!rows.length) { setStatus(`Probe '${probe}' not found in ${chr}.`, true); return; }
    const { pos, ...rest } = rows[0];
    const data = meltRow(rest, "meth_site", "discovery");
    Plotly.newPlot(document.getElementById("meth-plot"), violinTraces(data),
      violinLayout(`${probe} methylation by subtype`, "M-value", data), PLOTCFG);
    document.getElementById("meth-metrics").innerHTML = `
      <div class="metric"><div class="v">${chr}:${Number(pos).toLocaleString()}</div>
        <div class="k">Position (hg19)</div></div>`;
    document.getElementById("meth-note").textContent =
      "M-values, normalized and batch-corrected; discovery cohort. " +
      "Violin: subtype distribution with inner quartile box; p: Kruskal–Wallis across subtypes.";
    const dl = document.getElementById("meth-dl");
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(data.map((d) => ({ probe, chr, pos,
      sample: d.sample, subtype: d.subtype, m_value: d.value })), `${probe}_methylation.csv`);
  }

  // Median M-value per sample across all probes in the window. Median (not
  // mean) so a single extreme probe cannot drag a sample's summary.
  async function runRange(locus) {
    const { chr, start, end } = locus;
    setStatus(`Reading ${chr}:${start.toLocaleString()}-${end.toLocaleString()}…`);
    let rows;
    try {
      rows = await q(`SELECT * FROM read_parquet('${pqUrl("meth_" + chr)}')
                      WHERE pos BETWEEN ${start} AND ${end}`);
    } catch (err) {
      setStatus(`No methylation data for ${chr} (autosomes chr1–chr22 only).`, true); return;
    }
    setStatus("");
    if (!rows.length) {
      setStatus(`No probes in ${chr}:${start.toLocaleString()}-${end.toLocaleString()}.`, true); return;
    }
    // rows -> per-sample list of M-values, then median per sample
    const bySample = new Map();
    for (const r of rows) {
      for (const [k, v] of Object.entries(r)) {
        if (k === "meth_site" || k === "pos" || v == null || Number.isNaN(v)) continue;
        if (!bySample.has(k)) bySample.set(k, []);
        bySample.get(k).push(Number(v));
      }
    }
    const sm = subtypeOf("discovery");
    const data = [];
    for (const [sample, vals] of bySample) {
      if (!vals.length) continue;
      data.push({ sample, subtype: sm.get(sample) ?? null, value: S.median(vals), nProbes: vals.length });
    }
    const label = `${chr}:${start.toLocaleString()}-${end.toLocaleString()}`;
    Plotly.newPlot(document.getElementById("meth-plot"), violinTraces(data),
      violinLayout(`Median methylation, ${label}`, "Median M-value", data), PLOTCFG);
    const kb = ((end - start) / 1000).toFixed(1);
    document.getElementById("meth-metrics").innerHTML = `
      <div class="metric"><div class="v">${rows.length}</div><div class="k">CpG probes in region</div></div>
      <div class="metric"><div class="v">${kb} kb</div><div class="k">Region width</div></div>
      <div class="metric"><div class="v">${data.length}</div><div class="k">Samples</div></div>`;
    document.getElementById("meth-note").textContent =
      `Each sample summarized as the median M-value across the ${rows.length} CpG probe(s) in ` +
      `${label} (hg19). Median resists single-probe outliers; with few probes in a window the ` +
      `summary is correspondingly noisy. p: Kruskal–Wallis across subtypes.`;
    const dl = document.getElementById("meth-dl");
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(data.map((d) => ({ region: label, sample: d.sample,
      subtype: d.subtype, median_m_value: d.value, n_probes: d.nProbes })),
      `${chr}_${start}_${end}_median_methylation.csv`);
  }

  async function run() {
    document.getElementById("meth-metrics").innerHTML = "";
    document.getElementById("meth-note").textContent = "";
    document.getElementById("meth-dl").style.display = "none";
    if (modeSel.value === "range") {
      const locus = parseLocus(document.getElementById("meth-range").value);
      if (!locus) { setStatus("Expected a region like chr7:5,527,000-5,532,000.", true); return; }
      await runRange(locus);
    } else {
      await runProbe(document.getElementById("meth-probe").value.trim());
    }
  }

  modeSel.onchange = () => { syncControls(); run(); };
  document.getElementById("meth-go").onclick = run;
  for (const id of ["meth-probe", "meth-range"]) {
    document.getElementById(id).addEventListener("keydown", (e) => { if (e.key === "Enter") run(); });
  }
  syncControls();
  run();
}

// ---------------------------------------------------------------- cell lines
let clMeta = null, clGenes = null;
async function loadCellLineMeta() {
  if (!clMeta) {
    clMeta = await q(`SELECT * FROM read_parquet('${pqUrl("cellline_meta")}')`);
  }
  return clMeta;
}

async function fetchCellLineGene(gene) {
  const rows = await q(`SELECT * FROM read_parquet('${pqUrl("tpm_cellline")}')
                        WHERE gene_name = '${gene.replace(/'/g, "''")}'`);
  if (!rows.length) return [];
  const meta = new Map((await loadCellLineMeta()).map((m) => [m.replicate, m]));
  const out = [];
  for (const [k, v] of Object.entries(rows[0])) {
    if (k === "gene_name" || v == null) continue;
    const m = meta.get(k);
    if (m) out.push({ replicate: k, sample: m.sample_id, type: m.sample_type, tpm: Number(v) });
  }
  return out;
}

const log10p1 = (v) => Math.log10(v + 1);
const median = (a) => {
  const s = [...a].sort((x, y) => x - y);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};

async function pageCellLines(root) {
  root.appendChild(el(`<div class="card">
    <div class="row">
      ${geneInput("g-ccl", "TBXT").replace(">Gene<", ">Gene(s) — comma or space separated<")}
      <div><label>Include PDXs</label><select id="ccl-pdx">
        <option value="no">Cell lines only</option><option value="yes">Cell lines + PDXs</option></select></div>
    </div>
    <div class="plot" id="ccl-plot"></div>
    <div class="note">Chordoma cell line / PDX RNA-seq (TPM). Bars: median of replicates,
      error bars: ±SD, points: individual replicates; y = log10(TPM+1). Multiple genes
      switch to a per-gene dot view grouped by sample. Data from the CCL Tracker resource.</div>
    <button class="dl" id="ccl-dl" style="display:none">Download CSV</button>
  </div>`));

  async function run(raw) {
    const genes = [...new Set(raw.split(/[\s,;]+/).map((g) => g.trim()).filter(Boolean))];
    if (!genes.length) return;
    const includePdx = document.getElementById("ccl-pdx").value === "yes";
    setStatus(`Querying ${genes.length} gene(s)…`);
    const all = (await Promise.all(genes.map(async (g) => (await fetchCellLineGene(g)).map((d) => ({ gene: g, ...d })))))
      .flat().filter((d) => includePdx || d.type === "Cell Line");
    setStatus("");
    const found = [...new Set(all.map((d) => d.gene))];
    if (!found.length) { setStatus(`No cell-line data for: ${genes.join(", ")}`, true); return; }
    const missing = genes.filter((g) => !found.includes(g));
    if (missing.length) setStatus(`Not found: ${missing.join(", ")}`, true);

    const holder = document.getElementById("ccl-plot");
    if (found.length === 1) {
      // single gene: median bar + SD error + replicate points, sorted by median (CCL Tracker style)
      const bySample = {};
      all.forEach((d) => (bySample[d.sample] = bySample[d.sample] || []).push(d.tpm));
      const stats = Object.entries(bySample).map(([sample, v]) => {
        const med = median(v);
        const sdv = v.length > 1 ? Math.sqrt(v.reduce((s, x) => s + (x - S.mean(v)) ** 2, 0) / (v.length - 1)) : 0;
        return { sample, med, lo: log10p1(Math.max(med - sdv, 0)), hi: log10p1(med + sdv), y: log10p1(med) };
      }).sort((a, b) => a.med - b.med);
      const order = stats.map((s) => s.sample);
      const traces = [
        { type: "bar", x: order, y: stats.map((s) => s.y),
          error_y: { type: "data", symmetric: false,
            array: stats.map((s) => s.hi - s.y), arrayminus: stats.map((s) => s.y - s.lo),
            color: "#333", thickness: 1.2, width: 4 },
          marker: { color: "#31688E", opacity: 0.75 }, name: "median", showlegend: false },
        { type: "scatter", mode: "markers",
          x: all.map((d) => d.sample), y: all.map((d) => log10p1(d.tpm)),
          text: all.map((d) => d.replicate),
          marker: { color: "#1a1a2e", size: 6, opacity: 0.65 }, name: "replicate", showlegend: false },
      ];
      Plotly.newPlot(holder, traces, {
        ...BASELAYOUT, title: { text: `${found[0]} — chordoma cell lines${includePdx ? " + PDXs" : ""}`, font: { size: 15 } },
        xaxis: { categoryorder: "array", categoryarray: order, tickangle: -45 },
        yaxis: { title: "log10(TPM + 1)", zeroline: false }, height: 520, margin: { ...BASELAYOUT.margin, b: 110 },
      }, PLOTCFG);
    } else {
      // multi gene: per-sample box of gene means + colored gene dots (CCL Tracker style)
      const perGeneSample = {};
      all.forEach((d) => {
        const k = d.gene + "\u0000" + d.sample;
        (perGeneSample[k] = perGeneSample[k] || { gene: d.gene, sample: d.sample, v: [] }).v.push(log10p1(d.tpm));
      });
      const pts = Object.values(perGeneSample).map((e) => ({ gene: e.gene, sample: e.sample, y: S.mean(e.v) }));
      const sampleMin = {};
      pts.forEach((p) => sampleMin[p.sample] = Math.min(sampleMin[p.sample] ?? Infinity, p.y));
      const order = Object.keys(sampleMin).sort((a, b) => sampleMin[a] - sampleMin[b]);
      const traces = [
        { type: "box", x: pts.map((p) => p.sample), y: pts.map((p) => p.y),
          fillcolor: "rgba(0,0,0,0.06)", line: { color: "#666", width: 1.2 },
          boxpoints: false, showlegend: false },
        ...found.map((g, i) => ({
          type: "scatter", mode: "markers",
          x: pts.filter((p) => p.gene === g).map((p) => p.sample),
          y: pts.filter((p) => p.gene === g).map((p) => p.y),
          name: g, marker: { size: 8, opacity: 0.85 },
        })),
      ];
      Plotly.newPlot(holder, traces, {
        ...BASELAYOUT, title: { text: `${found.length} genes — chordoma cell lines${includePdx ? " + PDXs" : ""}`, font: { size: 15 } },
        xaxis: { categoryorder: "array", categoryarray: order, tickangle: -45 },
        yaxis: { title: "mean log10(TPM + 1)", zeroline: false }, height: 520,
        margin: { ...BASELAYOUT.margin, b: 110 }, boxmode: "overlay",
      }, PLOTCFG);
    }
    const dl = document.getElementById("ccl-dl");
    dl.style.display = "inline-block";
    dl.onclick = () => downloadCSV(all, `celllines_${found.join("_").slice(0, 60)}.csv`);
  }

  // autocomplete against the cell-line gene list (differs from tumor cohorts)
  const input = document.getElementById("g-ccl");
  const list = document.getElementById("g-ccl-ac");
  if (!clGenes) {
    const g = await q(`SELECT gene_name FROM read_parquet('${pqUrl("tpm_cellline")}')`);
    clGenes = g.map((r) => r.gene_name);
  }
  input.addEventListener("input", () => {
    const parts = input.value.split(/[\s,;]+/);
    const last = parts[parts.length - 1].trim().toUpperCase();
    if (last.length < 2) { list.style.display = "none"; return; }
    const hits = clGenes.filter((g) => g.toUpperCase().startsWith(last)).slice(0, 12);
    list.innerHTML = hits.map((h) => `<div>${h}</div>`).join("");
    list.style.display = hits.length ? "block" : "none";
    [...list.children].forEach((c) => c.onclick = () => {
      parts[parts.length - 1] = c.textContent;
      input.value = parts.join(" ");
      list.style.display = "none";
      run(input.value);
    });
  });
  input.addEventListener("keydown", (e) => { if (e.key === "Enter") { list.style.display = "none"; run(input.value); } });
  input.addEventListener("blur", () => setTimeout(() => list.style.display = "none", 200));
  document.getElementById("ccl-pdx").onchange = () => run(input.value);
  run("TBXT");
}

// Wait for gate.js to confirm the passphrase before booting, so no data
// request is issued while the gate is up. gate.js fires this immediately when
// the session is already unlocked.
const startApp = () => boot().catch((e) => {
  // DuckDB-WASM surfaces a missing/unreadable Parquet as an opaque
  // WebAssembly.Exception, which says nothing useful. The overwhelmingly
  // common cause locally is a data file the server cannot find, so name it.
  const opaque = /WebAssembly|Invalid Error/i.test(String(e && e.message));
  setStatus(opaque
    ? `Failed to load data from ${ABS_BASE}/ — check the files exist and are ` +
      `reachable (open the Network tab for 404s). If a datamap.json is present ` +
      `in dev, remove it: local files use readable names.`
    : "Failed to initialize: " + (e && e.message), true);
  console.error(e);
});
// This module is deferred, so gate.js may have already signalled by now;
// check the latch before falling back to the event.
if (globalThis.ATLAS_UNLOCKED) startApp();
else document.addEventListener("atlas-unlock", startApp, { once: true });
