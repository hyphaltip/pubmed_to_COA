/* PubMed to COA - browser version.
 * Python (coa_lib.py) runs under Pyodide (WebAssembly). This file does the
 * network calls (NCBI E-utilities, optional LLM provider) and the UI.
 */
"use strict";

const EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
const EFETCH_BATCH = 100;
const CLEANUP_CHUNK = 150;          // rows per LLM request
const WINDOW_DAYS = { nsf: 4 * 365, usda: 3 * 365 };
const DEFAULT_MODEL = { anthropic: "claude-opus-5-5", openai: "gpt-5" };

const $ = (id) => document.getElementById(id);
let pyodide = null, lib = null, rows = null, fmt = "nsf", openpyxlReady = false, usdaTemplate = null;

function log(msg, cls) {
  const el = $("status");
  const line = document.createElement("div");
  if (cls) line.className = cls;
  line.textContent = msg;
  el.appendChild(line);
  el.scrollTop = el.scrollHeight;
}
function clearLog() { $("status").textContent = ""; }

function defaultDateFrom(f) {
  const d = new Date(Date.now() - WINDOW_DAYS[f] * 86400e3);
  return `${d.getFullYear()}/${String(d.getMonth() + 1).padStart(2, "0")}/${String(d.getDate()).padStart(2, "0")}`;
}

// ---------------------------------------------------------------- Pyodide
async function initPython() {
  pyodide = await loadPyodide();
  const src = await (await fetch("coa_lib.py", { cache: "no-cache" })).text();
  pyodide.FS.writeFile("coa_lib.py", src);
  lib = pyodide.pyimport("coa_lib");
  $("pyostate").textContent = "Python runtime ready.";
  $("run").disabled = false;
}

async function ensureOpenpyxl() {
  if (openpyxlReady) return;
  log("installing openpyxl into the browser Python (one time)…");
  await pyodide.loadPackage("micropip");
  const micropip = pyodide.pyimport("micropip");
  await micropip.install("openpyxl");
  openpyxlReady = true;
}

function pyRows(jsRows) { return pyodide.toPy(jsRows); }
function jsRows(pyList) { return pyList.toJs({ dict_converter: Object.fromEntries }); }

// ---------------------------------------------------------------- PubMed
async function fetchPubmed(query, email, dateFrom, max) {
  const term = `${query} AND ("${dateFrom}"[EDAT] : "3000"[EDAT])`;
  const u = new URL(EUTILS + "esearch.fcgi");
  u.search = new URLSearchParams({ db: "pubmed", term, retmax: String(max), retmode: "json", email, tool: "pubmed_to_COA" });
  const res = await fetch(u);
  if (!res.ok) throw new Error(`esearch HTTP ${res.status}`);
  const ids = (await res.json()).esearchresult.idlist;
  log(`${ids.length} PubMed records match: ${term}`);
  let text = "";
  for (let i = 0; i < ids.length; i += EFETCH_BATCH) {
    const chunk = ids.slice(i, i + EFETCH_BATCH);
    log(`fetching PMIDs ${chunk[0]} … ${chunk[chunk.length - 1]}`);
    const body = new URLSearchParams({ db: "pubmed", id: chunk.join(","), rettype: "medline", retmode: "text", email, tool: "pubmed_to_COA" });
    const r = await fetch(EUTILS + "efetch.fcgi", { method: "POST", body });
    if (!r.ok) throw new Error(`efetch HTTP ${r.status}`);
    text += "\n" + await r.text();
    if (i + EFETCH_BATCH < ids.length) await new Promise(ok => setTimeout(ok, 400)); // stay under 3 req/s
  }
  return text;
}

// ---------------------------------------------------------------- render
function render() {
  const hdr = lib.header(fmt).toJs();
  const t = $("table");
  t.innerHTML = "";
  const thead = t.createTHead().insertRow();
  hdr.forEach(h => { const th = document.createElement("th"); th.textContent = h; thead.appendChild(th); });
  const tb = t.createTBody();
  rows.forEach(r => {
    const vals = lib.row_values(pyodide.toPy(r), fmt).toJs();
    const tr = tb.insertRow();
    vals.forEach(v => { tr.insertCell().textContent = v; });
    if (r.truncated) tr.cells[fmt === "nsf" ? 2 : 5].classList.add("warn");
  });
  $("dl-tsv").disabled = $("dl-xlsx").disabled = false;
  $("dl-docx").disabled = fmt !== "usda";
  $("cleanup").disabled = false;
  $("pastehint").textContent = fmt === "nsf"
    ? "Paste the rows into Table 4 of the official NSF template (insert rows first). The template itself must not be altered, and affiliations should not contain street addresses."
    : "Download the filled NIFA Word form, check it, then convert it to PDF before attaching it to the application.";
  const merges = rows.filter(r => (r.merged_from || []).length > 1);
  $("mergesbox").hidden = merges.length === 0;
  $("merges").textContent = merges.map(r => `${r.name}  <=  ${r.merged_from.join(" | ")}`).join("\n");
}

function download(name, data, type) {
  const blob = data instanceof Blob ? data : new Blob([data], { type });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = name; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}
function baseName() { return ($("query").value.replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "") || "coi") + "_" + fmt.toUpperCase(); }

// ---------------------------------------------------------------- LLM cleanup
async function callAnthropic(key, model, system, payload) {
  const body = {
    model, max_tokens: 32000, system,
    messages: [{ role: "user", content: payload }],
    output_config: { format: { type: "json_schema", schema: lib.CLEANUP_SCHEMA.toJs({ dict_converter: Object.fromEntries }) } },
  };
  const headers = {
    "content-type": "application/json",
    "x-api-key": key,
    "anthropic-version": "2023-06-01",
    "anthropic-dangerous-direct-browser-access": "true",   // required for calls from a browser page
  };
  // first try with server-side fallback (re-runs on another model if this one declines)
  let res = await fetch("https://api.anthropic.com/v1/messages", {
    method: "POST", headers: { ...headers, "anthropic-beta": "server-side-fallback-2026-07-01" },
    body: JSON.stringify({ ...body, fallbacks: "default" }),
  });
  if (res.status === 400) {
    res = await fetch("https://api.anthropic.com/v1/messages", { method: "POST", headers, body: JSON.stringify(body) });
  }
  const data = await res.json();
  if (!res.ok) throw new Error(`Anthropic ${res.status}: ${data.error?.message || JSON.stringify(data)}`);
  if (data.stop_reason === "refusal") throw new Error("the model declined the request");
  if (data.stop_reason === "max_tokens") throw new Error("output cut off at max_tokens; try fewer rows");
  return JSON.parse(data.content.find(b => b.type === "text").text);
}

async function callOpenAI(key, model, system, payload) {
  const res = await fetch("https://api.openai.com/v1/chat/completions", {
    method: "POST",
    headers: { "content-type": "application/json", authorization: `Bearer ${key}` },
    body: JSON.stringify({
      model,
      messages: [{ role: "system", content: system }, { role: "user", content: payload }],
      response_format: { type: "json_schema", json_schema: { name: "coa_cleanup", strict: true, schema: lib.CLEANUP_SCHEMA.toJs({ dict_converter: Object.fromEntries }) } },
    }),
  });
  const data = await res.json();
  if (!res.ok) throw new Error(`OpenAI ${res.status}: ${data.error?.message || JSON.stringify(data)}`);
  const ch = data.choices[0];
  if (ch.finish_reason === "length") throw new Error("output cut off; try fewer rows");
  if (ch.message.refusal) throw new Error(`the model declined: ${ch.message.refusal}`);
  return JSON.parse(ch.message.content);
}

async function runCleanup() {
  const provider = $("provider").value, key = $("apikey").value.trim(), model = $("model").value.trim() || DEFAULT_MODEL[provider];
  if (!key) { log("enter an API key first", "warn"); return; }
  if ($("remember").checked) { try { localStorage.setItem("coa_key_" + provider, key); } catch (e) { /* storage blocked */ } }
  $("cleanup").disabled = true;
  try {
    const system = lib.CLEANUP_INSTRUCTIONS;
    let out = [];
    for (let i = 0; i < rows.length; i += CLEANUP_CHUNK) {
      const chunk = rows.slice(i, i + CLEANUP_CHUNK);
      log(`AI cleanup (${provider} ${model}): rows ${i + 1}–${i + chunk.length} of ${rows.length}…`);
      const payload = lib.cleanup_payload(pyRows(chunk));
      const result = provider === "anthropic" ? await callAnthropic(key, model, system, payload) : await callOpenAI(key, model, system, payload);
      out = out.concat(jsRows(lib.apply_cleanup(pyRows(chunk), pyodide.toPy(result))));
    }
    const before = rows.length;
    rows = out.sort((a, b) => a.name.localeCompare(b.name));
    render();
    log(`cleanup done: ${before} rows -> ${rows.length} people. Review the merges list.`, "ok");
  } catch (e) {
    log("cleanup failed: " + e.message, "warn");
  } finally { $("cleanup").disabled = false; }
}

// ---------------------------------------------------------------- wiring
$("format").addEventListener("change", () => { $("datefrom").value = defaultDateFrom($("format").value); });
$("provider").addEventListener("change", () => {
  $("model").value = DEFAULT_MODEL[$("provider").value];
  try { $("apikey").value = localStorage.getItem("coa_key_" + $("provider").value) || ""; } catch (e) { /* ignore */ }
});
$("forget").addEventListener("click", () => {
  $("apikey").value = ""; $("remember").checked = false;
  try { localStorage.removeItem("coa_key_anthropic"); localStorage.removeItem("coa_key_openai"); } catch (e) { /* ignore */ }
  log("key removed from this page and from localStorage");
});

$("run").addEventListener("click", async () => {
  const query = $("query").value.trim(), email = $("email").value.trim();
  if (!query || !email) { log("query and email are both required (NCBI asks for an email)", "warn"); return; }
  fmt = $("format").value;
  $("run").disabled = true; clearLog();
  try {
    const text = await fetchPubmed(query, email, $("datefrom").value.trim(), Number($("max").value) || 200);
    const authors = lib.parse_medline(text);
    rows = jsRows(lib.build_rows(authors));
    render();
    const cut = rows.filter(r => r.truncated).length;
    log(`${rows.length} co-authors. ${cut ? cut + " affiliations were cut to 255 characters (highlighted)." : ""}`, "ok");
  } catch (e) {
    log("failed: " + e.message, "warn");
  } finally { $("run").disabled = false; }
});

$("dl-tsv").addEventListener("click", () => download(baseName() + ".tsv", lib.rows_to_tsv(pyRows(rows), fmt), "text/tab-separated-values"));
$("dl-xlsx").addEventListener("click", async () => {
  try {
    await ensureOpenpyxl();
    const bytes = lib.rows_to_xlsx_bytes(pyRows(rows), fmt).toJs();
    download(baseName() + ".xlsx", new Blob([bytes], { type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet" }));
  } catch (e) { log("xlsx failed: " + e.message, "warn"); }
});
$("dl-docx").addEventListener("click", async () => {
  try {
    if (!usdaTemplate) usdaTemplate = new Uint8Array(await (await fetch("usda_COI.docx")).arrayBuffer());
    const bytes = lib.fill_usda_docx(pyodide.toPy(usdaTemplate), pyRows(rows), $("piname").value.trim()).toJs();
    download(baseName() + ".docx", new Blob([bytes], { type: "application/vnd.openxmlformats-officedocument.wordprocessingml.document" }));
    log("USDA form written. Convert it to PDF before attaching it to the application.", "ok");
  } catch (e) { log("docx failed: " + e.message, "warn"); }
});

$("datefrom").value = defaultDateFrom("nsf");
$("run").disabled = true;
try { $("apikey").value = localStorage.getItem("coa_key_anthropic") || ""; } catch (e) { /* ignore */ }
initPython().catch(e => { $("pyostate").textContent = "Python runtime failed to load: " + e.message; });
