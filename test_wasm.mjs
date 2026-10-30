// Headless test of coa_lib.py under Pyodide (WebAssembly), the same runtime the browser page uses.
// Setup once:  npm install pyodide@0.27.7     Run:  node test_wasm.mjs
// Needs network access (NCBI E-utilities, jsdelivr for micropip, PyPI for openpyxl).
import { loadPyodide } from "pyodide";
import fs from "node:fs";
const REPO = new URL(".", import.meta.url).pathname.replace(/\/$/, "");
const py = await loadPyodide();
py.FS.writeFile("coa_lib.py", fs.readFileSync(REPO + "/coa_lib.py", "utf8"));
const lib = py.pyimport("coa_lib");
const medline = await (await fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=42845469,42803773,42788761&rettype=medline&retmode=text")).text();
const rows = lib.build_rows(lib.parse_medline(medline)).toJs({ dict_converter: Object.fromEntries });
console.log("WASM parse:", rows.length, "rows; first:", rows[0].name, "|", rows[0].active);
console.log("TSV head:", lib.rows_to_tsv(py.toPy(rows), "usda").split("\n")[1]);
const tpl = new Uint8Array(fs.readFileSync(REPO + "/usda_COI.docx"));
const docx = lib.fill_usda_docx(py.toPy(tpl), py.toPy(rows), "Stajich, Jason E").toJs();
console.log("WASM docx bytes:", docx.length);
await py.loadPackage("micropip");
const micropip = py.pyimport("micropip");
await micropip.install("openpyxl");
const xlsx = lib.rows_to_xlsx_bytes(py.toPy(rows), "nsf").toJs();
console.log("WASM xlsx bytes:", xlsx.length);
const res = { people: [{ name: rows[0].name, affiliation: "Univ X", optional: "", source_indices: [0] }] };
const cleaned = lib.apply_cleanup(py.toPy(rows), py.toPy(res)).toJs({ dict_converter: Object.fromEntries });
console.log("WASM apply_cleanup:", cleaned.length, "rows; first affil:", cleaned[0].affiliation);
