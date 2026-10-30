# pubmed_to_COA — CLAUDE.md

## Project Overview

Python tool that queries PubMed via NCBI Entrez and produces a
**Conflict of Interest / Collaborating Authors (COA)** list for NSF and USDA
grant applications. Two front ends share one library:

| File | Purpose |
|---|---|
| `coa_lib.py` | MEDLINE parsing, row building, TSV/xlsx writers, USDA docx filler, LLM cleanup prompt/schema/merge. Stdlib only at import; openpyxl imported lazily. Must stay runnable under Pyodide. |
| `pubmed_authors_to_COA.py` | CLI: Entrez esearch + batched efetch (100 PMIDs per call), `--format nsf|usda`, `--xlsx`, `--docx`, `--cleanup`. |
| `index.html`, `app.js` | Browser version. Loads Pyodide 0.27.7 from jsdelivr, fetches `coa_lib.py`, calls E-utilities with `fetch`, optional LLM cleanup with a pasted key. Served by GitHub Pages from the repo root of `main`. |
| `usda_COI.docx` | NIFA Conflict of Interest Word template (from nifa.usda.gov). Filled by `fill_usda_docx`; also fetched by the page. |
| `test_wasm.mjs` | Node script that runs coa_lib under Pyodide: parse, TSV, docx, micropip openpyxl, xlsx, apply_cleanup. |
| `.nojekyll` | Keeps GitHub Pages from running Jekyll. |
| `README.md` | User documentation |

### Formats

**nsf** (default): Table 4 of https://nsf-gov-resources.nsf.gov/files/coa_template.xlsx.
Columns `COLLAB(A:), NAME, AFFILIATION, OPTIONAL, LASTACTIVE`. 48-month window.
Affiliation limit 255 chars; dates `mm/dd/yyyy`. Template must not be altered,
so output is pasted into Table 4; the script does not write into the template.

**usda**: NIFA Conflict of Interest List (conflict_of_interest.doc). One table:
`Name | Co-Author | Collaborator | Advisees/Advisors | Other – Specify Nature`,
alphabetical, last name first, `x` under the relationship. 3-year window.
TSV adds `AFFILIATION`, `LASTACTIVE` as extra trailing columns. `--docx` clones
the first empty data row of `usda_COI.docx` for each person (regex on
`word/document.xml`, paraIds stripped) and fills `Name:` from `--pi-name`.
Form must be converted to PDF before submission.

## Key Commands

```bash
# NSF
python pubmed_authors_to_COA.py -q 'Stajich JE' -e you@email.com -o Stajich_COA.tsv --xlsx Stajich_COA.xlsx -v
# USDA
python pubmed_authors_to_COA.py --format usda -q 'Stajich JE' -e you@email.com --pi-name 'Stajich, Jason E' -o Stajich_COI.tsv --docx Stajich_COI.docx
# LLM cleanup (Anthropic via ANTHROPIC_API_KEY / `ant auth login`; OpenAI via OPENAI_API_KEY)
python pubmed_authors_to_COA.py --cleanup-input Stajich_COA.tsv --cleanup auto -o Stajich_COA_clean.tsv
python pubmed_authors_to_COA.py --cleanup-input Stajich_COA.tsv --cleanup openai --cleanup-dry-run -o x.tsv
# Test the Python side under WebAssembly (same runtime as the page)
npm install pyodide@0.27.7 && node test_wasm.mjs
```

On UCR HPCC use `/usr/bin/python3.12` (has BioPython; `pip install --user anthropic openai openpyxl`).

## Parameters

| Flag | Default | Description |
|---|---|---|
| `-e`/`--email` | — | Required by NCBI unless `--cleanup-input` |
| `-q`/`--query` | — | PubMed query; date filter appended as `AND ("<date_from>"[EDAT] : "3000"[EDAT])` |
| `-m`/`--max` | 100 | Max PubMed records |
| `-o`/`--out` | `coi.tsv` | Output TSV |
| `--format` | `nsf` | `nsf` or `usda` (`--usda` is an alias) |
| `--xlsx` | — | Also write xlsx (openpyxl) |
| `--docx` | — | USDA only: fill the Word template |
| `--template` | `usda_COI.docx` next to the script | USDA template path |
| `--pi-name` | — | USDA `Name:` line |
| `--date_from` | 48 months (nsf) / 3 years (usda) ago | `YYYY/MM/DD` |
| `--cleanup` | — | `auto`, `anthropic`, `openai` |
| `--cleanup-model` | `claude-opus-5-5` / `gpt-5` | Model id |
| `--cleanup-input` | — | Read rows from an existing TSV instead of PubMed |
| `--cleanup-dry-run` | off | Write prompt payload to `<out>.cleanup_payload.json`, do not call a model |
| `-v`, `--debug` | off | Verbose / first record only |

## Architecture Notes

- `parse_medline(text)` handles many concatenated records; a `PMID-` line
  starts a record. `FAU` sets the current author (kept even with no `AD`),
  `AD ` lines append affiliations, 6-space-indented lines continue the last
  `AD`, `EDAT` gives the date. `active` is taken from the first record an
  author appears in (records assumed most-recent-first).
- `clean_affiliation` strips emails (moved to OPTIONAL), `Electronic address:`,
  and `RINGGOLD:`/`ROR:`/`GRID:`/`ISNI:` identifiers that newer PubMed records
  append.
- LLM cleanup is provider-neutral in `coa_lib` (`CLEANUP_INSTRUCTIONS`,
  `CLEANUP_SCHEMA`, `cleanup_payload`, `apply_cleanup`). Providers live in the
  CLI (`anthropic` SDK with `output_config` json_schema + beta
  `server-side-fallback-2026-07-01` / `fallbacks="default"`, retried without
  the beta on 400; `openai` chat completions with strict json_schema) and in
  `app.js` (raw fetch; Anthropic needs the
  `anthropic-dangerous-direct-browser-access: true` header). The page sends
  150 rows per request. **Provider calls have not been exercised against a
  live key from this repo** (no Anthropic credential on the HPCC; the
  `OPENAI_API_KEY` in the environment returned 401 on 2026-10-30).
- `apply_cleanup` never drops a row: indices the model omits are appended
  unchanged; newest date wins on merge; `merged_from` records provenance.
- Verified 2026-10-30: NCBI E-utilities return `access-control-allow-origin: *`
  (browser calls work without a proxy). Pyodide 0.27.7 lacks openpyxl in its
  package set, so the page installs it with micropip (pure-Python wheel).

## Known Limitations / TODOs

- Author deduplication across middle-initial variants is not handled
  (e.g., `Smith J` vs `Smith JA` treated as separate people).
- Affiliation is the PubMed AD text with emails and identifiers removed; no
  parsing of department vs. institution vs. city, and street addresses (which
  the NSF template says to omit) are not detected without `--cleanup`.
- The script does not write into the NSF template workbook (openpyxl row
  insertion does not move the template's data validations); rows are pasted
  into Table 4 manually. The USDA template is filled directly.
- Cleanup chunks of 150 rows (page) or one request (CLI) mean a name variant
  pair split across chunks is not merged.
- `LASTACTIVE` ordering assumes Entrez returns records sorted
  most-recent-first; this is the typical behavior but is not enforced.
