# pubmed_to_COA

Build a co-author **Conflict of Interest / Collaborators and Other Affiliations** list from a PubMed query, in the layout the NSF and USDA NIFA forms expect.

Two ways to use it:

* **Browser page (no install):** https://hyphaltip.github.io/pubmed_to_COA/ . The same Python code runs in the browser through [Pyodide](https://pyodide.org) (WebAssembly). Nothing is sent anywhere except NCBI and, only if you turn it on, the AI provider you choose.
* **Command line:** `pubmed_authors_to_COA.py` (Python 3, BioPython).

Both are built on `coa_lib.py`, which holds the MEDLINE parsing and the output writers.

## Output formats

### NSF (`--format nsf`, the default)

Matches **Table 4** of the [NSF COA template](https://nsf-gov-resources.nsf.gov/files/coa_template.xlsx). The default window is 48 months, as Table 4 requires.

| Template column | Output column | Source |
|---|---|---|
| A (`A:` / `C:`) | `COLLAB` | always `A:` (co-author) |
| Name (Last, First MI) | `NAME` | PubMed `FAU` field, already `Lastname, Firstname M` |
| Organizational Affiliation | `AFFILIATION` | PubMed `AD` field(s), joined with `; `, emails and RINGGOLD/ROR identifiers removed, cut to the template's 255-character limit |
| Optional (email, Department) | `OPTIONAL` | email address(es) found in `AD`, if any |
| Last Active | `LASTACTIVE` | `EDAT` of the most recent record for that author, `mm/dd/yyyy` |

`--xlsx FILE` writes the same five columns with real Excel date cells and the template's column widths. NSF says the template layout must not be altered, so insert rows in Table 4 of the official template and paste these rows in. The template also says affiliations must not include street addresses; PubMed addresses are used as-is, so check them (or run the AI cleanup).

### USDA NIFA (`--format usda`)

Matches the [NIFA Conflict of Interest List](https://www.nifa.usda.gov/sites/default/files/resource/conflict_of_interest.doc): one table with the columns **Name | Co-Author | Collaborator | Advisees/Advisors | Other – Specify Nature**, names alphabetical with last name first, and an `x` under the relationship. The form asks for the past **three years**, which is the default window for this format.

* The TSV has those five columns, then `AFFILIATION` and `LASTACTIVE` as extra columns for your own records (delete them before pasting).
* `--docx FILE` fills the NIFA Word template (`usda_COI.docx`, bundled) directly: one row per co-author with `x` under Co-Author, and your name on the `Name:` line if `--pi-name` is given. NIFA requires the form to be converted to PDF before it is attached.
* `--xlsx FILE` writes the five NIFA columns.

Name variants (with and without a middle initial) are **not** merged automatically. Use the AI cleanup step or edit the list by hand. Authors with no affiliation in PubMed are kept with an empty affiliation.

## Command line

```
pip install biopython openpyxl        # openpyxl only for --xlsx
```

```bash
# NSF: TSV plus an xlsx to paste into the COA template
./pubmed_authors_to_COA.py -q 'Stajich JE' -e you@email.com -o Stajich_COA.tsv --xlsx Stajich_COA.xlsx

# USDA: fill the NIFA Word form
./pubmed_authors_to_COA.py --format usda -q 'Stajich JE' -e you@email.com \
    --pi-name 'Stajich, Jason E' -o Stajich_COI.tsv --docx Stajich_COI.docx

# Extend the window with --date_from (YYYY/MM/DD)
./pubmed_authors_to_COA.py -q '(Jarvis Erich OR Jarvis ED)' --date_from 2020/01/01 -e you@email.com

# Quick test on the first record only
./pubmed_authors_to_COA.py --debug -q 'Treseder KK' -e you@email.com
```

Records are fetched in batches of 100 PMIDs. `-m/--max` caps the number of records (default 100).

### Optional AI cleanup (`--cleanup`)

An LLM pass that merges name variants that are the same person, shortens each affiliation to institution plus city/state/country, and drops departments and street addresses. Every input row ends up in exactly one output person; rows the model does not mention are kept unchanged. Merges are written to `<out>.merges.txt` so you can review them. The newest `LASTACTIVE` among merged rows is kept.

Credentials are read from the environment; nothing is stored by the script.

| `--cleanup` | Credentials used | Default model |
|---|---|---|
| `anthropic` | `ANTHROPIC_API_KEY`, `ANTHROPIC_AUTH_TOKEN`, or the profile from `ant auth login` (Claude login) | `claude-opus-5-5` |
| `openai` | `OPENAI_API_KEY` | `gpt-5` |
| `auto` | Anthropic if any of its credentials are present, else OpenAI | |

```bash
pip install anthropic            # and/or: pip install openai
./pubmed_authors_to_COA.py -q 'Stajich JE' -e you@email.com --cleanup auto -o Stajich_COA.tsv

# Clean an existing TSV without touching PubMed
./pubmed_authors_to_COA.py --cleanup-input Stajich_COA.tsv --cleanup anthropic -o Stajich_COA_clean.tsv

# See exactly what would be sent, without calling a model
./pubmed_authors_to_COA.py --cleanup-input Stajich_COA.tsv --cleanup openai --cleanup-dry-run -o x.tsv
```

Override the model with `--cleanup-model`. The Anthropic call uses structured JSON output and the server-side fallback option, so a declined request is retried on another model inside the same call. **The provider calls were written against the SDK signatures but have not been run against a live key in this repository.**

## Browser page

`index.html` and `app.js` load Pyodide from jsdelivr, fetch `coa_lib.py` and run it unchanged. The page:

1. queries E-utilities directly (NCBI sends `Access-Control-Allow-Origin: *`, so no proxy is needed),
2. renders the table, and offers TSV, `.xlsx` (openpyxl installed into the browser Python on first use) and, for USDA, the filled Word form,
3. optionally runs the same cleanup through the Anthropic or OpenAI API with a key you paste.

Safe use of a key in a browser page:

* The page is static with no backend. The key goes only to the provider's API, straight from the tab. The page's `Content-Security-Policy` lists every host it may contact.
* The key is kept in memory only, unless you tick "remember", which stores it in `localStorage` for this site. "Forget key" clears both.
* Use a dedicated key with a spending limit and revoke it afterwards. Do not use the remember option on shared machines.
* Browser calls to the Anthropic API need the `anthropic-dangerous-direct-browser-access: true` header, which the page sends. The name is deliberate: it is acceptable for a personal tool where you paste your own key, not for a page that would embed a key for other people.

To host your own copy, enable GitHub Pages on the repository root of the `main` branch; no build step is needed. For local development serve the directory with any static server (`python3 -m http.server`); opening `index.html` as a `file://` URL does not work because the page fetches `coa_lib.py`.

## Files

| File | Purpose |
|---|---|
| `coa_lib.py` | Parsing and output writers shared by the CLI and the page (stdlib only; openpyxl loaded lazily) |
| `pubmed_authors_to_COA.py` | Command-line tool: PubMed download, formats, optional LLM cleanup |
| `index.html`, `app.js` | Browser version (Pyodide) |
| `usda_COI.docx` | NIFA Conflict of Interest template, filled by `--docx` and by the page |

## Author

[Jason Stajich](http://lab.stajich.org) [@hyphaltip](https://github.com/hyphaltip)
