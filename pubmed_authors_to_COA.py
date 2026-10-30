#!/usr/bin/env python3
"""
pubmed_authors_to_COA.py - build a co-author Conflict of Interest list from PubMed.

Output formats
  nsf   Table 4 of the NSF COA template (coa_template.xlsx):
        A: | Name | Organizational Affiliation | Optional | Last Active   (48-month window)
  usda  USDA NIFA Conflict of Interest List (conflict_of_interest.doc):
        Name | Co-Author (x) | Collaborator | Advisees/Advisors | Other      (3-year window)

Parsing and formatting live in coa_lib.py, which is shared with the browser
version (index.html). This file adds the PubMed download and the optional
LLM cleanup step.

Author: Jason Stajich (jason.stajich[AT]ucr.edu) / @hyphaltip
"""

import argparse
import json
import os
import sys
from datetime import date, timedelta

from Bio import Entrez

import coa_lib

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_USDA_TEMPLATE = os.path.join(HERE, "usda_COI.docx")
WINDOW_DAYS = {"nsf": 4 * 365, "usda": 3 * 365}
DEFAULT_MODEL = {"anthropic": "claude-opus-5-5", "openai": "gpt-5"}
EFETCH_BATCH = 100


def default_date_from(fmt):
    return (date.today() - timedelta(days=WINDOW_DAYS[fmt])).strftime("%Y/%m/%d")


def parse_args():
    p = argparse.ArgumentParser(
        prog="pubmed_authors_to_COA.py",
        description=(
            "Generate a Conflict of Interest / Collaborating Authors table from a PubMed query.\n"
            "  --format nsf  : NSF COA template Table 4 columns (COLLAB, NAME, AFFILIATION, OPTIONAL, LASTACTIVE),\n"
            "                  48-month window. --xlsx writes a sheet to paste into the template.\n"
            "  --format usda : USDA NIFA COI list (NAME, CO-AUTHOR, COLLABORATOR, ADVISEES/ADVISORS, OTHER),\n"
            "                  3-year window. --docx fills the NIFA Word template directly."
        ),
        epilog=(
            "Examples:\n"
            "  # NSF: TSV plus an xlsx ready for the COA template\n"
            "  pubmed_authors_to_COA.py -q 'Stajich JE' -e you@email.com -o Stajich_COA.tsv --xlsx Stajich_COA.xlsx\n"
            "\n"
            "  # USDA: fill the NIFA Word form (convert to PDF before uploading)\n"
            "  pubmed_authors_to_COA.py --format usda -q 'Stajich JE' -e you@email.com \\\n"
            "      --pi-name 'Stajich, Jason E' -o Stajich_COI.tsv --docx Stajich_COI.docx\n"
            "\n"
            "  # Same, with an LLM pass to merge name variants and shorten affiliations.\n"
            "  # Uses Claude via ANTHROPIC_API_KEY or `ant auth login`, or OpenAI via OPENAI_API_KEY.\n"
            "  pubmed_authors_to_COA.py -q 'Stajich JE' -e you@email.com --cleanup auto -o Stajich_COA.tsv\n"
            "\n"
            "  # Clean up an existing TSV without touching PubMed\n"
            "  pubmed_authors_to_COA.py --cleanup-input Stajich_COA.tsv --cleanup anthropic -o Stajich_COA_clean.tsv\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-e", "--email", help="Email address required by NCBI Entrez")
    p.add_argument("-q", "--query", help="PubMed query; the --date_from filter is appended automatically")
    p.add_argument("-m", "--max", type=int, default=100, help="Max PubMed records (default 100)")
    p.add_argument("-o", "--out", default="coi.tsv", help="Output TSV (default coi.tsv)")
    p.add_argument("--format", choices=coa_lib.FORMATS, default="nsf", help="nsf (default) or usda")
    p.add_argument("--usda", action="store_true", help="Same as --format usda")
    p.add_argument("--xlsx", help="Also write an .xlsx (needs openpyxl)")
    p.add_argument("--docx", help="USDA only: fill the NIFA Word template and write it here")
    p.add_argument("--template", default=DEFAULT_USDA_TEMPLATE,
                   help=f"USDA Word template to fill (default {DEFAULT_USDA_TEMPLATE})")
    p.add_argument("--pi-name", default="", help="USDA only: your name for the 'Name:' line, 'Last, First M'")
    p.add_argument("--date_from", help="Earliest EDAT, YYYY/MM/DD (default: 48 months ago for nsf, 3 years for usda)")
    p.add_argument("--cleanup", choices=["auto", "anthropic", "openai"],
                   help="Run an LLM pass to merge name variants and shorten affiliations")
    p.add_argument("--cleanup-model", help="Model id (default claude-opus-5-5 or gpt-5)")
    p.add_argument("--cleanup-input", help="Skip PubMed; read rows from this TSV (from an earlier run)")
    p.add_argument("--cleanup-dry-run", action="store_true",
                   help="Write the cleanup prompt payload next to --out and exit without calling a model")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--debug", action="store_true", help="Process only the first record")
    a = p.parse_args()
    if a.usda:
        a.format = "usda"
    if a.date_from is None:
        a.date_from = default_date_from(a.format)
    if not a.cleanup_input and not (a.email and a.query):
        p.error("-e/--email and -q/--query are required unless --cleanup-input is given")
    if a.docx and a.format != "usda":
        p.error("--docx is only for --format usda")
    return a


# --------------------------------------------------------------------------
# PubMed
# --------------------------------------------------------------------------
def fetch_rows(args):
    Entrez.email = args.email
    term = f'{args.query} AND ("{args.date_from}"[EDAT] : "3000"[EDAT])'
    ids = Entrez.read(Entrez.esearch(db="pubmed", term=term, retmax=args.max))["IdList"]
    if args.debug:
        ids = ids[:1]
    if args.verbose:
        print(f"{len(ids)} PubMed records for: {term}", file=sys.stderr)
    text = []
    for i in range(0, len(ids), EFETCH_BATCH):
        chunk = ids[i:i + EFETCH_BATCH]
        if args.verbose:
            print(f"fetching PMIDs {chunk[0]} .. {chunk[-1]}", file=sys.stderr)
        text.append(Entrez.efetch(db="pubmed", id=",".join(chunk), rettype="medline", retmode="text").read())
    authors = coa_lib.parse_medline("\n".join(text))
    rows = coa_lib.build_rows(authors)
    if args.verbose:
        for r in rows:
            if r['truncated']:
                print(f"[WARN] affiliation for {r['name']} cut to {coa_lib.NSF_AFFIL_MAXLEN} chars", file=sys.stderr)
    return rows


# --------------------------------------------------------------------------
# LLM cleanup
# --------------------------------------------------------------------------
def pick_provider(choice):
    if choice != "auto":
        return choice
    if os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN") \
            or os.path.isdir(os.path.expanduser("~/.config/anthropic")):
        return "anthropic"
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"
    sys.exit("--cleanup auto: no credentials found. Set ANTHROPIC_API_KEY, run `ant auth login`, "
             "or set OPENAI_API_KEY.")


def cleanup_anthropic(payload, model, verbose):
    import anthropic
    client = anthropic.Anthropic()   # ANTHROPIC_API_KEY, ANTHROPIC_AUTH_TOKEN, or `ant auth login` profile
    kwargs = dict(
        model=model,
        max_tokens=64000,
        system=coa_lib.CLEANUP_INSTRUCTIONS,
        messages=[{"role": "user", "content": payload}],
        output_config={"format": {"type": "json_schema", "schema": coa_lib.CLEANUP_SCHEMA}},
    )
    try:
        # server-side fallback: if the model declines, the API re-runs on another model in the same call
        with client.beta.messages.stream(betas=["server-side-fallback-2026-07-01"],
                                         fallbacks="default", **kwargs) as s:
            msg = s.get_final_message()
    except anthropic.BadRequestError as e:
        if verbose:
            print(f"[cleanup] fallback parameter rejected ({e.message}); retrying without it", file=sys.stderr)
        with client.messages.stream(**kwargs) as s:
            msg = s.get_final_message()
    if msg.stop_reason == "refusal":
        sys.exit("cleanup: the model declined the request")
    if msg.stop_reason == "max_tokens":
        sys.exit("cleanup: output was cut off at max_tokens; try fewer rows per run")
    text = next(b.text for b in msg.content if b.type == "text")
    return json.loads(text)


def cleanup_openai(payload, model, verbose):
    import openai
    client = openai.OpenAI()         # OPENAI_API_KEY
    resp = client.chat.completions.create(
        model=model,
        messages=[{"role": "system", "content": coa_lib.CLEANUP_INSTRUCTIONS},
                  {"role": "user", "content": payload}],
        response_format={"type": "json_schema",
                         "json_schema": {"name": "coa_cleanup", "strict": True,
                                         "schema": coa_lib.CLEANUP_SCHEMA}},
    )
    choice = resp.choices[0]
    if choice.finish_reason == "length":
        sys.exit("cleanup: output was cut off; try fewer rows per run")
    if getattr(choice.message, "refusal", None):
        sys.exit(f"cleanup: the model declined: {choice.message.refusal}")
    return json.loads(choice.message.content)


def run_cleanup(rows, args):
    provider = pick_provider(args.cleanup)
    model = args.cleanup_model or DEFAULT_MODEL[provider]
    payload = coa_lib.cleanup_payload(rows)
    if args.cleanup_dry_run:
        path = args.out + ".cleanup_payload.json"
        with open(path, "w") as fh:
            fh.write(payload)
        print(f"dry run: wrote {len(rows)} rows of prompt payload to {path} (provider {provider}, model {model})")
        sys.exit(0)
    if args.verbose:
        print(f"[cleanup] {provider} {model}: sending {len(rows)} rows", file=sys.stderr)
    result = cleanup_anthropic(payload, model, args.verbose) if provider == "anthropic" \
        else cleanup_openai(payload, model, args.verbose)
    new_rows = coa_lib.apply_cleanup(rows, result)
    merges = [r for r in new_rows if len(r.get('merged_from', [])) > 1]
    log = args.out + ".merges.txt"
    with open(log, "w") as fh:
        for r in merges:
            fh.write(f"{r['name']}\t<= {' | '.join(r['merged_from'])}\n")
    print(f"cleanup: {len(rows)} rows -> {len(new_rows)} people; {len(merges)} merges listed in {log}",
          file=sys.stderr)
    return new_rows


# --------------------------------------------------------------------------
def main():
    args = parse_args()
    if args.cleanup_input:
        with open(args.cleanup_input) as fh:
            rows = coa_lib.rows_from_tsv(fh.read())
        if not args.cleanup:
            args.cleanup = "auto"
    else:
        rows = fetch_rows(args)

    if args.cleanup:
        rows = run_cleanup(rows, args)

    with open(args.out, "w") as fh:
        fh.write(coa_lib.rows_to_tsv(rows, args.format))
    if args.verbose:
        print(f"wrote {len(rows)} rows to {args.out} ({args.format} format)", file=sys.stderr)

    if args.xlsx:
        with open(args.xlsx, "wb") as fh:
            fh.write(coa_lib.rows_to_xlsx_bytes(rows, args.format))
        if args.verbose:
            print(f"wrote {args.xlsx}", file=sys.stderr)

    if args.docx:
        with open(args.template, "rb") as fh:
            template = fh.read()
        with open(args.docx, "wb") as fh:
            fh.write(coa_lib.fill_usda_docx(template, rows, pi_name=args.pi_name))
        print(f"wrote {args.docx}; convert it to PDF before attaching it to the application", file=sys.stderr)


if __name__ == "__main__":
    main()
