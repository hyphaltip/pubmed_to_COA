"""
coa_lib.py - pure-Python helpers shared by the CLI (pubmed_authors_to_COA.py)
and the browser page (index.html / app.js running under Pyodide).

No network access and no non-stdlib imports at module level, so the same file
runs unchanged in CPython and in Pyodide. openpyxl is imported lazily.

Formats
-------
nsf  : Table 4 of the NSF COA template
       https://nsf-gov-resources.nsf.gov/files/coa_template.xlsx
       columns: A ("A:"), Name, Organizational Affiliation,
                Optional (email, Department), Last Active (mm/dd/yyyy)
usda : USDA NIFA Conflict of Interest list; same data without the "A:" column.
"""

import csv
import io
import json
import re
from datetime import datetime

NSF_AFFIL_MAXLEN = 255          # stated limit in the NSF template
FORMATS = ("nsf", "usda")

_tagline = re.compile(r'^\S+\s*-\s+')
_email_re = re.compile(r'[\w.+-]+@[\w-]+(?:\.[\w-]+)+')


# --------------------------------------------------------------------------
# MEDLINE parsing
# --------------------------------------------------------------------------
def parse_medline(text):
    """Parse one or more MEDLINE flat-file records (as returned by efetch
    rettype=medline). Returns a dict keyed by full author name (FAU):
        {name: {'affils': [str, ...], 'active': 'YYYY/MM/DD', 'pmids': [..]}}
    Records are assumed to arrive most-recent-first, so 'active' is taken from
    the first record in which an author appears.
    """
    authors = {}
    author = None
    in_ad = False
    these = []
    active = ""

    def close_record():
        for a in these:
            if 'active' not in authors[a]:
                authors[a]['active'] = active

    for line in text.splitlines():
        if line.startswith("PMID-"):
            close_record()
            these, active, author, in_ad = [], "", None, False
            pmid = _tagline.sub('', line.strip())
            continue
        if line.startswith("FAU"):
            author = _tagline.sub('', line.strip())
            these.append(author)
            authors.setdefault(author, {'affils': [], 'pmids': []})
            authors[author]['pmids'].append(pmid)
            in_ad = False
        elif line.startswith("AD "):
            in_ad = True
            if author is not None:
                authors[author]['affils'].append(_tagline.sub('', line.strip()))
        elif line.startswith("EDAT"):
            active = _tagline.sub('', line.strip()).split()[0]
            in_ad = False
        elif in_ad and line.startswith("      "):
            authors[author]['affils'][-1] += " " + line.strip()
        else:
            in_ad = False
    close_record()
    return authors


def clean_affiliation(text):
    """Strip emails and organisation identifiers from an AD string.
    Returns (affiliation, [emails])."""
    emails = _email_re.findall(text)
    text = re.sub(r'Electronic address:\s*', '', text, flags=re.I)
    text = _email_re.sub('', text)
    text = re.sub(r'\b(RINGGOLD|GRID|ISNI):\s*\S+', '', text, flags=re.I)
    text = re.sub(r'\bROR:\s*https?://\S+', '', text, flags=re.I)
    text = re.sub(r'\s+', ' ', text).strip(' .;,')
    return text, emails


def build_rows(authors, maxlen=NSF_AFFIL_MAXLEN):
    """Turn parse_medline() output into sorted output rows:
    [{'name', 'affiliation', 'optional', 'active', 'truncated'}]"""
    rows = []
    for name in sorted(authors):
        seen, emails = [], []
        for ad in authors[name]['affils']:
            affil, found = clean_affiliation(ad)
            emails.extend(e for e in found if e not in emails)
            if affil and affil not in seen:
                seen.append(affil)
        affiliation = "; ".join(seen)
        truncated = False
        if len(affiliation) > maxlen:
            kept = []
            for s in seen:
                if len("; ".join(kept + [s])) > maxlen:
                    break
                kept.append(s)
            affiliation = "; ".join(kept) if kept else seen[0][:maxlen]
            truncated = True
        rows.append({'name': name, 'affiliation': affiliation,
                     'optional': "; ".join(emails),
                     'active': authors[name].get('active', ""),
                     'truncated': truncated})
    return rows


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def to_display_date(edat):
    """EDAT 'YYYY/MM/DD' -> 'mm/dd/yyyy' (the format both forms accept)."""
    try:
        return datetime.strptime(edat, "%Y/%m/%d").strftime("%m/%d/%Y")
    except (ValueError, TypeError):
        return edat or ""


def header(fmt):
    if fmt == "usda":
        return ["NAME", "CO-AUTHOR", "COLLABORATOR", "ADVISEES/ADVISORS", "OTHER", "AFFILIATION", "LASTACTIVE"]
    return ["COLLAB", "NAME", "AFFILIATION", "OPTIONAL", "LASTACTIVE"]


def row_values(r, fmt):
    if fmt == "usda":
        # NIFA columns first; affiliation and date are kept as extra columns for your own records
        return [r['name'], "x", "", "", "", r['affiliation'], to_display_date(r['active'])]
    return ["A:", r['name'], r['affiliation'], r['optional'], to_display_date(r['active'])]


def rows_to_tsv(rows, fmt="nsf"):
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(header(fmt))
    for r in rows:
        w.writerow(row_values(r, fmt))
    return buf.getvalue()


def rows_from_tsv(text):
    """Read a TSV written by rows_to_tsv (either format) back into rows."""
    rdr = csv.reader(io.StringIO(text), delimiter="\t")
    hdr = next(rdr)
    idx = {h: i for i, h in enumerate(hdr)}
    rows = []
    for rec in rdr:
        if not rec:
            continue
        active = rec[idx['LASTACTIVE']] if 'LASTACTIVE' in idx else ""
        try:
            active = datetime.strptime(active, "%m/%d/%Y").strftime("%Y/%m/%d")
        except ValueError:
            pass
        rows.append({'name': rec[idx['NAME']],
                     'affiliation': rec[idx['AFFILIATION']] if 'AFFILIATION' in idx else "",
                     'optional': rec[idx['OPTIONAL']] if 'OPTIONAL' in idx else "",
                     'active': active, 'truncated': False})
    return rows


def rows_to_xlsx_bytes(rows, fmt="nsf"):
    """Return an .xlsx (bytes) with the form's columns and real date cells."""
    import openpyxl
    wb = openpyxl.Workbook()
    ws = wb.active
    if fmt == "nsf":
        ws.title = "COA Table 4"
        ws.append(["", "Name:", "Organizational Affiliation",
                   "Optional  (email, Department)", "Last Active"])
        widths = (4.14, 37.86, 38.86, 29.14, 12.43)   # from the NSF template
    else:
        ws.title = "USDA COI"
        ws.append(USDA_COLUMNS)
        widths = (37.86, 12, 12, 16, 22)
    for r in rows:
        if fmt == "usda":
            ws.append(usda_table_rows([r])[0])
            continue
        try:
            active = datetime.strptime(r['active'], "%Y/%m/%d")
        except (ValueError, TypeError):
            active = r['active'] or ""
        vals = ["A:", r['name'], r['affiliation'], r['optional'], active]
        ws.append(vals)
        if isinstance(active, datetime):
            ws.cell(row=ws.max_row, column=len(vals)).number_format = "m/d/yyyy"
    for col, width in zip("ABCDE", widths):
        ws.column_dimensions[col].width = width
    out = io.BytesIO()
    wb.save(out)
    return out.getvalue()


# --------------------------------------------------------------------------
# LLM cleanup step (provider-neutral: prompt + schema + merge logic)
# --------------------------------------------------------------------------
CLEANUP_SCHEMA = {
    "type": "object",
    "properties": {
        "people": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string",
                             "description": "Canonical name as 'Last, First M' (most complete variant)"},
                    "affiliation": {"type": "string",
                                    "description": "Institution plus city/state/country; no department, no street address; <=255 chars"},
                    "optional": {"type": "string",
                                 "description": "Email or department if helpful for disambiguation, else empty"},
                    "source_indices": {"type": "array", "items": {"type": "integer"},
                                       "description": "0-based indices of the input rows merged into this person"}
                },
                "required": ["name", "affiliation", "optional", "source_indices"],
                "additionalProperties": False
            }
        }
    },
    "required": ["people"],
    "additionalProperties": False
}

CLEANUP_INSTRUCTIONS = (
    "You are cleaning a list of co-authors for a US federal grant conflict-of-interest form "
    "(NSF COA Table 4 / USDA NIFA COI). Input is a JSON array of rows with 0-based index, "
    "name ('Last, First M' from PubMed), affiliation (raw PubMed address text), optional (emails), "
    "and last-active date.\n"
    "Tasks:\n"
    "1. Merge rows that are the same person (e.g. 'Smith, John' and 'Smith, John A'; "
    "   accent or hyphen variants). Only merge when the surname and first name agree and the "
    "   affiliations are compatible. When unsure, keep rows separate.\n"
    "2. Choose the most complete name variant as the canonical name.\n"
    "3. Shorten the affiliation to institution plus city/state/country. Remove departments, "
    "   street addresses, postal codes, and duplicates. If a person lists several institutions, "
    "   join them with '; '. Keep under 255 characters.\n"
    "4. Put every input index into exactly one person's source_indices. Do not drop anyone.\n"
    "Return only JSON matching the schema."
)


def cleanup_payload(rows):
    """JSON text of the rows to send to the model."""
    return json.dumps([{"index": i, "name": r['name'], "affiliation": r['affiliation'],
                        "optional": r['optional'], "last_active": to_display_date(r['active'])}
                       for i, r in enumerate(rows)], ensure_ascii=False)


def apply_cleanup(rows, result):
    """Merge model output (dict matching CLEANUP_SCHEMA) onto rows.
    Last-active becomes the newest date among merged rows. Rows the model did
    not mention are appended unchanged, so no one is lost. Returns new rows."""
    out = []
    used = set()
    for p in result.get("people", []):
        idxs = [i for i in p.get("source_indices", []) if 0 <= i < len(rows) and i not in used]
        if not idxs:
            continue
        used.update(idxs)
        actives = [rows[i]['active'] for i in idxs if rows[i]['active']]
        out.append({'name': p.get("name") or rows[idxs[0]]['name'],
                    'affiliation': (p.get("affiliation") or "")[:NSF_AFFIL_MAXLEN],
                    'optional': p.get("optional") or "",
                    'active': max(actives) if actives else "",
                    'truncated': False,
                    'merged_from': [rows[i]['name'] for i in idxs]})
    for i, r in enumerate(rows):
        if i not in used:
            out.append(dict(r, merged_from=[r['name']]))
    out.sort(key=lambda r: r['name'])
    return out


# --------------------------------------------------------------------------
# USDA NIFA Conflict of Interest List (.docx template)
# --------------------------------------------------------------------------
# The NIFA form (conflict_of_interest.doc) is one table with the columns
#   Name | Co-Author | Collaborator | Advisees/Advisors | Other - Specify Nature
# Names are listed alphabetically, last name first, and the relationship is
# marked with an "x". The form asks for the past three years and must be
# converted to PDF before it is attached to the application.
USDA_COLUMNS = ["Name", "Co-Author", "Collaborator", "Advisees/Advisors", "Other - Specify Nature"]
USDA_MIN_ROWS = 23      # empty data rows in the blank template; keep at least this many


def usda_table_rows(rows):
    """Rows as they appear in the NIFA table: name, then an 'x' under Co-Author."""
    return [[r['name'], "x", "", "", ""] for r in rows]


def _xml_escape(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def fill_usda_docx(template_bytes, rows, pi_name=""):
    """Fill the NIFA conflict_of_interest .docx template with `rows` and
    return the new .docx as bytes. Stdlib only (zipfile + re), so it also runs
    under Pyodide. The template's first empty data row is cloned for each
    person; existing empty rows are replaced."""
    import zipfile
    zin = zipfile.ZipFile(io.BytesIO(template_bytes))
    xml = zin.read('word/document.xml').decode('utf-8')
    tbl_m = re.search(r'<w:tbl>.*?</w:tbl>', xml, flags=re.S)
    if not tbl_m:
        raise ValueError("template has no table")
    tbl = tbl_m.group(0)
    trs = [m for m in re.finditer(r'<w:tr[ >].*?</w:tr>', tbl, flags=re.S)]

    def cell_texts(tr):
        return [''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', c))
                for c in re.findall(r'<w:tc>.*?</w:tc>', tr, flags=re.S)]

    hdr_i = next((i for i, m in enumerate(trs)
                  if cell_texts(m.group(0))[:2] == ["Name", "Co-Author"]), None)
    if hdr_i is None:
        raise ValueError("could not find the Name / Co-Author header row")
    data_trs = trs[hdr_i + 1:]
    proto = data_trs[0].group(0)
    # drop paragraph ids so clones do not collide
    proto = re.sub(r'\s+w14:(paraId|textId)="[^"]*"', '', proto)

    def make_row(values):
        cells = re.findall(r'<w:tc>.*?</w:tc>', proto, flags=re.S)
        out = proto
        for cell, val in zip(cells, values):
            if not val:
                continue
            run = ('<w:r><w:rPr><w:rFonts w:ascii="Arial" w:hAnsi="Arial"/><w:sz w:val="18"/></w:rPr>'
                   f'<w:t xml:space="preserve">{_xml_escape(val)}</w:t></w:r>')
            new_cell = cell.replace('</w:p>', run + '</w:p>', 1)
            out = out.replace(cell, new_cell, 1)
        return out

    values = usda_table_rows(rows)
    while len(values) < USDA_MIN_ROWS:
        values.append(["", "", "", "", ""])
    new_rows = ''.join(make_row(v) for v in values)
    new_tbl = tbl[:data_trs[0].start()] + new_rows + tbl[data_trs[-1].end():]

    if pi_name:
        # the "Name:" cell near the top of the table
        name_tr = next((m.group(0) for m in trs if cell_texts(m.group(0)) == ["Name:"]), None)
        if name_tr:
            run = f'<w:r><w:rPr><w:sz w:val="18"/></w:rPr><w:t xml:space="preserve"> {_xml_escape(pi_name)}</w:t></w:r>'
            filled = name_tr.replace('</w:r>', '</w:r>' + run, 1)
            new_tbl = new_tbl.replace(name_tr, filled, 1)

    xml = xml[:tbl_m.start()] + new_tbl + xml[tbl_m.end():]
    out = io.BytesIO()
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = xml.encode('utf-8') if item.filename == 'word/document.xml' else zin.read(item.filename)
            zout.writestr(item, data)
    return out.getvalue()
