#!/usr/bin/env python3
import sys, csv, re
from pathlib import Path
from datetime import datetime
from Bio import SeqIO

alignment, metadata, tips_str, outdir = sys.argv[1:5]
keep = set(tips_str.split(","))
outdir = Path(outdir)
outdir.mkdir(parents=True, exist_ok=True)

records = [r for r in SeqIO.parse(alignment, "fasta") if r.id.split()[0] in keep]
print(f"Extracted {len(records)} sequences: {[r.id for r in records]}")
assert len(records) >= 2, "Need at least JN + one Ethiopian sequence"
SeqIO.write(records, outdir / "JN_ET.fasta", "fasta")

def normalize_date(raw, name):
    s = (raw or "").strip()
    if not s:
        return "2025.1" if name.startswith("ET_") else "2008.5"
    # already decimal year
    if re.fullmatch(r"\d{4}(\.\d+)?", s):
        return s if "." in s else f"{s}.5"
    # YYYY-MM-DD
    m = re.fullmatch(r"(\d{4})-(\d{1,2})-(\d{1,2})", s)
    if m:
        y, mo, d = map(int, m.groups())
        return f"{y + (datetime(y, mo, d) - datetime(y, 1, 1)).days / 365.25:.4f}"
    # DD/MM/YYYY or DD/MM/YY
    m = re.fullmatch(r"(\d{1,2})/(\d{1,2})/(\d{2,4})", s)
    if m:
        d, mo, y = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if y < 100:
            y += 2000
        return f"{y + (datetime(y, mo, d) - datetime(y, 1, 1)).days / 365.25:.4f}"
    return s

dates = {}
with open(metadata, newline="", encoding="utf-8-sig") as fh:
    reader = csv.DictReader(fh)
    fm = {k.lower().strip(): k for k in (reader.fieldnames or [])}
    nk = fm.get("name") or fm.get("sample") or fm.get("strain")
    dk = fm.get("date")
    for r in reader:
        name = (r.get(nk) or "").strip()
        if name in keep:
            dates[name] = normalize_date(r.get(dk, ""), name)

with open(outdir / "dates.csv", "w") as out:
    out.write("strain,date\n")
    for name in sorted(keep):
        d = dates.get(name, "2025.1" if name.startswith("ET_") else "2008.5")
        out.write(f"{name},{d}\n")
print("Dates written:")
print((outdir / "dates.csv").read_text())
