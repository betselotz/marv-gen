#!/usr/bin/env python3
import sys, csv, re
from datetime import datetime
from Bio import SeqIO

MONTHS = {m:i for i,m in enumerate("jan feb mar apr may jun jul aug sep oct nov dec".split(),1)}
MONTHS.update({m:i for i,m in enumerate(
    "january february march april may june july august september october november december".split(),1)})

def parse_date(s):
    s = str(s).strip().replace("\xa0","")
    if not s or s.lower() in ("nan","none","?"): return None
    if re.fullmatch(r"\d{4}", s): return datetime(int(s),7,1)
    m = re.fullmatch(r"([A-Za-z]+)-(\d{4})", s)
    if m and m.group(1).lower() in MONTHS:
        return datetime(int(m.group(2)), MONTHS[m.group(1).lower()], 15)
    for fmt in ("%Y-%m-%d","%d/%m/%Y","%d/%m/%y","%Y/%m/%d","%d-%b-%Y"):
        try: return datetime.strptime(s, fmt)
        except ValueError: pass
    return None

def to_decimal(dt):
    start = datetime(dt.year,1,1)
    return dt.year + (dt-start).days / (366.0 if dt.year%4==0 else 365.0)

aln, meta, out_fa, out_tsv, out_csv = sys.argv[1:6]
rows = {}
with open(meta, newline="", encoding="utf-8-sig") as fh:
    r = csv.DictReader(fh)
    fm = {k.lower().strip(): k for k in (r.fieldnames or [])}
    nk = fm.get("name") or fm.get("sample") or fm.get("strain")
    dk = fm.get("date")
    for row in r:
        name = (row.get(nk) or "").strip().replace("\xa0","")
        dt = parse_date(row.get(dk) or "")
        if name and dt:
            rows[name] = dt
            rows[name.split(".")[0]] = dt

kept, dates, clean = [], [], []
for rec in SeqIO.parse(aln, "fasta"):
    rid = rec.id.split()[0]
    dt = rows.get(rid) or rows.get(rid.split(".")[0])
    if not dt: continue
    dec = to_decimal(dt)
    kept.append(rec)
    dates.append((rid, dec))
    clean.append({"name": rid, "date": dt.strftime("%Y-%m-%d"), "decimal": f"{dec:.6f}"})

SeqIO.write(kept, out_fa, "fasta")
open(out_tsv,"w").write("name\tdecimal_year\n"+"".join(f"{n}\t{d:.6f}\n" for n,d in dates))
with open(out_csv,"w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name","date","decimal"]); w.writeheader(); w.writerows(clean)
print(f"Prepared {len(kept)} tips")
