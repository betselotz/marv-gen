process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/13_multiqc", mode: params.publish_dir_mode
    publishDir "${params.outdir}/13_multiqc/summaries", mode: params.publish_dir_mode,
        pattern: "summary_*.csv"

    input:
    path results_dir

    output:
    path "**/*.html", optional: true, emit: reports
    path "summary_*.csv", optional: true, emit: summaries
    path "versions.yml", emit: versions

    script:
    """
    set +e
    RD=${results_dir}

    mkdir -p 01_fastp 02_host 03_taxonomy 05_qualimap 06_variants \\
             08_coverage 09_nextclade 10_lineage 11_phylogeny 12_timetree 13_pipeline

    # ========== 1) Aggregate CSVs (easy for MultiQC + Excel) ==========
    python3 - << 'PY'
from pathlib import Path
import csv, re

rd = Path("${results_dir}")

def write_csv(path, rows, fieldnames):
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("Wrote", path, "rows=", len(rows))

# --- Kraken2: one row per sample (top taxa + unclassified) ---
kraken_rows = []
for rep in sorted((rd / "03_taxonomy").glob("**/*")) if (rd / "03_taxonomy").exists() else []:
    if rep.suffix not in (".txt", ".report", "") and "report" not in rep.name.lower():
        continue
    if not rep.is_file():
        continue
    sample = rep.stem.replace(".kraken2", "").replace("_report", "").replace(".report", "")
    unclass = root = virus = marburg = ""
    top = []
    try:
        for line in rep.read_text(errors="replace").splitlines():
            p = line.split("\\t")
            if len(p) < 6:
                continue
            pct, reads, rank, taxid, name = p[0].strip(), p[1].strip(), p[3].strip(), p[4].strip(), p[5].strip()
            if "unclassified" in name.lower():
                unclass = pct
            if rank == "R" or name.lower() == "root":
                root = pct
            if "viruses" in name.lower() and rank in ("D", "R1", "R"):
                virus = pct
            if "marburg" in name.lower():
                marburg = pct
            if rank in ("S", "S1", "G", "F"):
                top.append((float(pct), name, reads))
        top = sorted(top, reverse=True)[:5]
        kraken_rows.append({
            "sample": sample,
            "pct_unclassified": unclass,
            "pct_viruses": virus,
            "pct_marburg": marburg,
            "top1_taxon": top[0][1] if top else "",
            "top1_pct": top[0][0] if top else "",
            "top2_taxon": top[1][1] if len(top)>1 else "",
            "top2_pct": top[1][0] if len(top)>1 else "",
            "report_file": str(rep.name),
        })
    except Exception as e:
        print("kraken parse skip", rep, e)
if kraken_rows:
    write_csv("summary_kraken2.csv", kraken_rows,
              ["sample","pct_unclassified","pct_viruses","pct_marburg","top1_taxon","top1_pct","top2_taxon","top2_pct","report_file"])

# --- Nextclade: merge all sample TSVs ---
nc_rows = []
nc_dir = rd / "09_nextclade"
if nc_dir.exists():
    for f in sorted(nc_dir.glob("*.tsv")):
        with f.open(newline="", encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\\t")
            for row in r:
                row["_source_file"] = f.name
                nc_rows.append(row)
if nc_rows:
    # keep useful columns if present
    prefer = ["seqName","clade","lineage","qc.overallStatus","qc.overallScore",
              "totalSubstitutions","totalMissing","coverage","totalNonACGTNs"]
    fields = [c for c in prefer if c in nc_rows[0]] + [c for c in nc_rows[0] if c not in prefer]
    write_csv("summary_nextclade.csv", nc_rows, fields)

# --- Variants (calling) ---
var_rows = []
vdir = rd / "06_variants" / "calling"
if not vdir.exists():
    vdir = rd / "06_variants"
if vdir.exists():
    for f in sorted(vdir.glob("*.tsv")):
        with f.open(newline="", encoding="utf-8") as fh:
            lines = fh.read().splitlines()
        if not lines:
            continue
        header = lines[0].split("\\t")
        for line in lines[1:]:
            cols = line.split("\\t")
            d = dict(zip(header, cols))
            d["sample"] = f.stem.replace("_variants","")
            var_rows.append(d)
if var_rows:
    fields = ["sample"] + [k for k in var_rows[0] if k != "sample"]
    write_csv("summary_variants.csv", var_rows, fields)

# --- Diversity ---
div_rows = []
ddir = rd / "06_variants" / "diversity"
if ddir.exists():
    for f in sorted(ddir.glob("*.tsv")):
        with f.open(newline="", encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\\t")
            for row in r:
                row["sample"] = f.stem.split(".")[0]
                div_rows.append(row)
if div_rows:
    fields = list(div_rows[0].keys())
    write_csv("summary_diversity.csv", div_rows, fields)

# --- Lineage SNPs ---
ls_rows = []
ldir = rd / "10_lineage" / "snps"
if ldir.exists():
    for f in sorted(ldir.glob("*.tsv")):
        with f.open(newline="", encoding="utf-8") as fh:
            r = csv.DictReader(fh, delimiter="\\t")
            for row in r:
                row["_file"] = f.name
                ls_rows.append(row)
if ls_rows:
    write_csv("summary_lineage_snps.csv", ls_rows, list(ls_rows[0].keys()))

# --- Sample QC ---
sq = rd / "08_coverage" / "sample_qc" / "sample_qc_summary.tsv"
if sq.exists():
    rows = []
    with sq.open(newline="", encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\\t")
        rows = list(r)
    if rows:
        write_csv("summary_sample_qc.csv", rows, list(rows[0].keys()))

print("Aggregation done")
PY

    # ========== 2) MultiQC native modules ==========
    run_mqc() {
        name=\$1; out=\$2; shift 2
        dirs=""
        for d in "\$@"; do [ -d "\$d" ] && dirs="\$dirs \$d"; done
        [ -z "\$dirs" ] && return 0
        multiqc -f -n multiqc_\${name}.html --title "MARV-GEN \$name" -o "\$out" \$dirs || true
    }

    run_mqc fastp      01_fastp      \$RD/01_fastp
    run_mqc host       02_host       \$RD/02_nonhuman_reads/stats \$RD/02_nonhuman_reads
    run_mqc taxonomy   03_taxonomy   \$RD/03_taxonomy
    run_mqc qualimap   05_qualimap   \$RD/05_mapping_qc
    run_mqc coverage   08_coverage   \$RD/08_coverage
    run_mqc nextclade  09_nextclade  \$RD/09_nextclade
    run_mqc pipeline   13_pipeline   \$RD/pipeline_info

    # MultiQC on summary CSVs dir (custom content if formatted)
    mkdir -p summaries_mqc
    cp -f summary_*.csv summaries_mqc/ 2>/dev/null || true
    multiqc -f -n multiqc_summaries.html --title "MARV-GEN Summaries" -o summaries_mqc summaries_mqc || true

    # ========== 3) HTML tables for key summaries ==========
    python3 - << 'PY'
from pathlib import Path
import csv, html

def csv_html(csv_path, title, out_html):
    p = Path(csv_path)
    if not p.exists():
        return
    out = Path(out_html)
    out.parent.mkdir(parents=True, exist_ok=True)
    with p.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return
    cols = list(rows[0].keys())
    parts = [f"<html><head><meta charset=utf-8><title>{title}</title>",
             "<style>body{font-family:sans-serif;margin:24px}table{border-collapse:collapse;width:100%;font-size:13px}",
             "th,td{border:1px solid #ccc;padding:4px 8px;text-align:left}th{background:#eee}</style></head><body>",
             f"<h1>{html.escape(title)}</h1><table><tr>" + "".join(f"<th>{html.escape(c)}</th>" for c in cols) + "</tr>"]
    for row in rows:
        parts.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(c,'')))}</td>" for c in cols) + "</tr>")
    out.write_text("\\n".join(parts) + "</table></body></html>")
    print("HTML", out)

csv_html("summary_kraken2.csv", "Kraken2 Summary", "03_taxonomy/multiqc_kraken2_table.html")
csv_html("summary_nextclade.csv", "Nextclade Summary", "09_nextclade/multiqc_nextclade_table.html")
csv_html("summary_variants.csv", "Variants Summary", "06_variants/multiqc_variants_table.html")
csv_html("summary_diversity.csv", "Diversity Summary", "06_variants/multiqc_diversity_table.html")
csv_html("summary_lineage_snps.csv", "Lineage SNPs", "10_lineage/multiqc_lineage_table.html")
csv_html("summary_sample_qc.csv", "Sample QC", "08_coverage/multiqc_sample_qc_table.html")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version 2>/dev/null | sed 's/multiqc, version //' || echo unknown)
    END_VERSIONS
    exit 0
    """
}
