process INTRA_HOST_DIVERSITY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/06_variants/diversity", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(ivar_tsv)

    output:
    tuple val(meta), path("*.minor_variants.tsv"), emit: minors
    tuple val(meta), path("*.diversity_summary.tsv"), emit: summary

    script:
    def prefix = meta.id
    def maf = params.minor_allele_freq
    """
    python3 - << 'PY'
import csv
from pathlib import Path

tsv = Path("${ivar_tsv}")
maf = float("${maf}")
rows = []
with tsv.open() as fh:
    reader = csv.DictReader(fh, delimiter="\\t")
    for r in reader:
        try:
            freq = float(r.get("ALT_FREQ", r.get("ALT_FREQ", 0)) or 0)
            depth = float(r.get("TOTAL_DP", r.get("ALT_DP", 0)) or 0)
        except Exception:
            continue
        if maf <= freq < 0.5 and depth >= 10:
            rows.append(r)

out = Path("${prefix}.minor_variants.tsv")
if rows:
    with out.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\\t")
        w.writeheader()
        w.writerows(rows)
else:
    out.write_text("REGION\\tPOS\\tREF\\tALT\\tALT_FREQ\\tTOTAL_DP\\n")

n_minors = len(rows)
mean_freq = sum(float(r.get("ALT_FREQ",0) or 0) for r in rows)/n_minors if n_minors else 0
with open("${prefix}.diversity_summary.tsv","w") as fh:
    fh.write("sample\\tn_minor_variants\\tmean_minor_freq\\tmaf_threshold\\n")
    fh.write(f"${prefix}\\t{n_minors}\\t{mean_freq:.4f}\\t{maf}\\n")
PY
    """
}
