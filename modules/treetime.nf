process TREETIME {
    conda "${projectDir}/envs/treetime_env.yaml"
    label 'process_medium'
    label 'process_long'
    publishDir "${params.outdir}/12_timetree/treetime", mode: params.publish_dir_mode

    input:
    path tree
    path alignment
    path metadata

    output:
    path "treetime_output", emit: results
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p treetime_output

    python3 - << 'PY'
import csv
from pathlib import Path
src = Path("${metadata}")
rows = []
with src.open(newline="", encoding="utf-8-sig") as fh:
    r = csv.DictReader(fh)
    fm = {k.lower().strip(): k for k in (r.fieldnames or [])}
    nk = fm.get("name") or fm.get("sample") or fm.get("strain")
    dk = fm.get("date")
    for row in r:
        name = (row.get(nk) or "").strip().replace("\\xa0", "")
        date = (row.get(dk) or "").strip()
        if name and date:
            rows.append({"name": name, "date": date})
with open("dates.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["name", "date"])
    w.writeheader()
    w.writerows(rows)
print(f"Wrote dates.csv with {len(rows)} tips")
PY

    if command -v treetime >/dev/null 2>&1; then
        TT=treetime
      elif python3 -c 'import treetime' 2>/dev/null; then
        TT="python3 -m treetime"
      else
        echo "ERROR: treetime not installed in conda env"
        conda list 2>/dev/null | head -50 || true
        exit 127
      fi
      \$TT --tree ${tree} --aln ${alignment} --dates dates.csv \\
        --clock-filter 0 --outdir treetime_output \\
        2>&1 | tee treetime_output/treetime.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        treetime: \$(treetime --version 2>&1 | head -1 || echo unknown)
    END_VERSIONS
    """
}
