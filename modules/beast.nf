process PREPARE_BEAST {
    label 'process_low'
    conda "${projectDir}/envs/python_analysis_env.yaml"
    publishDir "${params.outdir}/12_timetree/beast", mode: params.publish_dir_mode
    input:
    path alignment
    path metadata
    output:
    path "beast_alignment.fasta", emit: alignment
    path "beast_dates_decimal.tsv", emit: decimal
    path "beast_metadata_clean.csv", emit: meta_clean
    path "versions.yml", emit: versions
    script:
    """
    python3 ${projectDir}/bin/prepare_beast_dates.py \\
        ${alignment} ${metadata} \\
        beast_alignment.fasta beast_dates_decimal.tsv beast_metadata_clean.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}

process BEAST_RUN {
    tag "run${run_id}"
    label 'process_beast'
    conda "${projectDir}/envs/beast_env.yaml"
    publishDir "${params.outdir}/12_timetree/beast/runs", mode: params.publish_dir_mode
    input:
    path xml
    val run_id
    val seed
    output:
    tuple val(run_id), path("run${run_id}.log"), emit: log
    tuple val(run_id), path("run${run_id}.trees"), emit: trees
    path "versions.yml", emit: versions
    script:
    def chain = params.beast_chain_length ?: 50000000
    def logev = params.beast_log_every ?: 50000
    """
    set -euo pipefail
    cp -f ${xml} run_input.xml
    sed -i 's/chainLength="[0-9]*"/chainLength="${chain}"/g' run_input.xml
    sed -i 's/logEvery="[0-9]*"/logEvery="${logev}"/g' run_input.xml
    grep -E 'chainLength|logEvery' run_input.xml | head -5
    beast -beagle_SSE -beagle_threading_off -threads 1 -seed ${seed} -overwrite \
        -prefix run${run_id} run_input.xml 2>&1 | tee run${run_id}.stdout.log
    LOG=\$(ls -1 run${run_id}*.log 2>/dev/null | grep -v stdout | head -1 || true)
    TREES=\$(ls -1 run${run_id}*.trees 2>/dev/null | head -1 || true)
    if [ -z "\$LOG" ] || [ -z "\$TREES" ]; then
        echo "ERROR: missing BEAST outputs for run${run_id}"
        ls -la
        tail -80 run${run_id}.stdout.log || true
        exit 1
    fi
    [ "\$LOG" = "run${run_id}.log" ] || cp -f "\$LOG" run${run_id}.log
    [ "\$TREES" = "run${run_id}.trees" ] || cp -f "\$TREES" run${run_id}.trees
    NLINES=\$(wc -l < run${run_id}.log)
    echo "BEAST run${run_id}: lines=\$NLINES"
    if [ "\$NLINES" -lt 20 ]; then
        echo "ERROR: BEAST log almost empty"
        tail -80 run${run_id}.stdout.log || true
        exit 1
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beast: \$(beast -version 2>&1 | head -1 || echo unknown)
    END_VERSIONS
    """
}

process BEAST_POSTPROCESS {
    label 'process_medium'
    conda "${projectDir}/envs/beast_env.yaml"
    publishDir "${params.outdir}/12_timetree/beast/postprocess", mode: params.publish_dir_mode

    input:
    path run1_log
    path run2_log
    path run1_trees
    path run2_trees

    output:
    path "marburg_combined.trees", emit: combined_trees
    path "marburg_MCC.tree", emit: mcc
    path "marburg_combined.log", emit: combined_log, optional: true
    path "ESS_summary.txt", emit: ess
    path "TMRCA_summary.txt", emit: tmrca
    path "versions.yml", emit: versions

    script:
    // burn-in states (default 10% of 50M)
    def burnin = params.beast_burnin_states ?: (params.beast_burnin ?: 5000000)
    """
    set -euo pipefail

    echo "Inputs:"
    ls -lh ${run1_log} ${run2_log} ${run1_trees} ${run2_trees}

    # BEAST 1.10: output is LAST positional arg (no -o)
    logcombiner -trees -burnin ${burnin} \\
        ${run1_trees} ${run2_trees} marburg_combined.trees

    logcombiner -burnin ${burnin} \\
        ${run1_log} ${run2_log} marburg_combined.log || true

    treeannotator -burnin 0 -heights mean \\
        marburg_combined.trees marburg_MCC.tree

    {
        echo "BEAST postprocess"
        echo "burnin_states=${burnin}"
        echo "Open in Tracer for ESS (target > 200):"
        echo "  ${run1_log}"
        echo "  ${run2_log}"
        echo "  marburg_combined.log (if produced)"
        echo "MCC tree: marburg_MCC.tree"
        wc -l ${run1_log} ${run2_log} marburg_combined.trees 2>/dev/null || true
    } > ESS_summary.txt

    # Automatic TMRCA for Ethiopia-only and JN+Ethiopia
    python3 - << 'PY'
from Bio import Phylo
from pathlib import Path
tree = Phylo.read("marburg_MCC.tree", "nexus")
et_tips = [t for t in tree.get_terminals() if t.name and t.name.startswith("ET_MARV")]
jn = next((t for t in tree.get_terminals() if t.name and "JN408064" in (t.name or "")), None)
def height(clade):
    return max((tree.distance(clade, t) for t in clade.get_terminals()), default=0.0)
lines = []
if et_tips:
    et_mrca = tree.common_ancestor(et_tips)
    et_h = height(et_mrca)
    lines.append("Ethiopia_only_height_yr = %.4f" % et_h)
if et_tips and jn:
    jn_et_mrca = tree.common_ancestor(et_tips + [jn])
    jn_et_h = height(jn_et_mrca)
    lines.append("JN_plus_Ethiopia_height_yr = %.4f" % jn_et_h)
    lines.append("path_ET_MRCA_to_JN_ET_MRCA_yr = %.4f" % (jn_et_h - et_h))
if not lines:
    lines.append("no_ET_or_JN_tips_found")
text = chr(10).join(lines) + chr(10)
Path("TMRCA_summary.txt").write_text(text)
print(text, end="")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        logcombiner: \$(logcombiner -help 2>&1 | head -1 || echo beast-1.10)
        treeannotator: \$(treeannotator -help 2>&1 | head -1 || echo beast-1.10)
    END_VERSIONS
    """
}
