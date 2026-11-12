process FOCUSED_TREETIME {
    label 'process_medium'
    conda "${projectDir}/envs/focused_treetime_env.yaml"
    publishDir "${params.outdir}/12_focused_treetime", mode: params.publish_dir_mode

    input:
    path alignment
    path metadata

    output:
    path "focused_JN_ET", emit: results
    path "versions.yml", emit: versions

    when:
    params.focused_treetime

    script:
    def tips = params.focused_tips ?: "JN408064.1,ET_MARV_45,ET_MARV_31,ET_MARV_32,ET_MARV_23,ET_MARV_60,ET_MARV_262,ET_MARV_261"
    """
    set -euo pipefail
    mkdir -p focused_JN_ET

    python3 ${projectDir}/bin/prepare_focused_treetime.py \\
        ${alignment} ${metadata} "${tips}" focused_JN_ET

    iqtree -s focused_JN_ET/JN_ET.fasta -m MFP -B 1000 -T AUTO -pre focused_JN_ET/JN_ET

    treetime --tree focused_JN_ET/JN_ET.treefile \\
             --aln focused_JN_ET/JN_ET.fasta \\
             --dates focused_JN_ET/dates.csv \\
             --reroot JN408064.1 \\
             --clock-filter 0 \\
             --outdir focused_JN_ET/treetime_out

    echo "=== Focused TreeTime dates ==="
    cat focused_JN_ET/treetime_out/dates.tsv || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree -version 2>&1 | head -1 || echo unknown)
        treetime: \$(treetime --version 2>&1 | head -1 || echo unknown)
    END_VERSIONS
    """
}
