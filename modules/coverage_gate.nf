process COVERAGE_GATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/05_mapping_qc/coverage_gate", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path(bam), path(bai), path("${meta.id}.coverage_gate.tsv"), emit: all
    path "${meta.id}.coverage_gate.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    // Breadth = % reference bases with depth >= gate_breadth_min_depth (default 1 = any coverage)
    def min_pct   = params.gate_min_cov_pct ?: 50
    def min_depth = params.gate_breadth_min_depth ?: 1
    """
    set -euo pipefail

    # Breadth: fraction of ref bases with depth >= min_depth (usually 1)
    STATS=\$(samtools depth -a ${bam} | awk -v d=${min_depth} '
        { n++; s+=\$3; if (\$3+0 >= d) c++ }
        END {
            if (n==0) { print "0\\t0\\t0"; exit }
            printf "%.4f\\t%.2f\\t%d", 100.0*c/n, s/n, n
        }
    ')
    BREADTH=\$(echo "\$STATS" | cut -f1)
    MEAN=\$(echo "\$STATS" | cut -f2)
    NBASES=\$(echo "\$STATS" | cut -f3)

    # Also report % at 20x for information only (not used for PASS/FAIL)
    PCT20=\$(samtools depth -a ${bam} | awk '
        { n++; if (\$3+0 >= 20) c++ }
        END { if (n==0) print 0; else printf "%.4f", 100.0*c/n }
    ')

    STATUS="PASS"
    awk -v p="\$BREADTH" -v m="${min_pct}" 'BEGIN { if (p+0 < m+0) exit 1; exit 0 }' || STATUS="FAIL"

    echo -e "sample\\tbreadth_pct\\tmean_depth\\tpct_bases_ge_20x\\tthreshold_breadth_pct\\tbreadth_min_depth\\tn_ref_bases\\tstatus" > ${meta.id}.coverage_gate.tsv
    echo -e "${meta.id}\\t\${BREADTH}\\t\${MEAN}\\t\${PCT20}\\t${min_pct}\\t${min_depth}\\t\${NBASES}\\t\${STATUS}" >> ${meta.id}.coverage_gate.tsv

    echo "COVERAGE_GATE ${meta.id}: breadth=\${BREADTH}% (>=${min_depth}x) mean=\${MEAN} pct20=\${PCT20}% -> \${STATUS}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
