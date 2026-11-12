process SAMPLE_QC {
    label 'process_low'
    publishDir "${params.outdir}/08_coverage/sample_qc", mode: params.publish_dir_mode

    input:
    path coverage_tsvs

    output:
    path "sample_qc_summary.tsv", emit: summary
    path "passed_samples.txt", emit: passed
    path "failed_samples.txt", emit: failed

    script:
    def min_depth = params.min_mean_depth
    def min_cov   = params.min_coverage_pct
    def max_n     = params.max_n_pct
    """
    set +e
    echo -e "Sample\\tMeanDepth\\tCov10x\\tPercent_N\\tStatus\\tReason" > sample_qc_summary.tsv
    : > passed_samples.txt
    : > failed_samples.txt

    for f in ${coverage_tsvs}; do
        [ -f "\$f" ] || continue
        line=\$(tail -n +2 "\$f" | head -1)
        [ -z "\$line" ] && continue

        sample=\$(echo "\$line" | cut -f1)
        mean=\$(echo "\$line" | cut -f2)
        cov10=\$(echo "\$line" | cut -f7)
        pctn=\$(echo "\$line" | cut -f13)

        status="PASS"
        reason="."

        # numeric checks (do not abort script on non-zero)
        awk -v m="\$mean" -v c="\$cov10" -v n="\$pctn" \\
            -v min_d="${min_depth}" -v min_c="${min_cov}" -v max_n="${max_n}" 'BEGIN {
            if ((m+0) < (min_d+0)) exit 1
            if ((c+0) < (min_c+0)) exit 2
            if ((n+0) > (max_n+0)) exit 3
            exit 0
        }'
        rc=\$?

        if [ "\$rc" -eq 1 ]; then status="FAIL"; reason="mean_depth<${min_depth}"; fi
        if [ "\$rc" -eq 2 ]; then status="FAIL"; reason="cov10x<${min_cov}%"; fi
        if [ "\$rc" -eq 3 ]; then status="FAIL"; reason="percent_N>${max_n}%"; fi

        echo -e "\${sample}\\t\${mean}\\t\${cov10}\\t\${pctn}\\t\${status}\\t\${reason}" >> sample_qc_summary.tsv
        if [ "\$status" = "PASS" ]; then
            echo "\$sample" >> passed_samples.txt
        else
            echo "\$sample" >> failed_samples.txt
        fi
    done

    exit 0
    """
}
