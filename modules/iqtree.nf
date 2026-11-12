process IQTREE {
    label 'process_high'
    label 'process_long'
    publishDir "${params.outdir}/11_phylogeny/iqtree/iqtree", mode: params.publish_dir_mode

    input:
    path alignment

    output:
    path "*.treefile", emit: tree
    path "*.iqtree", emit: report
    path "*.log", emit: log
    path "versions.yml", emit: versions

    script:
    """
    # Prefer iqtree2 if present, else iqtree
    if command -v iqtree2 >/dev/null 2>&1; then
        IQTREE_BIN=iqtree2
    elif command -v iqtree >/dev/null 2>&1; then
        IQTREE_BIN=iqtree
    else
        echo "ERROR: iqtree/iqtree2 not found in PATH"
        exit 1
    fi

    \$IQTREE_BIN -s ${alignment} -m ${params.iqtree_model} \\
        -B ${params.iqtree_bootstrap} -T ${task.cpus} --prefix marv_iqtree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(\$IQTREE_BIN --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}
