process MSA {
    label 'process_high'
    label 'process_long'
    publishDir "${params.outdir}/11_phylogeny/iqtree/msa", mode: params.publish_dir_mode

    input:
    path fastas
    path reference
    path compare_fastas

    output:
    path "alignment.fasta", emit: alignment
    path "versions.yml", emit: versions

    script:
    """
    cat ${reference} ${fastas} > combined.fasta
    if ls ${compare_fastas} 1>/dev/null 2>&1; then
        cat ${compare_fastas} >> combined.fasta || true
    fi
    mafft --auto --thread ${task.cpus} combined.fasta > alignment.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | head -1)
    END_VERSIONS
    """
}
