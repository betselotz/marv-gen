process QUALIMAP {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/05_mapping_qc", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.id}", emit: report
    path "versions.yml", emit: versions

    script:
    """
    qualimap bamqc -bam ${bam} -outdir ${meta.id} -nt ${task.cpus} \\
        --java-mem-size=${task.memory.toGiga()}G
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: "2.3"
    END_VERSIONS
    """
}
