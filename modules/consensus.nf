process CONSENSUS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/07_consensus", mode: params.publish_dir_mode, pattern: "*.fa"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    samtools mpileup -A -d 0 -Q ${params.ivar_min_quality} -f ${reference} ${bam} \\
        | ivar consensus -p ${prefix} -q ${params.ivar_min_quality} \\
            -t ${params.ivar_min_frequency} -m ${params.ivar_min_depth} -n N
    if [ -f ${prefix}.fa ]; then sed -i "1s/.*/>${prefix}/" ${prefix}.fa; fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | head -1 | sed 's/iVar version //')
    END_VERSIONS
    """
}
