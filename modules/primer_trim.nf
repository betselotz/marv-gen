process PRIMER_TRIM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/02b_primer_trim", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed

    output:
    tuple val(meta), path("*.trimmed.bam"), path("*.trimmed.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    ivar trim -i ${bam} -b ${primer_bed} -p ${prefix}.trimmed -e -q 20 -m 30
    samtools sort -@ ${task.cpus} -o ${prefix}.trimmed.bam ${prefix}.trimmed.bam
    samtools index ${prefix}.trimmed.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | head -1 | sed 's/iVar version //')
    END_VERSIONS
    """
}
