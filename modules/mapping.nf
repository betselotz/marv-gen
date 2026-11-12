process MAPPING {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/04_mapping_results/bam", mode: params.publish_dir_mode, pattern: "*.{bam,bai}"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    minimap2 -ax ${params.minimap2_preset} -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} \\
        | samtools view -b -q 30 -F 4 - \\
        | samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam -
    samtools index ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
