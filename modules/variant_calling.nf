process VARIANT_CALLING {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/06_variants/calling", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("*_variants.tsv"), emit: tsv
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.*"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    samtools mpileup -A -d 1000000 -B -Q 0 -f ${reference} ${bam} \\
        | ivar variants -r ${reference} -p ${prefix}_variants -q ${params.ivar_min_quality}

    bcftools mpileup -Ou -f ${reference} ${bam} \\
        | bcftools call -mv -Oz --ploidy ${params.bcftools_ploidy} -o ${prefix}.vcf.gz
    bcftools index ${prefix}.vcf.gz || tabix -p vcf ${prefix}.vcf.gz || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | head -1 | sed 's/iVar version //')
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
