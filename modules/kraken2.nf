process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/03_taxonomy", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*.kraken2.report.txt"), emit: report
    tuple val(meta), path("*.kraken2.out.txt"), optional: true, emit: output
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --report ${prefix}.kraken2.report.txt \\
        --output ${prefix}.kraken2.out.txt \\
        --paired ${reads[0]} ${reads[1]} \\
        --gzip-compressed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
