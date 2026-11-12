process HOST_REMOVAL {
    tag "$meta.id"
    label 'process_high'
    conda "${projectDir}/envs/host_removal_env.yaml"
    publishDir "${params.outdir}/02_nonhuman_reads", mode: params.publish_dir_mode,
        pattern: "*nonhost*.fastq.gz"
    publishDir "${params.outdir}/02_nonhuman_reads/stats", mode: params.publish_dir_mode,
        pattern: "*.stats.txt"
    input:
    tuple val(meta), path(r1), path(r2)
    val host_index_prefix
    output:
    tuple val(meta), path("*_nonhost_R1.fastq.gz"), path("*_nonhost_R2.fastq.gz"), emit: reads
    path "*.stats.txt", emit: stats
    path "versions.yml", emit: versions
    script:
    def prefix = meta.id
    """
    set -euo pipefail
    IDX="${host_index_prefix}"
    if [ ! -f "\${IDX}.1.bt2" ] && [ ! -f "\${IDX}.1.bt2l" ]; then
        echo "ERROR: Bowtie2 index not found at: \$IDX"
        ls -la "\$(dirname "\$IDX")" 2>/dev/null | head -20 || true
        exit 1
    fi
    if [ ! -s "${r1}" ] || [ ! -s "${r2}" ]; then
        echo "ERROR: missing/empty trimmed reads: ${r1} ${r2}"
        ls -la
        exit 1
    fi
    echo "Using index: \$IDX"
    echo "R1=${r1} R2=${r2}"
    bowtie2 \\
        -x "\$IDX" \\
        -1 ${r1} \\
        -2 ${r2} \\
        --un-conc-gz ${prefix}_tmp_R%.fastq.gz \\
        --threads ${task.cpus} \\
        -S /dev/null \\
        2> ${prefix}.host_removal.stats.txt
    mv ${prefix}_tmp_R1.fastq.gz ${prefix}_nonhost_R1.fastq.gz
    mv ${prefix}_tmp_R2.fastq.gz ${prefix}_nonhost_R2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    END_VERSIONS
    """
}
