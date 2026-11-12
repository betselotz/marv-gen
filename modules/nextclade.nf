process NEXTCLADE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/09_nextclade", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(consensus)
    path dataset

    output:
    tuple val(meta), path("*.nextclade.tsv"), emit: tsv
    path "*.nextclade.json", optional: true, emit: json
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    set -euo pipefail

    DS="${dataset}"
    if [ -d "\$DS/marburg" ]; then DS="\$DS/marburg"; fi
    if [ ! -f "\$DS/pathogen.json" ] && [ ! -f "\$DS/reference.fasta" ]; then
        echo "ERROR: invalid Nextclade dataset at \$DS"
        find "${dataset}" -maxdepth 3 -type f | head -40
        exit 1
    fi

    nextclade run \\
        --input-dataset "\$DS" \\
        --output-tsv ${prefix}.nextclade.tsv \\
        --output-json ${prefix}.nextclade.json \\
        ${consensus}

    head -3 ${prefix}.nextclade.tsv || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1 | head -1 | sed 's/nextclade //' || echo unknown)
    END_VERSIONS
    """
}
