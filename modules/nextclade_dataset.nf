process NEXTCLADE_DATASET {
    label 'process_low'
    storeDir "${params.nextclade_dataset_dir}"

    output:
    path "marburg", emit: dataset
    path "versions.yml", emit: versions

    script:
    def name = params.nextclade_dataset_name
    """
    set -euo pipefail
    mkdir -p marburg

    echo "Downloading Nextclade dataset: ${name}"
    nextclade dataset get \\
        --name "${name}" \\
        --output-dir marburg

    if [ ! -f marburg/pathogen.json ] && [ ! -f marburg/reference.fasta ]; then
        echo "ERROR: dataset incomplete"
        ls -la marburg/
        exit 1
    fi
    ls -lh marburg/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version 2>&1 | head -1 | sed 's/nextclade //' || echo unknown)
        dataset: ${name}
    END_VERSIONS
    """
}
