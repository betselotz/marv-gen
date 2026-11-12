process KRAKEN2_BUILD {
    label 'process_high'
    label 'process_long'
    storeDir "${projectDir}/database"

    output:
    path "kraken2_viral", emit: db
    path "versions.yml", emit: versions

    script:
    def url = params.kraken2_prebuilt_url ?: "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20250714.tar.gz"
    """
    set -euo pipefail
    mkdir -p kraken2_viral
    DB="\$PWD/kraken2_viral"

    if [ -f "\$DB/hash.k2d" ] && [ -f "\$DB/taxo.k2d" ]; then
        echo "DB already present under storeDir"
    else
        echo "=== Download prebuilt viral Kraken2 DB ==="
        OK=0
        if command -v wget >/dev/null 2>&1; then
            wget -c -O db.tar.gz "${url}" && OK=1 || OK=0
        else
            curl -L -C - -o db.tar.gz "${url}" && OK=1 || OK=0
        fi
        if [ "\$OK" = "1" ] && [ -s db.tar.gz ]; then
            tar -xzf db.tar.gz -C kraken2_viral --strip-components=1 2>/dev/null \\
              || tar -xzf db.tar.gz -C kraken2_viral
            rm -f db.tar.gz
            if [ ! -f "\$DB/hash.k2d" ]; then
                SUB=\$(find kraken2_viral -name hash.k2d | head -1)
                if [ -n "\$SUB" ]; then
                    DIR=\$(dirname "\$SUB")
                    mv "\$DIR"/* kraken2_viral/ 2>/dev/null || true
                fi
            fi
        fi

        if [ ! -f "\$DB/hash.k2d" ]; then
            echo "=== Fallback: kraken2-build viral from NCBI ==="
            kraken2-build --download-taxonomy --db "\$DB" --threads ${task.cpus}
            kraken2-build --download-library viral --db "\$DB" --threads ${task.cpus}
            kraken2-build --build --db "\$DB" --threads ${task.cpus}
            kraken2-build --clean --db "\$DB" || true
        fi
    fi

    test -f "\$DB/hash.k2d"
    test -f "\$DB/taxo.k2d"
    ls -lh "\$DB"/hash.k2d "\$DB"/taxo.k2d

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version 2>&1 | head -1 || echo unknown)
    END_VERSIONS
    """
}
