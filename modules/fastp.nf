process FASTP {
    tag "$meta.id"
    label 'process_medium'
    conda "${projectDir}/envs/fastp_env.yaml"
    publishDir "${params.outdir}/01_fastp", mode: params.publish_dir_mode,
        pattern: "*.{html,json,log,txt}"
    publishDir "${params.outdir}/01_fastp/trimmed", mode: params.publish_dir_mode,
        pattern: "*trimmed*.fastq.gz"
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("*_R1.trimmed.fastq.gz"), path("*_R2.trimmed.fastq.gz"), emit: reads, optional: true
    path "*.json",  emit: json,   optional: true
    path "*.html",  emit: html,   optional: true
    path "*.log",   emit: log,    optional: true
    path "*qc_status.txt", emit: status, optional: true
    path "versions.yml", emit: versions
    script:
    def prefix = meta.id
    """
    set +e
    set -u
    PREFIX="${prefix}"
    R1="${reads[0]}"
    R2="${reads[1]}"
    LOG="\${PREFIX}.fastp.log"
    : > "\$LOG"
    STATUS="OK"
    REASON="."
    IN1="\$R1"
    IN2="\$R2"

    is_readable_fastq() {
        f="\$1"
        if echo "\$f" | grep -qi '\\.gz\$'; then
            gzip -t "\$f" 2>/dev/null && return 0
            return 1
        fi
        # plain fastq: must have at least one record
        [ -s "\$f" ] || return 1
        head -1 "\$f" 2>/dev/null | grep -q '^@'
    }

    salvage_gz() {
        in="\$1"
        out="\$2"
        ( gzip -dc "\$in" 2>/dev/null || true ) | awk '
            BEGIN { n=0 }
            {
              buf[n % 4] = \$0
              n++
              if (n % 4 == 0) {
                if (buf[0] ~ /^@/ && buf[2] ~ /^\\+/) {
                  print buf[0]
                  print buf[1]
                  print buf[2]
                  print buf[3]
                }
              }
            }
        ' | gzip -c > "\$out"
    }

    count_fq_gz() {
        gzip -dc "\$1" 2>/dev/null | awk 'END { print int(NR/4) }'
    }

    need_recover=0
    is_readable_fastq "\$R1" || need_recover=1
    is_readable_fastq "\$R2" || need_recover=1

    if [ "\$need_recover" -eq 1 ]; then
        echo "WARN: input not fully readable for \${PREFIX} — attempting recovery" | tee -a "\$LOG"
        STATUS="RECOVERED"
        # normalize each mate to gz
        if echo "\$R1" | grep -qi '\\.gz\$'; then
            salvage_gz "\$R1" "\${PREFIX}_R1.salvaged.fastq.gz"
        else
            cat "\$R1" 2>/dev/null | gzip -c > "\${PREFIX}_R1.salvaged.fastq.gz"
        fi
        if echo "\$R2" | grep -qi '\\.gz\$'; then
            salvage_gz "\$R2" "\${PREFIX}_R2.salvaged.fastq.gz"
        else
            cat "\$R2" 2>/dev/null | gzip -c > "\${PREFIX}_R2.salvaged.fastq.gz"
        fi
        N1=\$(count_fq_gz "\${PREFIX}_R1.salvaged.fastq.gz")
        N2=\$(count_fq_gz "\${PREFIX}_R2.salvaged.fastq.gz")
        N1=\${N1:-0}
        N2=\${N2:-0}
        if [ "\$N1" -lt "\$N2" ]; then N="\$N1"; else N="\$N2"; fi
        echo "Recovered R1=\$N1 R2=\$N2 paired_to=\$N" | tee -a "\$LOG"
        if [ "\$N" -le 0 ]; then
            echo "ERROR: zero recoverable pairs for \${PREFIX}" | tee -a "\$LOG"
            printf "sample\\tstatus\\treason\\n%s\\tFAIL\\tzero_recoverable\\n" "\$PREFIX" > "\${PREFIX}.qc_status.txt"
            cat > versions.yml << VERS
"${task.process}":
    fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //' || echo unknown)
VERS
            exit 0
        fi
        gzip -dc "\${PREFIX}_R1.salvaged.fastq.gz" 2>/dev/null | awk -v n="\$N" 'NR <= n*4' | gzip -c > "\${PREFIX}_R1.rec.fastq.gz"
        gzip -dc "\${PREFIX}_R2.salvaged.fastq.gz" 2>/dev/null | awk -v n="\$N" 'NR <= n*4' | gzip -c > "\${PREFIX}_R2.rec.fastq.gz"
        IN1="\${PREFIX}_R1.rec.fastq.gz"
        IN2="\${PREFIX}_R2.rec.fastq.gz"
        REASON="gzip_recovered_pairs=\${N}"
    fi

    fastp \\
        -i "\$IN1" \\
        -I "\$IN2" \\
        -o "\${PREFIX}_R1.trimmed.fastq.gz" \\
        -O "\${PREFIX}_R2.trimmed.fastq.gz" \\
        -q ${params.fastp_qualified_quality ?: 20} \\
        -l ${params.fastp_min_length ?: 50} \\
        --thread ${task.cpus} \\
        -j "\${PREFIX}.fastp.json" \\
        -h "\${PREFIX}.fastp.html" \\
        2>> "\$LOG"
    fp_rc=\$?
    if [ "\$fp_rc" -ne 0 ]; then
        echo "ERROR: fastp failed rc=\$fp_rc" | tee -a "\$LOG"
        printf "sample\\tstatus\\treason\\n%s\\tFAIL\\tfastp_rc=%s\\n" "\$PREFIX" "\$fp_rc" > "\${PREFIX}.qc_status.txt"
        cat > versions.yml << VERS
"${task.process}":
    fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //' || echo unknown)
VERS
        exit 0
    fi

    printf "sample\\tstatus\\treason\\n%s\\t%s\\t%s\\n" "\$PREFIX" "\$STATUS" "\$REASON" > "\${PREFIX}.qc_status.txt"
    cat > versions.yml << VERS
"${task.process}":
    fastp: \$(fastp --version 2>&1 | sed -e 's/fastp //' || echo unknown)
VERS
    exit 0
    """
}
