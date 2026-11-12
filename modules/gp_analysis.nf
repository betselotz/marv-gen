process GP_ANALYSIS {
    label 'process_low'
    publishDir "${params.outdir}/10_lineage/gp_analysis", mode: params.publish_dir_mode

    input:
    path protein_faas
    path domains_csv
    path reference
    path gff
    path comparators_faa

    output:
    path "gp_similarity.tsv", emit: similarity
    path "gp_mr191_sites.tsv", emit: mr191
    path "gp_escape_sites.tsv", emit: escape
    path "gp_mutations_vs_musoke.tsv", emit: mutations
    path "gp_domain_summary.tsv", emit: domain_summary
    path "gp_quality.tsv", emit: quality
    path "gp_identity_matrix.tsv", emit: identity_matrix
    path "gp_summary.md", emit: summary
    path "gp_musoke_reference.faa", emit: musoke_gp
    path "gp_all_for_tree.faa", emit: alignment_faa
    path "versions.yml", emit: versions

    script:
    def xmax = params.gp_max_x_fraction ?: 0.15
    """
    set -euo pipefail
    COMP=""
    if [ -f "${comparators_faa}" ] && [ "${comparators_faa}" != "NO_FILE" ]; then
        COMP="${comparators_faa}"
    fi
    python3 ${projectDir}/bin/gp_analysis.py \\
        ${domains_csv} ${reference} ${gff} "\$COMP" ${xmax}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}

process GP_TREE {
    label 'process_medium'
    publishDir "${params.outdir}/10_lineage/gp_analysis/tree", mode: params.publish_dir_mode

    input:
    path gp_faa

    output:
    path "gp_aligned.fasta", emit: alignment
    path "gp_iqtree.treefile", emit: tree
    path "gp_iqtree*", emit: iqtree_files
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail
    python3 - << 'PY'
name = None
seq = []
recs = []
def flush():
    global name, seq
    if not name:
        return
    s = "".join(seq).upper()
    x = s.count("X") / max(len(s), 1)
    if x <= 0.15 and len(s) >= 600:
        recs.append((name, s))
    name = None
    seq = []
with open("${gp_faa}") as fh:
    for line in fh:
        if line.startswith(">"):
            flush()
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    flush()
with open("gp_for_align.faa", "w") as out:
    for n, s in recs:
        out.write(f">{n}\\n")
        for i in range(0, len(s), 60):
            out.write(s[i:i+60] + "\\n")
print("tree input seqs:", len(recs))
PY
    mafft --auto gp_for_align.faa > gp_aligned.fasta
    if command -v iqtree2 >/dev/null 2>&1; then
        IQT=iqtree2
    elif command -v iqtree >/dev/null 2>&1; then
        IQT=iqtree
    elif command -v iqtree3 >/dev/null 2>&1; then
        IQT=iqtree3
    else
        echo "ERROR: iqtree not found"; exit 1
    fi
    \$IQT -s gp_aligned.fasta -m MFP -B 1000 -T ${task.cpus} --prefix gp_iqtree
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(\$IQT --version 2>&1 | head -1)
        mafft: \$(mafft --version 2>&1 | head -1 || echo unknown)
    END_VERSIONS
    """
}
