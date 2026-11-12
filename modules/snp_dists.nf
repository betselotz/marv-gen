process SNP_DISTS {
    label 'process_medium'
    publishDir "${params.outdir}/11_phylogeny/iqtree/snp_dists", mode: params.publish_dir_mode

    input:
    path alignment

    output:
    path "snp_distances.tsv", emit: matrix
    path "snp_clusters.tsv", emit: clusters
    path "versions.yml", emit: versions

    script:
    def thr = params.snp_cluster_threshold ?: 5
    """
    snp-dists -c ${alignment} > snp_distances.tsv

    python3 - << 'PY'
import csv
thr = int("${thr}")
with open("snp_distances.tsv") as fh:
    reader = csv.reader(fh, delimiter="\\t")
    header = next(reader)
    samples = header[1:]
    mat = {row[0]: {samples[i]: int(row[i+1]) for i in range(len(samples))} for row in reader}

# simple single-linkage clusters
parent = {s: s for s in samples}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
def union(a,b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[rb] = ra

for a in samples:
    for b in samples:
        if a < b and mat.get(a,{}).get(b, 999) <= thr:
            union(a,b)

clusters = {}
for s in samples:
    clusters.setdefault(find(s), []).append(s)
with open("snp_clusters.tsv","w") as fh:
    fh.write("cluster\\tn_samples\\tsamples\\n")
    for i, (k, members) in enumerate(sorted(clusters.items(), key=lambda x: -len(x[1])), 1):
        fh.write(f"C{i}\\t{len(members)}\\t{','.join(sorted(members))}\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snp-dists: \$(snp-dists -v 2>&1 | head -1 || echo "0.8.2")
    END_VERSIONS
    """
}
