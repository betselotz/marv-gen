process LINEAGE_SNPS {
    label 'process_medium'
    publishDir "${params.outdir}/10_lineage/snps", mode: params.publish_dir_mode

    input:
    path consensus_fastas
    path reference
    path gff

    output:
    path "lineage_snps.tsv", emit: snps
    path "lineage_snp_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python3 - << 'PY'
from Bio import SeqIO
from pathlib import Path
import collections

ref = str(list(SeqIO.parse("${reference}", "fasta"))[0].seq).upper()
seqs = {}
for f in Path(".").glob("*.fa"):
    rec = next(SeqIO.parse(f, "fasta"))
    seqs[rec.id] = str(rec.seq).upper()

# positions where any sample differs from ref
snp_rows = []
for pos in range(min(len(ref), min(len(s) for s in seqs.values()) if seqs else 0)):
    ref_base = ref[pos]
    if ref_base == "N":
        continue
    alleles = {sid: seq[pos] for sid, seq in seqs.items() if pos < len(seq) and seq[pos] not in "N-"}
    if not alleles:
        continue
    alt_alleles = {a for a in alleles.values() if a != ref_base}
    if not alt_alleles:
        continue
    for sid, base in alleles.items():
        if base != ref_base:
            snp_rows.append((sid, pos+1, ref_base, base))

with open("lineage_snps.tsv","w") as fh:
    fh.write("sample\\tpos\\tref\\talt\\n")
    for r in snp_rows:
        fh.write("\\t".join(map(str,r))+"\\n")

# summary: positions shared by >=2 samples
pos_count = collections.Counter(r[1] for r in snp_rows)
with open("lineage_snp_summary.tsv","w") as fh:
    fh.write("pos\\tref\\tn_samples\\tsamples\\n")
    by_pos = collections.defaultdict(list)
    for sid, pos, refb, alt in snp_rows:
        by_pos[pos].append((sid, refb, alt))
    for pos in sorted(by_pos):
        samples = by_pos[pos]
        fh.write(f"{pos}\\t{samples[0][1]}\\t{len(samples)}\\t{','.join(s[0] for s in samples)}\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
