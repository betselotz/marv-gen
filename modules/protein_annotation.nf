process PROTEIN_ANNOTATION {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/06_variants/protein", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(consensus)
    path gff
    path reference

    output:
    tuple val(meta), path("*.proteins.faa"), emit: proteins
    tuple val(meta), path("*.cds.fna"), emit: cds
    path "*.protein_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    python3 << 'PY'
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import re

sample = "${prefix}"
cons = str(SeqIO.read("${consensus}", "fasta").seq).upper().replace("U", "T")
ref  = str(SeqIO.read("${reference}", "fasta").seq).upper().replace("U", "T")
gff  = Path("${gff}")

features = []
with gff.open() as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\\n").split("\\t")
        if len(p) < 9:
            continue
        ftype = p[2]
        if ftype not in ("CDS", "gene", "mature_peptide"):
            continue
        start, end = int(p[3]), int(p[4])
        strand = p[6]
        attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
        name = attrs.get("gene") or attrs.get("Name") or attrs.get("product") or attrs.get("ID") or ftype
        product = attrs.get("product") or name
        features.append(dict(type=ftype, start=start, end=end, strand=strand, name=name, product=product))

cds_feats = [f for f in features if f["type"] == "CDS"] or [f for f in features if f["type"] == "gene"]

prot_records, cds_records, rows = [], [], []
for f in cds_feats:
    s, e = f["start"] - 1, f["end"]
    if s < 0 or s >= len(cons):
        continue
    e = min(e, len(cons))
    seq = cons[s:e]
    if f["strand"] == "-":
        seq = str(Seq(seq).reverse_complement())
    aa = str(Seq(seq).translate(to_stop=False)).replace("*", "X")
    gene = re.sub(r"[^A-Za-z0-9_.-]", "_", f["name"])
    pid = "%s|%s|%d-%d" % (sample, gene, f["start"], f["end"])
    prot_records.append(">%s %s\\n%s" % (pid, f["product"], "\\n".join(aa[j:j+60] for j in range(0, len(aa), 60))))
    cds_records.append(">%s\\n%s" % (pid, "\\n".join(seq[j:j+60] for j in range(0, len(seq), 60))))
    rseq = ref[s:min(e, len(ref))]
    if f["strand"] == "-":
        rseq = str(Seq(rseq).reverse_complement())
    raa = str(Seq(rseq).translate(to_stop=False)).replace("*", "X") if rseq else ""
    # Only count changes where both AAs come from clean codons (no N in codon)
    diffs = []
    for pos in range(min(len(aa), len(raa))):
        a, b = aa[pos], raa[pos]
        if a in "X*-" or b in "X*-":
            continue
        if a != b:
            diffs.append("%s%d%s" % (b, pos+1, a))
    # Fraction of protein that is callable (not X)
    callable_aa = sum(1 for a in aa if a not in "X*-")
    rows.append((sample, f["name"], f["product"], f["start"], f["end"], f["strand"], len(seq), len(aa), callable_aa, len(diffs), ",".join(diffs[:50])))

Path("%s.proteins.faa" % sample).write_text("\\n".join(prot_records) + ("\\n" if prot_records else ""))
Path("%s.cds.fna" % sample).write_text("\\n".join(cds_records) + ("\\n" if cds_records else ""))
with open("%s.protein_summary.tsv" % sample, "w") as out:
    out.write("sample\\tgene\\tproduct\\tstart\\tend\\tstrand\\tcds_len\\taa_len\\tn_aa_changes\\taa_changes\\n")
    for r in rows:
        out.write("\\t".join(map(str, r)) + "\\n")
print("%s: %d features" % (sample, len(rows)))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c 'import Bio; print(Bio.__version__)')
    END_VERSIONS
    """
}
