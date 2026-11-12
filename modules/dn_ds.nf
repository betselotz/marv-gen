process DN_DS {
    label 'process_medium'
    publishDir "${params.outdir}/10_lineage/dn_ds", mode: params.publish_dir_mode

    input:
    path consensus_fastas
    path reference
    path gff

    output:
    path "selective_pressure.tsv", emit: results
    path "versions.yml", emit: versions

    script:
    """
    python3 - << 'PY'
from Bio import SeqIO, Align
from Bio.Seq import Seq
from pathlib import Path
import collections

ref_rec = list(SeqIO.parse("${reference}", "fasta"))[0]
ref = str(ref_rec.seq).upper()

# CDS from GFF
cds = []
gff = Path("${gff}")
if gff.exists():
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip(): continue
        p = line.split("\\t")
        if len(p) < 9 or p[2] != "CDS": continue
        start, end = int(p[3])-1, int(p[4])
        name = "CDS"
        for key in ("Name=", "gene=", "ID="):
            if key in p[8]:
                name = p[8].split(key)[1].split(";")[0].strip('"')
                break
        cds.append((name, start, end, p[6]))

seqs = {}
for f in Path(".").glob("*.fa"):
    r = next(SeqIO.parse(f, "fasta"))
    seqs[r.id] = str(r.seq).upper()

def count_ns(ref_codon, alt_codon):
    if "N" in ref_codon or "N" in alt_codon or "-" in ref_codon+alt_codon:
        return 0, 0
    if len(ref_codon) != 3 or len(alt_codon) != 3:
        return 0, 0
    try:
        ra = str(Seq(ref_codon).translate())
        aa = str(Seq(alt_codon).translate())
    except Exception:
        return 0, 0
    n_diff = sum(a!=b for a,b in zip(ref_codon, alt_codon))
    if n_diff == 0:
        return 0, 0
    if ra == aa:
        return 0, n_diff  # synonymous
    return n_diff, 0  # nonsyn

rows = []
for gene, start, end, strand in cds:
    N_tot = S_tot = 0
    sites = 0
    for sid, seq in seqs.items():
        if end > len(seq) or end > len(ref):
            continue
        for i in range(start, end-2, 3):
            rc = ref[i:i+3]
            ac = seq[i:i+3]
            if len(rc) < 3 or len(ac) < 3:
                continue
            n, s = count_ns(rc, ac)
            N_tot += n; S_tot += s; sites += 1
    dn = N_tot / max(sites,1)
    ds = S_tot / max(sites,1)
    omega = (dn/ds) if ds > 0 else float("nan")
    rows.append((gene, N_tot, S_tot, f"{dn:.6f}", f"{ds:.6f}", f"{omega:.4f}" if omega==omega else "NA", sites))

with open("selective_pressure.tsv","w") as fh:
    fh.write("gene\\tN\\tS\\tdN\\tdS\\tomega\\tsites\\n")
    for r in rows:
        fh.write("\\t".join(map(str,r))+"\\n")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
