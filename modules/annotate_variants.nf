process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/06_variants/annotated", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path gff
    path reference

    output:
    tuple val(meta), path("*.annotated.tsv"), emit: tsv
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    # Gene-level annotation from GFF using bcftools + simple awk mapping
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\n' ${vcf} > ${prefix}.raw.tsv

    python3 - << 'PY'
import sys
from pathlib import Path

# Parse GFF for CDS features
genes = []
gff = Path("${gff}")
if gff.exists():
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\\t")
        if len(parts) < 9:
            continue
        if parts[2] not in ("CDS", "gene", "mature_peptide"):
            continue
        start, end = int(parts[3]), int(parts[4])
        attrs = parts[8]
        name = "unknown"
        for key in ("Name=", "gene=", "ID=", "product="):
            if key in attrs:
                name = attrs.split(key)[1].split(";")[0].strip('"')
                break
        genes.append((start, end, name, parts[6]))

def annotate(pos, ref, alt):
    pos = int(pos)
    gene, region = "intergenic", "noncoding"
    for start, end, name, strand in genes:
        if start <= pos <= end:
            gene = name
            # simple synonymous heuristic: length-preserving SNP in CDS
            if len(ref) == 1 and len(alt) == 1:
                region = "missense_or_synonymous"
            else:
                region = "indel"
            break
    return gene, region

out = open("${prefix}.annotated.tsv", "w")
out.write("sample\\tchrom\\tpos\\tref\\talt\\tqual\\tgene\\tregion\\n")
with open("${prefix}.raw.tsv") as fh:
    for line in fh:
        chrom, pos, ref, alt, qual = line.strip().split("\\t")
        gene, region = annotate(pos, ref, alt)
        out.write(f"${prefix}\\t{chrom}\\t{pos}\\t{ref}\\t{alt}\\t{qual}\\t{gene}\\t{region}\\n")
out.close()
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
