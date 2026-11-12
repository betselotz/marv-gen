#!/usr/bin/env python3
import sys, csv, re
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord

MR191_SITES = sorted(set(
    list(range(62, 65)) + [67] + list(range(70, 73)) + [91, 95] +
    list(range(97, 101)) + [123] + list(range(125, 132)) + [154]
))
NPC1_RBS_SITES = sorted(set([129] + list(range(125, 132)) + [91, 95, 97, 98, 99, 100]))
ESCAPE = {
    129: ("N", "S", "NPC1 RBS / N129S"),
    220: ("S", "P", "glycan cap / S220P"),
    226: ("C", "Y", "glycan cap / C226Y"),
    455: ("P", "L", "mucin-like / P455L"),
}

def load_domains(path):
    domains = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            domains.append({
                "name": row["Name"].strip().strip('"'),
                "type": row["Type"].strip(),
                "start": int(row["Minimum"]),
                "end": int(row["Maximum"]),
            })
    return domains

def domain_of(domains, pos):
    hits = [d["name"] for d in domains if d["start"] <= pos <= d["end"]]
    return ";".join(hits) if hits else "other"

def in_core_ss(domains, pos):
    return any(d["type"] in ("helix", "strand") and d["start"] <= pos <= d["end"] for d in domains)

def extract_musoke_gp(reference, gff):
    ref = str(SeqIO.read(reference, "fasta").seq).upper().replace("U", "T")
    candidates = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip().split("\t")
            if len(p) < 9 or p[2] not in ("CDS", "gene"):
                continue
            attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            blob = " ".join([attrs.get("gene",""), attrs.get("Name",""), attrs.get("product",""), attrs.get("ID","")]).upper()
            if "GP" not in blob and "GLYCOPROTEIN" not in blob:
                continue
            start, end, strand = int(p[3]), int(p[4]), p[6]
            nuc = ref[start-1:end]
            if strand == "-":
                nuc = str(Seq(nuc).reverse_complement())
            if len(nuc) % 3:
                nuc = nuc[:len(nuc)//3*3]
            aa = str(Seq(nuc).translate(to_stop=True)).replace("*", "")
            candidates.append((abs(len(aa)-681), len(aa), start, end, strand, aa, blob))
    candidates.sort()
    if candidates and 500 <= candidates[0][1] <= 750:
        return candidates[0][5], candidates[0][2], candidates[0][3], candidates[0][6]
    start, end = 5941, 7986
    nuc = ref[start-1:end]
    if len(nuc) % 3:
        nuc = nuc[:len(nuc)//3*3]
    aa = str(Seq(nuc).translate(to_stop=True)).replace("*", "")
    return aa, start, end, "fallback_5941_7986"

def load_sample_gps(workdir):
    samples = {}
    for faa in Path(workdir).glob("*.proteins.faa"):
        sample = faa.name.replace(".proteins.faa", "")
        for rec in SeqIO.parse(faa, "fasta"):
            rid = rec.id.upper()
            desc = (rec.description or "").upper()
            if "|GP|" in rid or "GLYCOPROTEIN" in desc or re.search(r"(^|[|_])GP([|_]|$)", rid):
                seq = str(rec.seq).replace("*", "")
                while seq.endswith("X"):
                    seq = seq[:-1]
                samples[sample] = seq
                break
    return samples

def load_comparators(path):
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return {}
    return {r.id.split()[0]: str(r.seq).replace("*","").upper() for r in SeqIO.parse(path, "fasta")}

def align_map(ref_aa, qry_aa):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aln = aligner.align(ref_aa, qry_aa)[0]
    try:
        a, b = aln[0], aln[1]
    except Exception:
        parts = str(aln).split("\n")
        a, b = parts[0], parts[2]
    matches = tot = ref_pos = 0
    pos_map = {}
    for x, y in zip(a, b):
        if x != "-":
            ref_pos += 1
        if x != "-" and y != "-":
            tot += 1
            if x == y:
                matches += 1
            pos_map[ref_pos] = y
        elif x != "-":
            pos_map[ref_pos] = "-"
    ident = 100.0 * matches / tot if tot else 0.0
    return ident, pos_map, matches, tot

def write_tsv(path, rows):
    if not rows:
        Path(path).write_text("")
        return
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

def main():
    domains_csv, reference, gff = sys.argv[1], sys.argv[2], sys.argv[3]
    comparators_faa = sys.argv[4] if len(sys.argv) > 4 else ""
    max_x_frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.15
    domains = load_domains(domains_csv)
    musoke, gstart, gend, gsrc = extract_musoke_gp(reference, gff)
    print("Musoke GP len=%d %d-%d %s start=%s" % (len(musoke), gstart, gend, gsrc, musoke[:20]))
    Path("gp_musoke_reference.faa").write_text(
        ">Musoke_GP\n" + "\n".join(musoke[i:i+60] for i in range(0, len(musoke), 60)) + "\n"
    )
    samples = load_sample_gps(".")
    comps = load_comparators(comparators_faa)
    all_gps = {"Musoke_ref": musoke}
    all_gps.update(comps)
    all_gps.update(samples)
    sim_rows, mut_rows, mr191_rows, esc_rows = [], [], [], []
    domain_count_rows, quality_rows = [], []
    faa_for_tree = [SeqRecord(Seq(musoke), id="Musoke_ref", description="ref")]

    for sample, gp in sorted(samples.items()):
        fracX = gp.count("X") / max(len(gp), 1)
        qc = "PASS" if fracX <= max_x_frac and len(gp) >= 500 else "WARN_low_quality_GP"
        quality_rows.append({"sample": sample, "gp_len": len(gp), "fraction_X": "%.4f" % fracX,
            "max_x_frac_threshold": max_x_frac, "qc_flag": qc, "gp_start": gp[:15], "musoke_start": musoke[:15]})
        ident, pos_map, matches, tot = align_map(musoke, gp)
        sim_rows.append({"sample": sample, "comparator": "Musoke_ref", "gp_len": len(gp),
            "comparator_len": len(musoke), "pct_identity": "%.2f" % ident, "aligned_matches": matches,
            "aligned_columns": tot, "fraction_X": "%.4f" % fracX, "qc_flag": qc})
        for cname, cseq in comps.items():
            ident_c, _, m_c, t_c = align_map(cseq, gp)
            sim_rows.append({"sample": sample, "comparator": cname, "gp_len": len(gp),
                "comparator_len": len(cseq), "pct_identity": "%.2f" % ident_c, "aligned_matches": m_c,
                "aligned_columns": t_c, "fraction_X": "%.4f" % fracX, "qc_flag": qc})
        bucket = defaultdict(int)
        n_mut = 0
        for pos in range(1, len(musoke)+1):
            ra, sa = musoke[pos-1], pos_map.get(pos, "-")
            if sa in ("-", "X") or ra == "X" or sa == ra:
                continue
            n_mut += 1
            dom = domain_of(domains, pos)
            core = in_core_ss(domains, pos)
            mut_rows.append({"sample": sample, "aa_pos": pos, "musoke": ra, "sample_aa": sa,
                "mutation": "%s%d%s" % (ra, pos, sa), "domains": dom,
                "in_mr191_site": "yes" if pos in MR191_SITES else "no",
                "in_npc1_rbs_region": "yes" if pos in NPC1_RBS_SITES else "no",
                "in_core_helix_strand": "yes" if core else "no", "qc_flag": qc})
            if core: bucket["core_helix_strand"] += 1
            if "mucin" in dom.lower(): bucket["mucin"] += 1
            if "signal" in dom.lower(): bucket["signal"] += 1
            if "wing" in dom.lower(): bucket["wing"] += 1
            if "cap" in dom.lower(): bucket["cap"] += 1
            if "fusion" in dom.lower(): bucket["fusion"] += 1
            if "mr191" in dom.lower(): bucket["mr191_annotated"] += 1
        domain_count_rows.append({"sample": sample, "total_aa_diffs_vs_musoke": n_mut,
            "core_helix_strand": bucket["core_helix_strand"], "signal": bucket["signal"],
            "mucin": bucket["mucin"], "wing": bucket["wing"], "cap": bucket["cap"],
            "fusion": bucket["fusion"], "mr191_annotated": bucket["mr191_annotated"], "qc_flag": qc})
        for pos in MR191_SITES:
            if pos > len(musoke): continue
            ra, sa = musoke[pos-1], pos_map.get(pos, "X")
            cons = "gap" if sa == "-" else ("missing" if sa == "X" else ("yes" if sa == ra else "no"))
            mr191_rows.append({"sample": sample, "aa_pos": pos, "musoke_aa": ra, "sample_aa": sa,
                "conserved": cons, "domains": domain_of(domains, pos), "qc_flag": qc})
        for pos, (wt, mut, note) in ESCAPE.items():
            sa = pos_map.get(pos, "X")
            if sa in ("X", "-"): status = "missing"
            elif sa == wt: status = "wildtype"
            elif sa == mut: status = "escape_allele"
            else: status = "other:%s" % sa
            esc_rows.append({"sample": sample, "aa_pos": pos, "wildtype": wt, "escape_allele": mut,
                "sample_aa": sa, "status": status, "note": note, "domains": domain_of(domains, pos), "qc_flag": qc})
        faa_for_tree.append(SeqRecord(Seq(gp.replace("X", "-")), id=sample, description=qc))

    ids = list(all_gps.keys())
    matrix_rows = []
    for a in ids:
        row = {"id": a}
        for b in ids:
            if a == b: row[b] = "100.00"
            else:
                ident, _, _, _ = align_map(all_gps[a], all_gps[b])
                row[b] = "%.2f" % ident
        matrix_rows.append(row)

    write_tsv("gp_similarity.tsv", sim_rows)
    write_tsv("gp_mutations_vs_musoke.tsv", mut_rows)
    write_tsv("gp_mr191_sites.tsv", mr191_rows)
    write_tsv("gp_escape_sites.tsv", esc_rows)
    write_tsv("gp_domain_summary.tsv", domain_count_rows)
    write_tsv("gp_quality.tsv", quality_rows)
    write_tsv("gp_identity_matrix.tsv", matrix_rows)
    SeqIO.write(faa_for_tree, "gp_all_for_tree.faa", "fasta")

    lines = ["# GP in-depth analysis\n", "Musoke GP length: %d (%s)\n" % (len(musoke), gsrc)]
    for sample in sorted(samples):
        sim = next(r for r in sim_rows if r["sample"]==sample and r["comparator"]=="Musoke_ref")
        cons = [r for r in mr191_rows if r["sample"]==sample]
        n_ok = sum(1 for r in cons if r["conserved"]=="yes")
        esc = [r for r in esc_rows if r["sample"]==sample]
        dom = next(r for r in domain_count_rows if r["sample"]==sample)
        q = next(r for r in quality_rows if r["sample"]==sample)
        lines.append("## %s [%s]\n" % (sample, q["qc_flag"]))
        lines.append("- Identity vs Musoke: **%s%%** (X=%s)\n" % (sim["pct_identity"], sim["fraction_X"]))
        lines.append("- MR191 conserved: **%d/%d**\n" % (n_ok, len(cons)))
        lines.append("- Escape: " + ", ".join("%s=%s" % (r["aa_pos"], r["status"]) for r in esc) + "\n")
        lines.append("- Domain diffs total=%s core=%s mucin=%s signal=%s wing=%s\n" % (
            dom["total_aa_diffs_vs_musoke"], dom["core_helix_strand"], dom["mucin"], dom["signal"], dom["wing"]))
        lines.append("")
    Path("gp_summary.md").write_text("\n".join(lines))
    print("Done samples=", len(samples))

if __name__ == "__main__":
    main()
