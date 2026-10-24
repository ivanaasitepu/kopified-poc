#!/usr/bin/env python3
"""Shared helpers for KopiFied marker discovery."""
import gzip, os, subprocess

BASE = "/work/ivanasitepu/genohack_Kopicert/build"
GENOMES = f"{BASE}/data/genomes"

MAP = {
    "GCF_003713205": ("C. eugenioides", "eugenioides"),
    "GCA_036785875": ("C. eugenioides", "eugenioides"),
    "GCA_900059795": ("C. canephora", "canephora"),
    "GCA_036785865": ("C. canephora", "canephora"),
    "GCF_003713225": ("C. arabica", "arabica"),
    "GCF_036785885": ("C. arabica", "arabica"),
    "GCA_030873655": ("C. arabica", "arabica"),
    "GCA_001990825": ("A. carbonarius", "mold_target"),
    "GCF_002849105": ("A. steynii", "mold_target"),
    "GCA_001307345": ("A. westerdijkiae", "mold_target"),
    "GCA_026259865": ("A. westerdijkiae", "mold_target"),
    "GCA_019976475": ("A. westerdijkiae", "mold_target"),
    "GCA_004849945": ("A. ochraceus", "mold_target"),
    "GCA_025594125": ("A. ochraceus", "mold_target"),
    "GCA_005784425": ("A. ochraceus", "mold_target"),
    "GCA_026108165": ("A. ochraceus", "mold_target"),
    "GCA_037127345": ("A. ochraceus", "mold_target"),
    "GCA_051529875": ("A. niger", "mold_nontarget"),
    "GCA_047651775": ("A. niger", "mold_nontarget"),
    "GCF_009017415": ("A. flavus", "mold_nontarget"),
    "GCA_000184455": ("A. oryzae", "mold_nontarget"),
    "GCA_988135455": ("A. tubingensis", "mold_nontarget"),
    "GCF_028828655": ("P. verrucosum", "mold_nontarget"),
    "GCA_051861825": ("S. cerevisiae", "microflora"),
    "GCF_061378395": ("L. plantarum", "microflora"),
}

def find_genomes():
    out = {}
    for f in sorted(os.listdir(GENOMES)):
        if f.endswith(".fna.gz") and os.path.exists(os.path.join(GENOMES, f + ".DONE")):
            parts = f.split("_genomic")[0].split("_")
            key = (parts[0] + "_" + parts[1]).rsplit(".", 1)[0]
            if key in MAP:
                out[f] = MAP[key]
    return out

def iter_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf)

def run_minimap2(ref_fasta, qry_fasta, paf_out, threads=4):
    if os.path.exists(paf_out) and os.path.getsize(paf_out) > 0:
        return paf_out
    subprocess.run(["minimap2", "-x", "asm5", "-f", "0.02", "-I", "2000M", "-t", str(threads), "-o", paf_out, ref_fasta, qry_fasta], check=True)
    return paf_out

def paf_ref_coverage(paf_path, min_mapq=20, min_alen=100):
    cov = {}
    with open(paf_path) as fh:
        for line in fh:
            f = line.rstrip().split("\t")
            if len(f) < 12 or line.startswith("-"):
                continue
            if int(f[11]) < min_mapq or int(f[10]) < min_alen:
                continue
            cov.setdefault(f[5], []).append((int(f[7]), int(f[8])))
    return {k: merge_intervals(v) for k, v in cov.items()}

def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]

def subtract_intervals(ref_len_map, covered):
    out = {}
    for name, L in ref_len_map.items():
        prev = 0
        for s, e in covered.get(name, []):
            if s > prev:
                out.setdefault(name, []).append((prev, s))
            prev = max(prev, e)
        if prev < L:
            out.setdefault(name, []).append((prev, L))
    return out

def intersect_intervals(a_map, b_map, min_ov=50):
    out = {}
    for name, ivs in a_map.items():
        for s, e in ivs:
            for bs, be in b_map.get(name, []):
                os_, oe = max(s, bs), min(e, be)
                if oe - os_ >= min_ov:
                    out.setdefault(name, []).append((os_, oe))
                    break
    return {k: merge_intervals(v) for k, v in out.items()}

def extract_regions(fasta_gz, regions, min_len=400, max_len=2500, gc_min=0.25, gc_max=0.65, max_n_frac=0.05):
    """regions: name->[(s,e)]; returns list of dicts with chunked sequences"""
    cands = []
    reg = {k: v for k, v in regions.items()}
    for name, seq in iter_fasta(fasta_gz):
        if name not in reg:
            continue
        for s, e in reg[name]:
            seg = seq[s:e]
            # split into windows of max_len
            i = 0
            while i < len(seg):
                w = seg[i:i+max_len]
                i += max_len
                if len(w) < min_len:
                    continue
                gc = (w.upper().count("G") + w.upper().count("C")) / len(w)
                nf = w.upper().count("N") / len(w)
                if not (gc_min <= gc <= gc_max) or nf > max_n_frac:
                    continue
                cands.append({"chrom": name, "start": s + i - max_len if i >= max_len else s,
                              "end": min(s + i, e), "len": len(w), "gc": round(gc, 3), "seq": w})
    return cands
