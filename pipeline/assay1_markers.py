#!/usr/bin/env python3
"""Assay 1 — Arabica authenticity markers.
Logic: regions of the C. eugenioides genome that are
  (a) NOT covered by any C. canephora assembly (absent from pure Robusta),
  (b) covered by every available C. arabica assembly (present in Arabica),
  (c) covered by the second C. eugenioides assembly (conserved within species).
Then BLASTN-verify each candidate locus against all coffee genomes.
"""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay1"
os.makedirs(OUT, exist_ok=True)
G = find_genomes()
by_group = {}
for f, (sp, grp) in G.items():
    by_group.setdefault(grp, []).append(f)
print("genomes:", {k: len(v) for k, v in by_group.items()})

need = ["eugenioides", "canephora", "arabica"]
for n in need:
    assert by_group.get(n), f"missing group {n}"

REF = f"{GENOMES}/{sorted(by_group['eugenioides'])[0]}"   # Ceug_1.0 preferred (GCF sorts after GCA; sort puts GCA first!) -> choose explicitly
eug_sorted = sorted(by_group["eugenioides"])
REF = f"{GENOMES}/{[f for f in eug_sorted if f.startswith('GCF_003713205')][0]}"
print("reference:", REF)

ref_len = {n: len(s) for n, s in iter_fasta(REF)}
print("reference size:", sum(ref_len.values())/1e6, "Mb")

# coverage maps on the eugenioides reference
cov = {}
for grp in ["canephora", "arabica", "eugenioides"]:
    for f in sorted(by_group[grp]):
        if f == os.path.basename(REF):
            continue
        paf = f"{OUT}/aln_{f.split('_genomic')[0]}_vs_ref.paf"
        print("aligning", f)
        run_minimap2(REF, f"{GENOMES}/{f}", paf, threads=8)
        cov[f] = paf_ref_coverage(paf)

# unique vs canephora: NOT covered by ANY canephora assembly
merged_can = {}
for f in by_group["canephora"]:
    for name, ivs in cov[f].items():
        merged_can[name] = merge_intervals(merged_can.get(name, []) + ivs)
uniq = subtract_intervals(ref_len, merged_can)
n_uniq = sum(e-s for ivs in uniq.values() for s, e in ivs)
print(f"eugenioides-unique vs canephora: {n_uniq/1e6:.1f} Mb")

# must be present in ALL arabica assemblies
sel = uniq
for f in sorted(by_group["arabica"]):
    sel = intersect_intervals(sel, cov[f], min_ov=100)
    print(f"after requiring {f[:22]}: {sum(e-s for ivs in sel.values() for s,e in ivs)/1e6:.2f} Mb")

# conserved in second eugenioides assembly
for f in sorted(by_group["eugenioides"]):
    if f == os.path.basename(REF):
        continue
    sel = intersect_intervals(sel, cov[f], min_ov=100)
    print(f"after requiring {f[:22]}: {sum(e-s for ivs in sel.values() for s,e in ivs)/1e6:.2f} Mb")

cands = extract_regions(REF, sel)
print("candidate loci:", len(cands))

# BLASTN verification vs every coffee genome
candfa = f"{OUT}/assay1_candidates.fasta"
with open(candfa, "w") as fh:
    for i, c in enumerate(cands):
        fh.write(f">A1_{i:04d}|{c['chrom']}:{c['start']}-{c['end']}|len={c['len']}|gc={c['gc']}\n{c['seq']}\n")

rows = []
COFFEE = [f for f, (sp, grp) in sorted(G.items()) if grp in ("eugenioides", "canephora", "arabica")]
for f in COFFEE:
    sp, grp = G[f]
    db = f"{OUT}/db_{f.split('_genomic')[0]}"
    if not os.path.exists(db + ".nsq"):
        subprocess.run(["makeblastdb", "-in", f"{GENOMES}/{f}", "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
    bl = subprocess.run(["blastn", "-query", candfa, "-db", db, "-outfmt",
                         "6 qseqid pident length evalue", "-evalue", "1e-20", "-max_target_seqs", "1",
                         "-num_threads", "8", "-dust", "no"],
                        check=True, capture_output=True, text=True)
    best = {}
    for line in bl.stdout.splitlines():
        q, p, l, e = line.split("\t"); q = q.split("|")[0]
        if q not in best or float(e) < float(best[q][2]):
            best[q] = (float(p), int(l), e)
    rows.append((f, sp, grp, best))

hdr = ["marker", "len", "gc"] + [f"{G[f][0]}" for f in COFFEE]
with open(f"{OUT}/assay1_marker_panel.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    for i, c in enumerate(cands):
        qid = f"A1_{i:04d}"
        row = [qid, c["len"], c["gc"]]
        for f, sp, grp, best in rows:
            if qid in best:
                p, l, e = best[qid]
                row.append(f"{p:.1f}%/{l}bp")
            else:
                row.append("-")
        w.writerow(row)
print("panel saved:", f"{OUT}/assay1_marker_panel.csv")
