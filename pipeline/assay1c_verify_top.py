#!/usr/bin/env python3
"""Verifikasi BLAST 200 kandidat M1-Ce teratas vs 6 DB kopi (DB sudah dibangun).
Output: output/assay1/assay1_marker_panel_top200.csv"""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay1"
G = find_genomes()
COFFEE = [f for f, (sp, grp) in sorted(G.items()) if grp in ("eugenioides", "canephora", "arabica")
          and os.path.exists(f"{OUT}/db_{f.split('_genomic')[0]}.nsq")]

cands = []
for name, seq in iter_fasta(f"{OUT}/assay1_candidates.fasta"):
    mid = name.split("|")[0]
    ln = int(name.split("len=")[1].split("|")[0])
    gc = float(name.split("gc=")[1].split("|")[0])
    if 0.38 <= gc <= 0.58:
        cands.append((mid, ln, gc, seq))
cands.sort(key=lambda x: -x[1])
top = cands[:200]
fa = f"{OUT}/assay1_top200.fasta"
with open(fa, "w") as fh:
    for mid, ln, gc, seq in top:
        fh.write(f">{mid}\n{seq}\n")
print("200 kandidat teratas ditulis")

results = {}
for f in COFFEE:
    db = f"{OUT}/db_{f.split('_genomic')[0]}"
    assert os.path.exists(db + ".nsq"), f"DB hilang: {db}"
    bl = subprocess.run(["blastn", "-query", fa, "-db", db, "-outfmt", "6 qseqid pident length evalue",
                         "-evalue", "1e-20", "-max_target_seqs", "1", "-num_threads", "8", "-dust", "no"],
                        check=True, capture_output=True, text=True)
    best = {}
    for line in bl.stdout.splitlines():
        q, p, l, e = line.split("\t")
        q = q.split("|")[0]
        if q not in best or float(e) < float(best[q][2]):
            best[q] = (float(p), int(l), e)
    results[f] = best
    print(f"{f[:44]:46} hit: {len(best)}/200")

hdr = ["marker", "len", "gc"] + [G[f][0] for f in COFFEE]
with open(f"{OUT}/assay1_marker_panel_top200.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    for mid, ln, gc, seq in top:
        row = [mid, ln, gc]
        for f in COFFEE:
            b = results[f].get(mid)
            row.append(f"{b[0]:.1f}%/{b[1]}bp" if b else "-")
        w.writerow(row)

n_perfect = 0
for mid, ln, gc, seq in top:
    ara = sum(1 for f in COFFEE if G[f][1] == "arabica" and mid in results[f])
    can = sum(1 for f in COFFEE if G[f][1] == "canephora" and mid in results[f])
    eug = sum(1 for f in COFFEE if G[f][1] == "eugenioides" and mid in results[f])
    if ara == 2 and can == 0 and eug == 2:
        n_perfect += 1
print(f"M1-Ce sempurna (2/2 arabica, 0 canephora, 2/2 eugenioides): {n_perfect}/200")
print("panel tersimpan:", f"{OUT}/assay1_marker_panel_top200.csv")
