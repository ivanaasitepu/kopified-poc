#!/usr/bin/env python3
"""Assay 2 — OTA-mold early-warning markers.
Per target species S: regions of the representative assembly that are
  (a) NOT covered by any other fungal genome (other targets + A. niger),
  (b) conserved across ALL assemblies of S.
BLASTN-verify vs all fungal genomes + coffee (C. canephora as plant background).
"""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay2"
os.makedirs(OUT, exist_ok=True)
G = find_genomes()

by_species = {}
coffee = None
for f, (sp, grp) in G.items():
    if grp == "mold_target":
        by_species.setdefault(sp, []).append(f)
    elif grp == "mold_nontarget":
        by_species.setdefault(sp, [])
    elif grp == "canephora":
        coffee = f
print({k: len(v) for k, v in by_species.items() if v}, "| coffee bg:", coffee)

PREF = {  # representative (reference) assembly per species
    "A. ochraceus": "GCA_004849945",
    "A. westerdijkiae": "GCA_001307345",
    "A. carbonarius": "GCA_001990825",
    "A. steynii": "GCF_002849105",
}
all_fungi = [f for f, (sp, grp) in G.items() if grp in ("mold_target", "mold_nontarget")]
all_panels = []

for sp, files in by_species.items():
    if not files:
        continue
    ref_file = [f for f in files if f.startswith(PREF[sp])][0]
    REF = f"{GENOMES}/{ref_file}"
    others = [f for f in all_fungi if f not in files]
    ref_len = {n: len(s) for n, s in iter_fasta(REF)}
    print(f"\n=== {sp}: ref={ref_file[:30]}, size={sum(ref_len.values())/1e6:.1f} Mb, assemblies={len(files)}")

    # non-target coverage (other species)
    merged_nt = {}
    for f in others:
        paf = f"{OUT}/aln_{f.split('_genomic')[0]}_vs_{PREF[sp]}.paf"
        run_minimap2(REF, f"{GENOMES}/{f}", paf, threads=8)
        for name, ivs in paf_ref_coverage(paf).items():
            merged_nt[name] = merge_intervals(merged_nt.get(name, []) + ivs)
    uniq = subtract_intervals(ref_len, merged_nt)
    print(f"unique vs other fungi: {sum(e-s for ivs in uniq.values() for s,e in ivs)/1e6:.2f} Mb")

    # conserved across own assemblies
    sel = uniq
    for f in sorted(files):
        if f == ref_file:
            continue
        paf = f"{OUT}/aln_{f.split('_genomic')[0]}_vs_{PREF[sp]}.paf"
        run_minimap2(REF, f"{GENOMES}/{f}", paf, threads=8)
        sel = intersect_intervals(sel, paf_ref_coverage(paf), min_ov=100)
        print(f"  + own assembly {f[:24]}: {sum(e-s for ivs in sel.values() for s,e in ivs)/1e6:.2f} Mb")

    cands = extract_regions(REF, sel, min_len=400, max_len=2000)
    print(f"  candidates: {len(cands)}")

    candfa = f"{OUT}/cand_{PREF[sp]}.fasta"
    with open(candfa, "w") as fh:
        for i, c in enumerate(cands):
            fh.write(f">M_{PREF[sp]}_{i:04d}|{c['chrom']}:{c['start']}-{c['end']}|len={c['len']}\n{c['seq']}\n")

    # BLASTN vs all fungi + coffee background
    dbs = all_fungi + ([coffee] if coffee else [])
    panel = {}
    for f in dbs:
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
            q, p, l, e = line.split("\t")
            if q not in best or float(e) < float(best[q][2]):
                best[q] = (float(p), int(l), e)
        panel[f] = best

    for i, c in enumerate(cands):
        qid = f"M_{PREF[sp]}_{i:04d}"
        row = {"marker": qid, "species": sp, "len": c["len"], "gc": c["gc"]}
        for f in dbs:
            b = panel[f].get(qid)
            row[f] = f"{b[0]:.1f}%/{b[1]}bp" if b else "-"
        all_panels.append(row)

hdr = ["marker", "species", "len", "gc"] + [f for f in (all_fungi + ([coffee] if coffee else []))]
with open(f"{OUT}/assay2_marker_panel.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    for r in all_panels:
        w.writerow([r.get(h, "-") for h in hdr])
print("\npanel saved:", f"{OUT}/assay2_marker_panel.csv", "rows:", len(all_panels))
