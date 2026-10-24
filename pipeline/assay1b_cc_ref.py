#!/usr/bin/env python3
"""Assay 1b — marker M1-Cc dan REF sesuai doc Drey.
M1-Cc: wilayah genom canephora yang TIDAK tertutup eugenioides, TAPI tertutup semua arabica
       (kopi arabika membawa subgenom canephora -> bukti komponen robusta; pasangan rasio).
REF  : wilayah eugenioides yang tertutup KEDUA canephora + KEDUA arabica + eugenioides-2
       (genom inti konservatif, kontrol internal DNA kopi)."""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay1"
os.makedirs(OUT, exist_ok=True)
G = find_genomes()
by_group = {}
for f, (sp, grp) in G.items():
    by_group.setdefault(grp, []).append(f)

# ---------- M1-Cc ----------
CANREF = f"{GENOMES}/GCA_900059795.1_AUK_PRJEB4211_v1_genomic.fna.gz"  # canephora kromosom-level
print("ref canephora:", os.path.basename(CANREF))
can_len = {n: len(s) for n, s in iter_fasta(CANREF)}
cov = {}
targets = by_group["eugenioides"] + sorted(by_group["arabica"])[:2]
for f in targets:
    paf = f"{OUT}/aln_{f.split('_genomic')[0]}_vs_canref.paf"
    print("aligning", f)
    run_minimap2(CANREF, f"{GENOMES}/{f}", paf, threads=4)
    cov[f] = paf_ref_coverage(paf)

merged_eug = {}
for f in by_group["eugenioides"]:
    for name, ivs in cov[f].items():
        merged_eug[name] = merge_intervals(merged_eug.get(name, []) + ivs)
uniq_cc = subtract_intervals(can_len, merged_eug)
print(f"canephora-unique vs eugenioides: {sum(e-s for ivs in uniq_cc.values() for s,e in ivs)/1e6:.1f} Mb")
sel = uniq_cc
for f in sorted(by_group["arabica"])[:2]:
    sel = intersect_intervals(sel, cov[f], min_ov=100)
    print(f"setelah syarat {f[:20]}: {sum(e-s for ivs in sel.values() for s,e in ivs)/1e6:.2f} Mb")

cands_cc = extract_regions(CANREF, sel)
print("kandidat M1-Cc:", len(cands_cc))

# ---------- REF ----------
eugref = f"{GENOMES}/GCF_003713205.1_Ceug_1.0_genomic.fna.gz"
eug_len = {n: len(s) for n, s in iter_fasta(eugref)}
def cov_of(paf):
    return paf_ref_coverage(paf)
# pakai PAF yang sudah ada dari run utama (vs eug ref)
core = {k: [(0, eug_len[k])] for k in eug_len}
for f in sorted(by_group["canephora"]) + sorted(by_group["arabica"])[:2] + [x for x in by_group["eugenioides"] if not x.startswith("GCF_003713205")]:
    paf = f"{OUT}/aln_{f.split('_genomic')[0]}_vs_ref.paf"
    if not os.path.exists(paf):
        print("PAF belum ada (alignment utama masih berjalan?):", paf); continue
    core = intersect_intervals(core, cov_of(paf), min_ov=100)
print(f"genom inti konservatif: {sum(e-s for ivs in core.values() for s,e in ivs)/1e6:.1f} Mb")
cands_ref = extract_regions(eugref, core)
print("kandidat REF:", len(cands_ref))

# ---------- simpan kandidat gabungan ----------
with open(f"{OUT}/candidates_cc_ref.fasta", "w") as fh:
    for i, c in enumerate(cands_cc):
        fh.write(f">Cc_{i:04d}|{c['chrom']}:{c['start']}-{c['end']}|len={c['len']}|gc={c['gc']}\n{c['seq']}\n")
    for i, c in enumerate(cands_ref):
        fh.write(f">REF_{i:04d}|{c['chrom']}:{c['start']}-{c['end']}|len={c['len']}|gc={c['gc']}\n{c['seq']}\n")
print("tersimpan:", f"{OUT}/candidates_cc_ref.fasta")
