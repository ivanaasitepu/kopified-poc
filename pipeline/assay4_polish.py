#!/usr/bin/env python3
"""Assay 4 — polishing panel v0.1:
1) kandidat REF cadangan yang juga mengamplifikasi Cara_1.0
2) pengganti marker lemah (A1_0008) dari kandidat M1-Ce terbaik berikutnya
3) cek heterodimer semua primer (calc_heterodimer, ambang dG < -9000 cal/mol = berisiko)
Menulis: results/primers_v0.1.tsv (diperbarui), results/multiplex_check.tsv
"""
import csv, itertools, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *
import primer3

OUT = f"{BASE}/output/assay1"
RES = f"{BASE}/results"
G = find_genomes()

RC = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(RC)[::-1]
def mismatches(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def insilico_pcr(seqs, fwd, rev, max_mm=2, min_amp=50, max_amp=1000):
    rrev = revcomp(rev); fa, ra = fwd[-12:], rrev[-12:]
    count = 0
    for genome_seq in seqs.values():
        i = genome_seq.find(fa)
        while i != -1:
            s = i + 12 - len(fwd)
            if s >= 0:
                sub = genome_seq[s:s+len(fwd)]
                if len(sub) == len(fwd) and mismatches(sub, fwd) <= max_mm and mismatches(sub[-3:], fwd[-3:]) == 0:
                    j = genome_seq.find(ra, s)
                    while j != -1 and j - s <= max_amp + len(rev):
                        rs = j + 12 - len(rev)
                        if rs >= 0:
                            rsub = genome_seq[rs:rs+len(rev)]
                            if len(rsub) == len(rev) and mismatches(rsub, rrev) <= max_mm and mismatches(rsub[-3:], rrev[-3:]) == 0:
                                if min_amp <= j + 12 - s <= max_amp:
                                    count += 1
                        j = genome_seq.find(ra, j + 1)
            i = genome_seq.find(fa, i + 1)
    return count

genome_cache = {}
def load_genome(f):
    if f not in genome_cache:
        genome_cache[f] = dict(iter_fasta(f"{GENOMES}/{f}"))
    return genome_cache[f]

BASE_ARGS = dict(PRIMER_NUM_RETURN=3, PRIMER_PRODUCT_SIZE_RANGE=[[70, 130]],
                 PRIMER_OPT_SIZE=20, PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=25,
                 PRIMER_OPT_TM=60.0, PRIMER_MIN_TM=57.0, PRIMER_MAX_TM=63.0,
                 PRIMER_OPT_GC_PERCENT=50.0, PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=60.0,
                 PRIMER_MAX_POLY_X=4, PRIMER_MAX_SELF_ANY_TH=45.0, PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0)

def design(mid, seq):
    for wp in (True, False):
        args = dict(BASE_ARGS)
        if wp:
            args.update(dict(PRIMER_PICK_INTERNAL_OLIGO=1, PRIMER_INTERNAL_OPT_TM=69.0,
                             PRIMER_INTERNAL_MIN_TM=63.0, PRIMER_INTERNAL_MAX_TM=72.0,
                             PRIMER_INTERNAL_MAX_POLY_X=3))
        else:
            args["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        res = primer3.bindings.design_primers({"SEQUENCE_ID": mid, "SEQUENCE_TEMPLATE": seq}, args)
        if res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0:
            return dict(marker_id=mid, fwd=res["PRIMER_LEFT_0_SEQUENCE"], rev=res["PRIMER_RIGHT_0_SEQUENCE"],
                        probe=res.get("PRIMER_INTERNAL_0_SEQUENCE", ""),
                        amplicon=res["PRIMER_PAIR_0_PRODUCT_SIZE"],
                        tm_f=round(res["PRIMER_LEFT_0_TM"], 1), tm_r=round(res["PRIMER_RIGHT_0_TM"], 1),
                        tm_p=round(res.get("PRIMER_INTERNAL_0_TM", 0), 1))
    return None

COFFEE = [f for f, (sp, grp) in sorted(G.items()) if grp in ("eugenioides", "canephora", "arabica")]
def test_coffee(mid, mtype, d, expect, forbid):
    obs = {}
    for f in COFFEE:
        grp = G[f][1]
        obs[grp] = obs.get(grp, 0) + insilico_pcr(load_genome(f), d["fwd"], d["rev"])
    obs_grp = {g for g, n in obs.items() if n > 0}
    lolos = all(e in obs_grp for e in expect) and not any(x in obs_grp for x in forbid)
    print(f"{mid:12} {mtype:6} obs={obs} lolos={lolos}")
    return lolos

panel = list(csv.DictReader(open(f"{RES}/primers_v0.1.tsv"), delimiter="\t"))
new_entries = []

# ---------- 1) REF cadangan ----------
print("=== REF cadangan ===")
n_ref_new = 0
for name, seq in iter_fasta(f"{OUT}/candidates_cc_ref.fasta"):
    if not name.startswith("REF_"):
        continue
    if n_ref_new >= 2:
        break
    mid = name.split("|")[0]
    if mid in {p["marker_id"] for p in panel}:
        continue
    gc = float(name.split("gc=")[1])
    if not 0.38 <= gc <= 0.58:
        continue
    d = design(mid, seq)
    if not d:
        continue
    if test_coffee(mid, "REF", d, {"arabica", "canephora", "eugenioides"}, set()):
        d["type"] = "REF"
        new_entries.append(d); n_ref_new += 1
print("REF baru lolos:", n_ref_new)

# ---------- 2) pengganti A1_0008 ----------
print("=== pengganti A1_0008 ===")
seqs_ce = {name.split("|")[0]: seq for name, seq in iter_fasta(f"{OUT}/assay1_candidates.fasta")}
rows = list(csv.reader(open(f"{OUT}/assay1_marker_panel_top200.csv")))
hdr = rows[0]
can_cols = [i for i, h in enumerate(hdr) if h == "C. canephora"]
ara_cols = [i for i, h in enumerate(hdr) if h == "C. arabica"]
eug_cols = [i for i, h in enumerate(hdr) if h == "C. eugenioides"]
current_ce = {p["marker_id"] for p in panel if p["type"] == "M1-Ce"}
alt = []
for row in rows[1:]:
    def maxid(cols):
        v = [float(row[i].split("%")[0]) for i in cols if row[i] != "-"]
        return max(v) if v else 0.0
    if row[0] in current_ce:
        continue
    if sum(row[i] != "-" for i in ara_cols) == len(ara_cols) and maxid(ara_cols) >= 95 and maxid(eug_cols) >= 90 and 0 < maxid(can_cols) <= 88:
        d = design(row[0], seqs_ce[row[0]])
        if not d:
            continue
        if test_coffee(row[0], "M1-Ce", d, {"arabica", "eugenioides"}, {"canephora"}):
            d["type"] = "M1-Ce"
            new_entries.append(d)
            break  # satu pengganti cukup
# buang A1_0008 dari panel
panel = [p for p in panel if p["marker_id"] != "A1_0008"]

panel += [{**d, "amplicon_bp": d["amplicon"], "tm_fwd": d["tm_f"], "tm_rev": d["tm_r"], "tm_probe": d["tm_p"]} for d in new_entries]
with open(f"{RES}/primers_v0.1.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["marker_id", "type", "fwd", "rev", "probe", "amplicon_bp", "tm_fwd", "tm_rev", "tm_probe"])
    for p in panel:
        w.writerow([p["marker_id"], p["type"], p["fwd"], p["rev"], p["probe"], p["amplicon_bp"],
                    p["tm_fwd"], p["tm_rev"], p["tm_probe"]])
print("panel sekarang:", len(panel), "pasangan")

# ---------- 3) cek heterodimer ----------
print("=== heterodimer ===")
oligos = []
for p in panel:
    oligos.append((p["marker_id"] + ":F", p["fwd"]))
    oligos.append((p["marker_id"] + ":R", p["rev"]))
    if p.get("probe"):
        oligos.append((p["marker_id"] + ":P", p["probe"]))
risky = []
for (n1, s1), (n2, s2) in itertools.combinations(oligos, 2):
    r = primer3.calc_heterodimer(s1, s2)
    if r.dg < -9000:
        risky.append((n1, n2, round(r.dg / 1000, 1)))
with open(f"{RES}/multiplex_check.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["oligo_1", "oligo_2", "dG_kcal_mol", "status"])
    for n1, n2, dg in risky:
        w.writerow([n1, n2, dg, "RISIKO (dG < -9)"])
print(f"heterodimer berisiko: {len(risky)} dari {len(oligos)*(len(oligos)-1)//2} pasangan oligo")
for r in risky[:10]:
    print(" ", r)
