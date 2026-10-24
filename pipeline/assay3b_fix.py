#!/usr/bin/env python3
"""Assay 3b — perbaikan: M1-Cc (probe opsional) + marker OTA (jendela konservatif + carbonarius-native).
Menambah baris ke results/primers_v0.1.tsv dan results/insilico_pcr.tsv."""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *
import primer3

OUT = f"{BASE}/output/assay1"
OUT2 = f"{BASE}/output/assay2"
RES = f"{BASE}/results"
G = find_genomes()

RC = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(RC)[::-1]
def mismatches(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def insilico_pcr(seqs, fwd, rev, max_mm=2, min_amp=50, max_amp=1000):
    rrev = revcomp(rev)
    fa, ra = fwd[-12:], rrev[-12:]
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
PROBE_ARGS = dict(PRIMER_PICK_INTERNAL_OLIGO=1, PRIMER_INTERNAL_OPT_TM=69.0,
                  PRIMER_INTERNAL_MIN_TM=63.0, PRIMER_INTERNAL_MAX_TM=72.0, PRIMER_INTERNAL_MAX_POLY_X=3)

def design(mid, seq):
    for with_probe in (True, False):
        args = dict(BASE_ARGS)
        if with_probe:
            args.update(PROBE_ARGS)
        else:
            args["PRIMER_PICK_INTERNAL_OLIGO"] = 0
        res = primer3.bindings.design_primers({"SEQUENCE_ID": mid, "SEQUENCE_TEMPLATE": seq}, args)
        if res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0:
            i = 0
            return dict(marker_id=mid, fwd=res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                        rev=res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                        probe=res.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", ""),
                        amplicon=res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                        tm_f=round(res[f"PRIMER_LEFT_{i}_TM"], 1), tm_r=round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
                        tm_p=round(res.get(f"PRIMER_INTERNAL_{i}_TM", 0), 1))
    return None

new_rows = []      # untuk primers_v0.1
ins_rows = []      # untuk insilico_pcr

def test_and_add(mid, mtype, fwd, rev, probe, amp, tmf, tmr, tmp, expect, forbid, groups_test):
    obs = {}
    for f, (sp, grp) in sorted(G.items()):
        if grp not in groups_test:
            continue
        obs[grp] = obs.get(grp, 0) + insilico_pcr(load_genome(f), fwd, rev)
    obs_grp = {g for g, n in obs.items() if n > 0}
    lolos = all(e in obs_grp for e in expect) and not any(f in obs_grp for f in forbid)
    print(f"{mid:12} {mtype:8} obs={obs_grp} expect={expect} lolos={lolos}")
    if lolos:
        new_rows.append([mid, mtype, fwd, rev, probe, amp, tmf, tmr, tmp])
        r = {"marker_id": mid, "type": mtype}
        for f, (sp, grp) in sorted(G.items()):
            if grp in groups_test:
                r[f"{sp}|{grp}"] = obs.get(grp, 0) if False else ""
        ins_rows.append((mid, mtype, obs))
    return lolos

# ---------- M1-Cc ----------
print("=== M1-Cc ===")
cc = []
for name, seq in iter_fasta(f"{OUT}/candidates_cc_ref.fasta"):
    if name.startswith("Cc_") and float(name.split("gc=")[1]) >= 0.38:
        cc.append((name.split("|")[0], seq, int(name.split("len=")[1].split("|")[0])))
cc.sort(key=lambda x: -x[2])
n_cc = 0
for mid, seq, ln in cc[:40]:
    if n_cc >= 3:
        break
    d = design(mid, seq)
    if not d:
        continue
    if test_and_add(mid, "M1-Cc", d["fwd"], d["rev"], d["probe"], d["amplicon"], d["tm_f"], d["tm_r"], d["tm_p"],
                    expect={"arabica", "canephora"}, forbid={"eugenioides", "mold_target", "mold_nontarget"},
                    groups_test={"arabica", "canephora", "eugenioides", "mold_nontarget", "microflora"}):
        n_cc += 1
print("M1-Cc lolos:", n_cc)

# ---------- OTA: jendela konservatif + carbonarius-native ----------
print("=== OTA ===")
GENEDIR = f"{BASE}/data/ota_genes"
def extract_homolog(gene_fa, db, min_id=65):
    """kembalikan sekuens region homolog terbaik dari genom di db"""
    bl = subprocess.run(["blastn", "-query", gene_fa, "-db", db, "-task", "blastn",
                         "-outfmt", "6 sseqid sstart send pident length", "-evalue", "1e-10",
                         "-max_target_seqs", "1", "-max_hsps", "1", "-num_threads", "4"],
                        check=True, capture_output=True, text=True)
    best = None
    for line in bl.stdout.splitlines():
        c, s, e, p, l = line.split("\t")
        if best is None or int(l) > best[1]:
            best = ((c, min(int(s), int(e)), max(int(s), int(e))), int(l), float(p))
    return best

producer_clades = {
    "circumdati": ["GCA_001307345.1_ASM130734v1", "GCF_002849105.1_Aspste1"],   # westerdijkiae + steynii
    "carbonarius": ["GCA_001990825.1_Aspca3"],
}
n_ota = 0
for gname in ["cluster_otaB", "cluster_otaD"]:
    fp = f"{GENEDIR}/{gname}.fasta"
    wseq = next(iter_fasta(fp))[1]
    # A) jendela di westerdijkiae yang konservatif di steynii (identitas >=75%)
    bl = subprocess.run(["blastn", "-query", fp, "-db", f"{OUT2}/db_GCF_002849105.1_Aspste1", "-task", "blastn",
                         "-outfmt", "6 qstart qend pident", "-evalue", "1e-10", "-num_threads", "4"],
                        check=True, capture_output=True, text=True)
    cov = {}
    for line in bl.stdout.splitlines():
        qs, qe, p = line.split("\t")
        if float(p) >= 75:
            for i in range(min(int(qs), int(qe)) - 1, max(int(qs), int(qe))):
                cov[i] = 1
    pos = sorted(cov)
    win = 300
    made = 0
    i = 0
    while i + win <= len(pos) and made < 2 and n_ota < 6:
        s = pos[i]
        j = i
        while j < len(pos) and pos[j] - s < win:
            j += 1
        if j - i >= win * 0.8:
            w = wseq[s:s+win]
            gc = (w.count("G") + w.count("C")) / len(w)
            if 0.38 <= gc <= 0.62:
                d = design(f"OTA-{gname[-1].upper()}-{made+1}", w)
                if d and test_and_add(d["marker_id"], f"OTA-{gname[-1].upper()}", d["fwd"], d["rev"], d["probe"],
                                      d["amplicon"], d["tm_f"], d["tm_r"], d["tm_p"],
                                      expect={"mold_target"}, forbid={"mold_nontarget", "microflora", "arabica", "canephora"},
                                      groups_test={"mold_target", "mold_nontarget", "microflora", "arabica", "canephora"}):
                    made += 1; n_ota += 1
        i = j
    # B) carbonarius-native: ekstrak region homolog, desain dari sekuens carbonarius
    db_c = f"{OUT2}/db_GCA_001990825.1_Aspca3"
    best = extract_homolog(fp, db_c)
    if best:
        (c, s, e), l, p = best
        gseqs = load_genome("GCA_001990825.1_Aspca3_genomic.fna.gz")
        region = gseqs[c][max(0, s-100):e+100]
        if p < 90:  # perlu varian primer sendiri
            mid = f"OTA-{gname[-1].upper()}-carb"
            d = design(mid, region)
            if d and test_and_add(mid, f"OTA-{gname[-1].upper()}", d["fwd"], d["rev"], d["probe"],
                                  d["amplicon"], d["tm_f"], d["tm_r"], d["tm_p"],
                                  expect={"mold_target"}, forbid={"mold_nontarget", "microflora", "arabica", "canephora"},
                                  groups_test={"mold_target", "mold_nontarget", "microflora", "arabica", "canephora"}):
                n_ota += 1

# ---------- append ke file hasil ----------
with open(f"{RES}/primers_v0.1.tsv", "a", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    for r in new_rows:
        w.writerow(r)

# perbarui insilico_pcr.tsv: tambah baris baru (kolom diselaraskan saat dibaca)
existing = list(csv.reader(open(f"{RES}/insilico_pcr.tsv"), delimiter="\t"))
cols = existing[0]
with open(f"{RES}/insilico_pcr.tsv", "a", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    for mid, mtype, obs in ins_rows:
        r = {"marker_id": mid, "type": mtype}
        for f, (sp, grp) in sorted(G.items()):
            key = f"{sp}|{grp}"
            if key in cols:
                r[key] = obs.get(grp, 0)
        w.writerow([r.get(c, "") for c in cols])
print(f"ditambahkan: {len(new_rows)} primer baru. Total di panel: {len(new_rows)} baru + yang lama")
