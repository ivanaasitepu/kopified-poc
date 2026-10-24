#!/usr/bin/env python3
"""Assay 3c — finalisasi panel primer.
1) Desain ulang marker OTA pada SEKUENS GENOMIK (rentang gen otaB/otaD di klaster MW526250,
   jendela dalam satu exon, konservatif di steynii; plus carbonarius-native).
2) In-silico PCR per-genom untuk SEMUA primer panel; tulis ulang insilico_pcr.tsv dengan benar.
"""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *
import primer3
from Bio import SeqIO

OUT = f"{BASE}/output/assay1"
OUT2 = f"{BASE}/output/assay2"
RES = f"{BASE}/results"
GENEDIR = f"{BASE}/data/ota_genes"
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
    for with_probe in (True, False):
        args = dict(BASE_ARGS)
        if with_probe:
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

# ---------- 1) marker OTA dari sekuens genomik klaster ----------
rec = SeqIO.read(f"{GENEDIR}/MW526250.gb", "genbank")
cluster = str(rec.seq)
# rentang exon otaB & otaD pada koordinat genomik
exons = {}
for feat in rec.features:
    if feat.type == "CDS" and feat.qualifiers.get("gene", [""])[0] in ("otaB", "otaD", "oatA"):
        gname = feat.qualifiers["gene"][0].replace("oatA", "otaA")
        exons[gname] = [(int(p.start), int(p.end)) for p in feat.location.parts]
print("exon:", {k: len(v) for k, v in exons.items()})

# jendela konservatif: blastn klaster vs steynii, ambil posisi klaster tercakup >=75% id
bl = subprocess.run(["blastn", "-query", f"{GENEDIR}/MW526250.fasta", "-db", f"{OUT2}/db_GCF_002849105.1_Aspste1",
                     "-task", "blastn", "-outfmt", "6 qstart qend pident", "-evalue", "1e-10", "-num_threads", "4"],
                    check=True, capture_output=True, text=True)
cov = set()
for line in bl.stdout.splitlines():
    qs, qe, p = line.split("\t")
    if float(p) >= 75:
        cov.update(range(min(int(qs), int(qe)), max(int(qs), int(qe))))

new_primers = []
for gname in ["otaB", "otaD"]:
    made = 0
    for es, ee in exons.get(gname, []):
        if made >= 1:
            break
        # jendela 300 bp di dalam exon yang tercakup steynii >=80%
        for s in range(es, ee - 300, 150):
            w = cluster[s:s+300]
            if sum(1 for i in range(s, s+300) if i in cov) < 240:
                continue
            gc = (w.count("G") + w.count("C")) / 300
            if not (0.38 <= gc <= 0.62):
                continue
            mid = f"OTA-{gname[-1].upper()}-w"
            d = design(mid, w)
            if d:
                d["type"] = f"OTA-{gname[-1].upper()}"
                new_primers.append(d)
                made += 1
                print(f"desain {mid} pada {gname} exon {es}-{ee}")
                break

# carbonarius-native untuk otaB (yang sebelumnya lolos dipertahankan dari panel lama)
existing = list(csv.DictReader(open(f"{RES}/primers_v0.1.tsv"), delimiter="\t"))
print("panel lama:", len(existing))

# ---------- 2) in-silico PCR per-genom untuk semua primer ----------
all_primers = existing + [{"marker_id": d["marker_id"], "type": d["type"], "fwd": d["fwd"], "rev": d["rev"],
                           "probe": d["probe"], "amplicon_bp": d["amplicon"], "tm_fwd": d["tm_f"],
                           "tm_rev": d["tm_r"], "tm_probe": d["tm_p"]} for d in new_primers]
# dedup marker_id
seen = set()
dedup = []
for p in all_primers:
    if p["marker_id"] not in seen:
        seen.add(p["marker_id"])
        dedup.append(p)
all_primers = dedup
# ---------- tulis ulang insilico_pcr.tsv per genom ----------
genomes_sorted = sorted(G.items())
def colkey(f, sp):
    return f"{sp}|" + f.split("_genomic")[0]
cols = ["marker_id", "type"] + [f"{sp}|{f.split(chr(95)+chr(103)+chr(101)+chr(110)+chr(111)+chr(109)+chr(105)+chr(99))[0]}" for f, (sp, grp) in genomes_sorted]
out_rows = []
for p in all_primers:
    row = {"marker_id": p["marker_id"], "type": p["type"]}
    for f, (sp, grp) in genomes_sorted:
        is_coffee_marker = p["type"] in ("M1-Ce", "M1-Cc", "REF")
        if is_coffee_marker and grp in ("mold_target",):
            row[colkey(f, sp)] = ""
            continue
        n = insilico_pcr(load_genome(f), p["fwd"], p["rev"])
        row[colkey(f, sp)] = n
    out_rows.append(row)
    obs = {k for k, v in row.items() if k not in ("marker_id", "type") and v not in ("", 0, "0")}
    print(p["marker_id"], p["type"], "->", len(obs), "genom positif")

with open(f"{RES}/insilico_pcr.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in out_rows:
        w.writerow([r.get(c, "") for c in cols])

# ---------- tulis ulang primers_v0.1.tsv (panel lama + baru) ----------
with open(f"{RES}/primers_v0.1.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["marker_id", "type", "fwd", "rev", "probe", "amplicon_bp", "tm_fwd", "tm_rev", "tm_probe"])
    for p in all_primers:
        w.writerow([p["marker_id"], p["type"], p["fwd"], p["rev"], p["probe"], p["amplicon_bp"],
                    p["tm_fwd"], p["tm_rev"], p["tm_probe"]])
print(f"panel final: {len(all_primers)} pasangan primer")
