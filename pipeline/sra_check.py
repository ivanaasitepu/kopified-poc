#!/usr/bin/env python3
"""Cek situs primer pada reads resekuensing arabika (bukti konservasi pada sampel independen).
Untuk tiap primer: hitung reads yang mengandung 3' 12-mer (exact, dua arah) dan primer penuh (<=1 mm)."""
import csv, os, sys

BASE = "/work/ivanasitepu/genohack_Kopicert/build"
RC = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(RC)[::-1]

primers = list(csv.DictReader(open(f"{BASE}/results/primers_v0.1.tsv"), delimiter="\t"))

def scan_fastq(path):
    """baca semua urutan read ke satu set counter — versi hemat memori: hitung per kmer"""
    anchors = {}
    for p in primers:
        for tag, seq in [("fwd", p["fwd"]), ("rev", p["rev"])]:
            a12 = seq[-12:]
            anchors[(p["marker_id"], tag)] = (a12, revcomp(a12), seq, revcomp(seq))
    counts = {k: [0, 0] for k in anchors}  # [kmer12 hit, full primer hit]
    n = 0
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h: break
            seq = fh.readline().strip().upper()
            fh.readline(); fh.readline()
            n += 1
            for k, (a, arc, full, frc) in anchors.items():
                if a in seq or arc in seq:
                    counts[k][0] += 1
                    # cek primer penuh dengan <=1 mismatch: cari anchor dulu
                    for fullseq in (full, frc):
                        idx = seq.find(fullseq[:12])
                        if idx != -1 and len(seq) >= idx + len(fullseq):
                            sub = seq[idx:idx+len(fullseq)]
                            if sum(1 for x, y in zip(sub, fullseq) if x != y) <= 1:
                                counts[k][1] += 1
                                break
    return counts, n

RUNS = ["SRR40981778", "SRR40981779"]
all_res = {}
for run in RUNS:
    fq = f"{BASE}/data/sra/{run}_1.fastq"
    if not os.path.exists(fq):
        print(run, "belum ada, lewati"); continue
    c, n = scan_fastq(fq)
    all_res[run] = (c, n)
    print(f"{run}: {n:,} reads dipindai")

# ringkasan per marker: primer dianggap "terbukti" jika 12-mer fwd DAN rev terdeteksi >=1x
hdr = ["marker", "type"] + [f"{r}_{m}" for r in all_res for m in ("fwd12", "rev12", "fwd_full", "rev_full")]
lines = ["\t".join(hdr)]
ok_count = {}
for p in primers:
    row = [p["marker_id"], p["type"]]
    for run, (c, n) in all_res.items():
        row += [str(c[(p["marker_id"], "fwd")][0]), str(c[(p["marker_id"], "rev")][0]),
                str(c[(p["marker_id"], "fwd")][1]), str(c[(p["marker_id"], "rev")][1])]
        ok = c[(p["marker_id"], "fwd")][0] > 0 and c[(p["marker_id"], "rev")][0] > 0
        ok_count.setdefault(run, [0, 0])
        ok_count[run][0] += ok; ok_count[run][1] += 1
    lines.append("\t".join(row))
open(f"{BASE}/output/sra_check/sra_primer_sites.tsv", "w").write("\n".join(lines) + "\n")
for run, (ok, tot) in ok_count.items():
    print(f"{run}: {ok}/{tot} marker dengan kedua situs primer terdeteksi di reads")
print("tersimpan:", f"{BASE}/output/sra_check/sra_primer_sites.tsv")
