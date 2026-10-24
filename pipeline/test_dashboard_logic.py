#!/usr/bin/env python3
"""Uji logika dashboard pada 7 sampel demo (CLI, sebelum demo live)."""
import csv, sys, time
sys.path.insert(0, "/work/ivanasitepu/genohack_Kopicert/build/dashboard")
from insilico import count_amplicons

BASE = "/work/ivanasitepu/genohack_Kopicert/build"
primers = list(csv.DictReader(open(f"{BASE}/results/primers_v0.1.tsv"), delimiter="\t"))
sim = list(csv.DictReader(open(f"{BASE}/results/blend_simulation.tsv"), delimiter="\t"))
r0 = [r for r in sim if r["robusta_pct"] == "0"][0]
r25 = [r for r in sim if r["robusta_pct"] == "25"][0]
THRESH = float(r0["ratio_Cc_Ce"]) + (float(r25["ratio_Cc_Ce"]) - float(r0["ratio_Cc_Ce"])) * 0.3
RATIO_EXCLUDE = {"A1_0008"}  # sama dengan dashboard/app.py
print(f"ambang rasio Cc/Ce: {THRESH:.2f}\n")

def read_fasta(path):
    seqs, buf = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if buf: seqs.append("".join(buf))
            buf = []
        elif line: buf.append(line)
    if buf: seqs.append("".join(buf))
    return seqs

import os
for fn in sorted(os.listdir(f"{BASE}/dashboard/samples")):
    if not fn.endswith(".fa"): continue
    t0 = time.time()
    seqs = read_fasta(f"{BASE}/dashboard/samples/{fn}")
    ce = sum(count_amplicons(seqs, p["fwd"], p["rev"]) for p in primers if p["type"] == "M1-Ce" and p["marker_id"] not in RATIO_EXCLUDE)
    cc = sum(count_amplicons(seqs, p["fwd"], p["rev"]) for p in primers if p["type"] == "M1-Cc")
    ref = sum(count_amplicons(seqs, p["fwd"], p["rev"]) for p in primers if p["type"] == "REF")
    ota = {p["marker_id"]: count_amplicons(seqs, p["fwd"], p["rev"]) for p in primers if p["type"].startswith("OTA")}
    n_ce = max(1, sum(1 for q in primers if q["type"]=="M1-Ce" and q["marker_id"] not in RATIO_EXCLUDE)); n_cc = max(1, sum(1 for q in primers if q["type"]=="M1-Cc"))
    ratio = (cc / n_cc) / (ce / n_ce) if ce else float("nan")
    if ref == 0 and ce == 0 and cc == 0: call1 = "INVALID"
    elif ce == 0: call1 = "FLAG bukan arabika"
    elif ratio == ratio and ratio > THRESH: call1 = f"FLAG campuran robusta (rasio {ratio:.2f})"
    else: call1 = "ARABICA AUTHENTIC"
    n_ota = sum(1 for n in ota.values() if n > 0)
    call2 = "FLAG OTA" if n_ota >= 2 else ("CHECK" if n_ota == 1 else "bersih")
    print(f"{fn:34} Ce={ce:2} Cc={cc:2} REF={ref:2} rasio={ratio:5.2f} | {call1:32} | OTA: {call2}  ({time.time()-t0:.1f}s)")
