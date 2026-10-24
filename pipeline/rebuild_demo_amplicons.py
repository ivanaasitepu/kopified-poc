#!/usr/bin/env python3
"""Bangun ulang fragmen amplikon sampel demo TANPA genom.

Fragmen latar (non-amp_) dan sekuens amplikon per marker diambil dari file sampel yang sudah ada;
hanya jumlah salinan amplikon yang dihitung ulang dengan model copy-number make_demo_samples2.py.
Pakai ini bila genom sumber tidak tersedia; dengan genom, jalankan make_demo_samples2.py saja.
"""
import csv, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(HERE, "..", "dashboard", "samples")
primers = list(csv.DictReader(open(os.path.join(HERE, "..", "results", "primers_v0.1.tsv")), delimiter="\t"))
K = 10  # sama dengan default amp_frags di make_demo_samples2.py

# (file, proporsi robusta, ada kapang OTA)
PLAN = [("arabica_caturra.fa", 0.0, False), ("blend_10pct.fa", 0.10, False), ("blend_25pct.fa", 0.25, False),
        ("blend_50pct.fa", 0.50, False), ("arabica_with_carbonarius.fa", 0.0, True), ("mold_free.fa", 0.0, False)]

def read_records(path):
    lines = open(path).read().split("\n")
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines) - 1, 2) if lines[i].startswith(">")]

AMP = {}
for f in os.listdir(S):
    if f.endswith(".fa"):
        for name, seq in read_records(os.path.join(S, f)):
            if name.startswith("amp_"):
                AMP.setdefault(name[4:].rsplit("_", 1)[0], seq)

def amp_frags(r_robusta, mold):
    out, n = [], 0
    for p in primers:
        mid, t = p["marker_id"], p["type"]
        if mid not in AMP:
            continue
        if t == "M1-Ce": reps = round(2 * (1 - r_robusta) * K)
        elif t == "M1-Cc": reps = round((2 * (1 - r_robusta) + 2 * r_robusta) * K)
        elif t == "REF": reps = round((4 * (1 - r_robusta) + 2 * r_robusta) * K)
        else: reps = 3 if mold else 0
        for _ in range(reps):
            out.append(f">amp_{mid}_{n}\n{AMP[mid]}")
            n += 1
    return out

for fn, r, mold in PLAN:
    path = os.path.join(S, fn)
    bg = [f">{name}\n{seq}" for name, seq in read_records(path) if not name.startswith("amp_")]
    parts = bg + amp_frags(r, mold)
    random.Random(99).shuffle(parts)
    with open(path, "w") as fh:
        fh.write("\n".join(parts) + "\n")
    print(f"{fn}: {len(parts) - len(bg)} amplikon + {len(bg)} latar")
