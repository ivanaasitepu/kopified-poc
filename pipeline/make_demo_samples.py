#!/usr/bin/env python3
"""Buat sampel FASTA demo dari genom asli (fragmen 150 bp, seed=42)."""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

S = f"{BASE}/dashboard/samples"
os.makedirs(S, exist_ok=True)
random.seed(42)

def load(prefix):
    for f in os.listdir(GENOMES):
        if f.startswith(prefix) and f.endswith(".fna.gz"):
            return dict(iter_fasta(f"{GENOMES}/{f}"))
    raise FileNotFoundError(prefix)

ARA = load("GCF_003713225")   # arabica Cara_1.0
CAN = load("GCA_900059795")   # canephora AUK
CARB = load("GCA_001990825")  # A. carbonarius

def frags(seqs, n, tag, seed):
    rng = random.Random(seed)
    chrs = sorted(seqs, key=lambda k: -len(seqs[k]))[:11]
    out = []
    for i in range(n):
        c = rng.choice(chrs)
        L = len(seqs[c])
        s = rng.randrange(0, L - 150)
        out.append(f">{tag}_{i}|{c}:{s}-{s+150}\n{seqs[c][s:s+150]}")
    return out

def write(name, parts):
    with open(f"{S}/{name}", "w") as fh:
        fh.write("\n".join(parts) + "\n")
    print(name, "->", sum(p.count(">") for p in parts), "fragmen,", os.path.getsize(f"{S}/{name}")//1024, "KB")

write("arabica_caturra.fa", frags(ARA, 4000, "ara", 1))
write("robusta.fa", frags(CAN, 4000, "rob", 2))
write("blend_10pct.fa", frags(ARA, 3600, "ara", 3) + frags(CAN, 400, "rob", 3))
write("blend_25pct.fa", frags(ARA, 3000, "ara", 4) + frags(CAN, 1000, "rob", 4))
write("arabica_with_carbonarius.fa", frags(ARA, 3800, "ara", 5) + frags(CARB, 200, "mold", 5))
write("mold_free.fa", frags(ARA, 4000, "ara", 6))
write("mold_ota_positive.fa", frags(CARB, 4000, "mold", 7))
