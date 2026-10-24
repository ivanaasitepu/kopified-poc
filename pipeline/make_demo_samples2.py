#!/usr/bin/env python3
"""Sampel demo v2: fragmen latar 150 bp + fragmen amplikon marker nyata.
Copy number mengikuti model: arabika Ce=2,Cc=2,REF=4 per sel; robusta Cc=2,REF=2; mold amplicon apa adanya."""
import csv, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

S = f"{BASE}/dashboard/samples"
RES = f"{BASE}/results"
primers = list(csv.DictReader(open(f"{RES}/primers_v0.1.tsv"), delimiter="\t"))

RC = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(RC)[::-1]

def find_amplicons(seqs, fwd, rev, max_mm=2, max_amp=1000):
    """kembalikan sekuens amplikon (fwd..rev) dari dict name->seq"""
    rrev = revcomp(rev); fa, ra = fwd[-12:], rrev[-12:]
    out = []
    for genome_seq in seqs.values():
        i = genome_seq.find(fa)
        while i != -1:
            s = i + 12 - len(fwd)
            if s >= 0 and genome_seq[s:s+len(fwd)]:
                j = genome_seq.find(ra, s)
                while j != -1 and j - s <= max_amp + len(rev):
                    rs = j + 12 - len(rev)
                    if rs >= 0:
                        amp = genome_seq[s:j+12]
                        if 50 <= len(amp) <= max_amp:
                            out.append(amp)
                    j = genome_seq.find(ra, j + 1)
            i = genome_seq.find(fa, i + 1)
    return out

def load(prefix):
    for f in os.listdir(GENOMES):
        if f.startswith(prefix) and f.endswith(".fna.gz"):
            return dict(iter_fasta(f"{GENOMES}/{f}"))

ARA = load("GCF_003713225"); CAN = load("GCA_900059795"); CARB = load("GCA_001990825"); WST = load("GCA_001307345"); ET39 = load("GCF_036785885")

# ekstrak amplikon per marker dari genom sumber yang relevan
AMP = {}
for p in primers:
    t, mid = p["type"], p["marker_id"]
    if t == "M1-Ce": src = ET39 if mid == "A1_0008" else ARA
    elif t == "M1-Cc": src = CAN
    elif t == "REF": src = CAN
    elif mid.endswith("-carb"): src = CARB
    else: src = WST
    amps = find_amplicons(src, p["fwd"], p["rev"])
    if amps:
        AMP[mid] = amps[0]
print("amplikon terekstrak:", {k: len(v) for k, v in AMP.items()})

def rand_frags(seqs, n, tag, seed):
    rng = random.Random(seed)
    chrs = sorted(seqs, key=lambda k: -len(seqs[k]))[:11]
    out = []
    for i in range(n):
        c = rng.choice(chrs)
        L = len(seqs[c])
        s = rng.randrange(0, L - 150)
        out.append(f">{tag}_{i}\n{seqs[c][s:s+150]}")
    return out

def amp_frags(r_robusta=0.0, k=10, mold=None, skip_coffee=False):
    """fragmen amplikon dgn copy number: Ce~2(1-r)k, Cc~2(1-r)k+2rk, REF~4(1-r)k+2rk; mold: tiap amplikon 3x.
    k=10 supaya 2(1-r)k bulat untuk r = 0,10/0,25/0,50 (k=5 → 7,5 dibulatkan ke 8 → blend 25% terbaca ≈80%)."""
    out = []
    n = 0
    for p in primers:
        mid, t = p["marker_id"], p["type"]
        if mid not in AMP:
            continue
        if t == "M1-Ce": reps = max(0, round(2 * (1 - r_robusta) * k))
        elif t == "M1-Cc": reps = round((2 * (1 - r_robusta) + 2 * r_robusta) * k)
        elif t == "REF": reps = round((4 * (1 - r_robusta) + 2 * r_robusta) * k)
        else: reps = 3 if mold else 0
        if skip_coffee and t != "OTA-B" and not t.startswith("OTA"):
            continue
        for _ in range(reps):
            out.append(f">amp_{mid}_{n}\n{AMP[mid]}")
            n += 1
    return out

def write(name, parts):
    random.Random(99).shuffle(parts)
    with open(f"{S}/{name}", "w") as fh:
        fh.write("\n".join(parts) + "\n")
    print(f"{name}: {sum(1 for p in parts if p.startswith('>amp_'))} amplikon + latar")

write("arabica_caturra.fa", rand_frags(ARA, 3000, "ara", 1) + amp_frags(0.0))
write("robusta.fa", rand_frags(CAN, 3000, "rob", 2) + [f">amp_{p['marker_id']}_{i}\n{AMP[p['marker_id']]}" for p in primers if p["type"] in ("M1-Cc", "REF") and p["marker_id"] in AMP for i in range(10)])
write("blend_10pct.fa", rand_frags(ARA, 2700, "ara", 3) + rand_frags(CAN, 300, "rob", 33) + amp_frags(0.10))
write("blend_25pct.fa", rand_frags(ARA, 2250, "ara", 4) + rand_frags(CAN, 750, "rob", 44) + amp_frags(0.25))
write("blend_50pct.fa", rand_frags(ARA, 1500, "ara", 5) + rand_frags(CAN, 1500, "rob", 55) + amp_frags(0.50))
write("arabica_with_carbonarius.fa", rand_frags(ARA, 2800, "ara", 6) + rand_frags(CARB, 200, "mold", 66) + amp_frags(0.0, mold=True))
write("mold_free.fa", rand_frags(ARA, 3000, "ara", 7) + amp_frags(0.0, mold=False))
write("mold_ota_positive.fa", rand_frags(CARB, 3000, "mold", 8) + amp_frags(0.0, mold=True, skip_coffee=True))
