#!/usr/bin/env python3
"""Assay 2b — OTA biosynthetic-pathway detection markers.
Rationale: detect the toxin pathway (ota genes), not just species — biologically the
correct early-warning signal; also solves the sister-species problem
(A. ochraceus vs A. westerdijkiae) by grouping all OTA producers.
"""
import csv, os, subprocess, sys, time, json, urllib.parse, urllib.request, ssl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay2b"
os.makedirs(OUT, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

def get(url, tries=5):
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, context=CTX, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2 + t * 2)
    return ""

# 1. fetch ota gene/cluster sequences from NCBI (skip if already curated)
ota_fa = f"{OUT}/ota_genes.fasta"
if not os.path.exists(ota_fa):
    IDS = ["2695619581", "2695596150", "259486486", "350640282", "350631947",
           "2026808933", "2026808927", "2026808921", "2026808916", "2026808910"]
    meta = []
    with open(ota_fa, "w") as fh:
        for uid in IDS:
            txt = get(f"{EU}/efetch.fcgi?db=nuccore&id={uid}&rettype=fasta&retmode=text")
            time.sleep(0.9)
            lines = [l for l in txt.splitlines() if l.startswith(">") or (l and not l.startswith("{"))]
            if not lines or not lines[0].startswith(">"):
                print(uid, "fetch failed"); continue
            fh.write("\n".join(lines) + "\n")
            meta.append(lines[0][1:100])
            print("fetched:", lines[0][:100])
    print("ota sequences:", len(meta))
else:
    print("using existing curated", ota_fa)

# 2. blastn vs all fungal genomes
G = find_genomes()
dbs = sorted(G)
res = {}
for f in dbs:
    sp, grp = G[f]
    db = f"{BASE}/output/assay2/db_{f.split('_genomic')[0]}"
    if not os.path.exists(db + ".nsq"):
        subprocess.run(["makeblastdb", "-in", f"{GENOMES}/{f}", "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
    bl = subprocess.run(["blastn", "-query", ota_fa, "-db", db, "-outfmt",
                         "6 qseqid pident length evalue", "-evalue", "1e-10", "-num_threads", "8", "-dust", "no"],
                        check=True, capture_output=True, text=True)
    hits = [l.split("\t") for l in bl.stdout.splitlines()]
    res[f] = hits
    print(f"{sp:22} {f[:26]:26} hits={len(hits)}")

# 3. presence matrix
queries = sorted({x[0] for f in dbs for x in res[f]})
with open(f"{OUT}/ota_presence_matrix.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ota_query", "species", "assembly", "n_hits", "best_pident", "best_len"])
    for f in dbs:
        sp, grp = G[f]
        byq = {}
        for q, p, l, e in res[f]:
            b = byq.get(q)
            if not b or float(e) < float(b[2]):
                byq[q] = (float(p), int(l), e)
        for q in queries:
            b = byq.get(q)
            if b:
                w.writerow([q, sp, f.split("_genomic")[0], sum(1 for x in res[f] if x[0] == q), f"{b[0]:.1f}", b[1]])
print("saved:", f"{OUT}/ota_presence_matrix.csv")
