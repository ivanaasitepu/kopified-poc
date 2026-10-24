#!/usr/bin/env python3
"""Assay 2c — tabel kehadiran GEN otaA/otaB/otaD per genom kapang (sesuai doc Drey Tahap 4).
Ambil sekuens gen dari NCBI, BLASTN ke semua genom fungi, tulis results/ota_gene_presence.tsv."""
import json, os, subprocess, sys, time, urllib.parse, urllib.request, ssl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay2"
GENEDIR = f"{BASE}/data/ota_genes"
os.makedirs(GENEDIR, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

def http(u, tries=4):
    for t in range(tries):
        try:
            with urllib.request.urlopen(u, context=CTX, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2 + t * 2)
    return ""

QUERIES = {
    "otaA_pks": 'Aspergillus[Organism] AND ochratoxin[All Fields] AND polyketide synthase[All Fields]',
    "otaB_nrps": 'Aspergillus[Organism] AND ochratoxin[All Fields] AND nonribosomal peptide synthetase[All Fields]',
    "otaD_hal": 'Aspergillus[Organism] AND ochratoxin[All Fields] AND halogenase[All Fields]',
}
genes = {}
for gname, q in QUERIES.items():
    fa = f"{GENEDIR}/{gname}.fasta"
    if os.path.exists(fa) and os.path.getsize(fa) > 500:
        genes[gname] = fa; continue
    txt = http(f"{EU}/esearch.fcgi?db=nuccore&term={urllib.parse.quote(q)}&retmax=5&retmode=json")
    time.sleep(1)
    ids = json.loads(txt)["esearchresult"]["idlist"]
    got = False
    for gid in ids:
        rec = http(f"{EU}/efetch.fcgi?db=nuccore&id={gid}&rettype=fasta&retmode=text")
        time.sleep(1)
        if not rec.startswith(">"):
            continue
        seq = "".join(l for l in rec.splitlines() if not l.startswith(">"))
        if 800 < len(seq) < 30000:  # panjang gen wajar
            open(fa, "w").write(rec)
            genes[gname] = fa
            print(f"{gname}: {rec.splitlines()[0][:80]} ({len(seq)} bp)")
            got = True
            break
    if not got:
        print(f"{gname}: TIDAK DITEMUKAN")

if len(genes) < 3:
    print("PERINGATAN: gen tidak lengkap", genes)

# BLASTN tiap gen ke tiap genom fungi
G = find_genomes()
FUNGI = [f for f, (sp, grp) in sorted(G.items()) if grp in ("mold_target", "mold_nontarget")]
rows = {}
for f in FUNGI:
    db = f"{OUT}/db_{f.split('_genomic')[0]}"
    if not os.path.exists(db + ".nsq"):
        subprocess.run(["makeblastdb", "-in", f"{GENOMES}/{f}", "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
    sp = G[f][0]
    rows.setdefault(sp, {})
    for gname, fa in genes.items():
        bl = subprocess.run(["blastn", "-query", fa, "-db", db, "-outfmt", "6 pident length evalue",
                             "-evalue", "1e-10", "-max_target_seqs", "3", "-num_threads", "4"],
                            check=True, capture_output=True, text=True)
        best = None
        for line in bl.stdout.splitlines():
            p, l, e = line.split("\t")
            if int(l) >= 300:
                if best is None or float(p) > best[0]:
                    best = (float(p), int(l))
        rows[sp][gname] = f"{best[0]:.0f}%/{best[1]}bp" if best else "-"

hdr = ["genome"] + list(genes.keys())
lines = ["\t".join(hdr)]
for sp, d in rows.items():
    lines.append("\t".join([sp] + [d.get(g, "-") for g in genes]))
os.makedirs(f"{BASE}/results", exist_ok=True)
open(f"{BASE}/results/ota_gene_presence.tsv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("\ntersimpan:", f"{BASE}/results/ota_gene_presence.tsv")
