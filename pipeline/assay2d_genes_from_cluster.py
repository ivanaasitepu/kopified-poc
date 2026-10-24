#!/usr/bin/env python3
"""Assay 2d — tabel kehadiran gen OTA yang benar: ambil CDS dari klaster MW526250.1
(A. westerdijkiae OTA cluster, annotated), ekstrak gen ota*, BLASTN ke semua genom fungi."""
import gzip, os, re, subprocess, sys, time, urllib.request, ssl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

OUT = f"{BASE}/output/assay2"
GENEDIR = f"{BASE}/data/ota_genes"
os.makedirs(GENEDIR, exist_ok=True)
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

def http(u, tries=4):
    for t in range(tries):
        try:
            with urllib.request.urlopen(u, context=CTX, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            time.sleep(2 + t * 2)
    return ""

gb_path = f"{GENEDIR}/MW526250.gb"
fa_path = f"{GENEDIR}/MW526250.fasta"
if not os.path.exists(gb_path):
    open(gb_path, "w").write(http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MW526250.1&rettype=gb&retmode=text"))
if not os.path.exists(fa_path):
    open(fa_path, "w").write(http("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MW526250.1&rettype=fasta&retmode=text"))

# sekuens klaster
cluster = "".join(l.strip() for l in open(fa_path) if not l.startswith(">"))
print("panjang klaster:", len(cluster))

# parse CDS dari GenBank (Biopython, dukung join/complement)
from Bio import SeqIO
rec = SeqIO.read(gb_path, "genbank")
genes = {}
for feat in rec.features:
    if feat.type != "CDS":
        continue
    label = feat.qualifiers.get("gene", [None])[0] or feat.qualifiers.get("product", [f"cds_{feat.location.start}"])[0][:40]
    prod = feat.qualifiers.get("product", [""])[0]
    seq = str(feat.extract(rec.seq))
    genes[label] = (seq, prod)
print("gen di klaster:")
for k, (sq, prod) in genes.items():
    print(f"  {k:12} {len(sq):6} bp  {prod[:60]}")

# pilih gen target: otaA/otaB/otaD (atau nama produk yang cocok)
TARGET_KEYS = []
for k in genes:
    kl = k.lower()
    if kl.startswith(("ota","oat")) and kl.rstrip("1234567890")[-1] in "abcdr" or kl in ("otapks","otanrps"):
        TARGET_KEYS.append(k)
if not TARGET_KEYS:
    TARGET_KEYS = [k for k in genes if any(w in genes[k][1].lower() for w in ["polyketide", "peptide synthetase", "halogenase", "p450", "cytochrome"])]
TARGET_KEYS = TARGET_KEYS[:6]
print("gen dipakai:", TARGET_KEYS)

# tulis fasta per gen
gen_files = {}
for k in TARGET_KEYS:
    fp = f"{GENEDIR}/cluster_{k}.fasta"
    sq = genes[k][0]
    open(fp, "w").write(f">{k}|MW526250|{genes[k][1][:60]}\n" + "\n".join(sq[i:i+80] for i in range(0, len(sq), 80)) + "\n")
    gen_files[k] = fp

# BLASTN ke semua fungi
G = find_genomes()
FUNGI = [f for f, (sp, grp) in sorted(G.items()) if grp in ("mold_target", "mold_nontarget")]
grp_of = {G[f][0]: G[f][1] for f in FUNGI}
rows = {}
for f in FUNGI:
    db = f"{OUT}/db_{f.split('_genomic')[0]}"
    if not os.path.exists(db + ".nsq"):
        subprocess.run(["makeblastdb", "-in", f"{GENOMES}/{f}", "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
    sp = G[f][0]
    rows.setdefault(sp, {})
    for gname, fp in gen_files.items():
        glen = len(genes[gname][0])
        bl = subprocess.run(["blastn", "-query", fp, "-db", db, "-outfmt", "6 pident length",
                             "-task", "blastn", "-evalue", "1e-10", "-max_target_seqs", "5", "-num_threads", "4"],
                            check=True, capture_output=True, text=True)
        tot, bestp = 0, 0.0
        for line in bl.stdout.splitlines():
            p, l = line.split("\t")
            if int(l) >= 150:
                tot += int(l); bestp = max(bestp, float(p))
        cov_pct = min(100, tot / glen * 100)
        rows[sp][gname] = f"{bestp:.0f}%id/{cov_pct:.0f}%cov" if tot > 0 else "-"

hdr = ["genome", "group"] + TARGET_KEYS
lines = ["\t".join(hdr)]
for sp, d in rows.items():
    lines.append("\t".join([sp, grp_of[sp]] + [d.get(g, "-") for g in TARGET_KEYS]))
open(f"{BASE}/results/ota_gene_presence.tsv", "w").write("\n".join(lines) + "\n")
print("\n" + "\n".join(lines))
