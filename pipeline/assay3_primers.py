#!/usr/bin/env python3
"""Assay 3 — desain primer + probe TaqMan (primer3) + in-silico PCR untuk semua marker.
Input: output/assay1/assay1_marker_panel.csv (M1-Ce), candidates_cc_ref.fasta (M1-Cc, REF),
       gen otaB/otaD dari klaster MW526250 (marker OTA, jendela konservatif).
Output: results/primers_v0.1.tsv, results/insilico_pcr.tsv
"""
import csv, os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *
import primer3

OUT = f"{BASE}/output/assay1"
RES = f"{BASE}/results"
os.makedirs(RES, exist_ok=True)
G = find_genomes()

def parse_hdr_fasta(path):
    """return list of (id, chrom, start, end, length, gc, seq)"""
    out = []
    for name, seq in iter_fasta(path):
        mid, loc = name.split("|")[0], name.split("|")[1]
        chrom, se = loc.split(":")
        s, e = se.split("-")
        ln = int(name.split("len=")[1].split("|")[0])
        gc = float(name.split("gc=")[1].split("|")[0])
        out.append(dict(id=mid, chrom=chrom, start=int(s), end=int(e), len=ln, gc=gc, seq=seq))
    return out

# ---------------- 1. pilih lokus ----------------
markers = []  # dict: id, type, seq, expect(groups yang harus amplifikasi), forbid(groups yang tidak boleh)

# --- M1-Ce: dari panel top200 (ara 2/2 id tinggi, canephora divergen <=88%, eug hit)
panel = f"{OUT}/assay1_marker_panel_top200.csv"
seqs_ce = {name.split("|")[0]: seq for name, seq in iter_fasta(f"{OUT}/assay1_candidates.fasta")}
good_ce = []
with open(panel) as fh:
    r = csv.reader(fh)
    hdr = next(r)
    can_cols = [i for i, h in enumerate(hdr) if h == "C. canephora"]
    ara_cols = [i for i, h in enumerate(hdr) if h == "C. arabica"]
    eug_cols = [i for i, h in enumerate(hdr) if h == "C. eugenioides"]
    for row in r:
        def maxid(cols):
            v = [float(row[i].split("%")[0]) for i in cols if row[i] != "-"]
            return max(v) if v else 0.0
        n_ara = sum(row[i] != "-" for i in ara_cols)
        if n_ara == len(ara_cols) and maxid(ara_cols) >= 95 and maxid(eug_cols) >= 90 and 0 < maxid(can_cols) <= 88:
            good_ce.append((row[0], maxid(can_cols)))
good_ce.sort(key=lambda x: x[1])  # identitas canephora terendah dulu
good_ce = [g[0] for g in good_ce[:6]]
for mid in good_ce:
    markers.append(dict(id=mid, type="M1-Ce", seq=seqs_ce[mid],
                        expect={"arabica", "eugenioides"}, forbid={"canephora", "mold_target", "mold_nontarget"}))
print("M1-Ce dipilih:", good_ce)

# --- M1-Cc & REF dari candidates_cc_ref.fasta
ccref = parse_hdr_fasta(f"{OUT}/candidates_cc_ref.fasta")
cc = sorted([c for c in ccref if c["id"].startswith("Cc_") and 0.38 <= c["gc"] <= 0.58 and c["len"] >= 800], key=lambda x: -x["len"])[:4]
rf = sorted([c for c in ccref if c["id"].startswith("REF_") and 0.38 <= c["gc"] <= 0.58 and c["len"] >= 800], key=lambda x: -x["len"])[:3]
for c in cc:
    markers.append(dict(id=c["id"], type="M1-Cc", seq=c["seq"],
                        expect={"arabica", "canephora"}, forbid={"eugenioides", "mold_target", "mold_nontarget"}))
for c in rf:
    markers.append(dict(id=c["id"], type="REF", seq=c["seq"],
                        expect={"arabica", "canephora", "eugenioides"}, forbid={"mold_target", "mold_nontarget"}))
print(f"M1-Cc: {[c['id'] for c in cc]}  REF: {[c['id'] for c in rf]}")

# --- OTA: jendela konservatif otaB & otaD (identitas >=80% pada semua produsen yang membawa gen)
def conserved_windows(gene_fa, producer_dbs, win=250, min_id=80.0):
    seq = next(iter_fasta(gene_fa))[1]
    covs = []
    for db in producer_dbs:
        bl = subprocess.run(["blastn", "-query", gene_fa, "-db", db, "-task", "blastn",
                             "-outfmt", "6 qstart qend pident", "-evalue", "1e-10",
                             "-num_threads", "4"], check=True, capture_output=True, text=True)
        cov = {}
        for line in bl.stdout.splitlines():
            qs, qe, p = line.split("\t")
            if float(p) >= min_id:
                for i in range(min(int(qs), int(qe)) - 1, max(int(qs), int(qe))):
                    cov[i] = 1
        covs.append(cov)
    # interseksi posisi
    inter = set(covs[0])
    for c in covs[1:]:
        inter &= set(c)
    pos = sorted(inter)
    wins, i = [], 0
    while i + win <= len(pos) and pos:
        s = pos[i]; j = i
        while j < len(pos) and pos[j] - s < win:
            j += 1
        if j - i >= win * 0.8:
            wins.append((s, s + win))
        i = j
    return seq, wins[:3]

GENEDIR = f"{BASE}/data/ota_genes"
OUT2 = f"{BASE}/output/assay2"
producer_dbs = [f"{OUT2}/db_GCA_001990825.1_Aspca3", f"{OUT2}/db_GCF_002849105.1_Aspste1"]  # carbonarius + steynii
for gname, mtype in [("cluster_otaB", "OTA-B"), ("cluster_otaD", "OTA-D")]:
    fp = f"{GENEDIR}/{gname}.fasta"
    if not os.path.exists(fp):
        continue
    seq, wins = conserved_windows(fp, producer_dbs)
    print(f"{gname}: {len(wins)} jendela konservatif", wins)
    for k, (s, e) in enumerate(wins[:2]):
        wseq = seq[s:e]
        gc = (wseq.count("G") + wseq.count("C")) / len(wseq)
        if 0.38 <= gc <= 0.62:
            markers.append(dict(id=f"{mtype}-{k+1}", type=mtype, seq=wseq,
                                expect={"mold_target"}, forbid={"mold_nontarget", "arabica", "canephora", "eugenioides", "microflora"}))

# ---------------- 2. desain primer ----------------
ARGS = dict(PRIMER_NUM_RETURN=3, PRIMER_PRODUCT_SIZE_RANGE=[[70, 130]],
            PRIMER_OPT_SIZE=20, PRIMER_MIN_SIZE=18, PRIMER_MAX_SIZE=25,
            PRIMER_OPT_TM=60.0, PRIMER_MIN_TM=57.0, PRIMER_MAX_TM=63.0,
            PRIMER_OPT_GC_PERCENT=50.0, PRIMER_MIN_GC=40.0, PRIMER_MAX_GC=60.0,
            PRIMER_MAX_POLY_X=4, PRIMER_MAX_SELF_ANY_TH=45.0, PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0,
            PRIMER_PICK_INTERNAL_OLIGO=1, PRIMER_INTERNAL_OPT_TM=69.0,
            PRIMER_INTERNAL_MIN_TM=67.0, PRIMER_INTERNAL_MAX_TM=71.0, PRIMER_INTERNAL_MAX_POLY_X=3)

designed = []
for m in markers:
    try:
        res = primer3.bindings.design_primers(
            {"SEQUENCE_ID": m["id"], "SEQUENCE_TEMPLATE": m["seq"]}, ARGS)
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        for i in range(min(n, 1)):  # ambil pasangan terbaik
            designed.append(dict(
                marker_id=m["id"], type=m["type"],
                fwd=res[f"PRIMER_LEFT_{i}_SEQUENCE"], rev=res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                probe=res.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", ""),
                amplicon=res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                tm_f=round(res[f"PRIMER_LEFT_{i}_TM"], 1), tm_r=round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
                tm_p=round(res.get(f"PRIMER_INTERNAL_{i}_TM", 0), 1),
                expect=m["expect"], forbid=m["forbid"]))
    except Exception as e:
        print(f"primer3 gagal {m['id']}: {e}")
print("pasangan primer terdesain:", len(designed))

# ---------------- 3. in-silico PCR ----------------
RC = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(RC)[::-1]
def mismatches(a, b): return sum(1 for x, y in zip(a, b) if x != y)

def insilico_pcr(seqs, fwd, rev, max_mm=2, min_amp=50, max_amp=1000):
    """kembalikan jumlah amplikon valid dalam genome dict name->seq"""
    rrev = revcomp(rev)
    fa = fwd[-12:]  # anchor 3' 12-mer
    ra = rrev[-12:]
    count = 0
    for name, genome_seq in seqs.items():
        L = len(genome_seq)
        i = genome_seq.find(fa)
        while i != -1:
            s = i + 12 - len(fwd)
            if s >= 0:
                sub = genome_seq[s:s+len(fwd)]
                if len(sub) == len(fwd) and mismatches(sub, fwd) <= max_mm and mismatches(sub[-3:], fwd[-3:]) == 0:
                    # cari primer rev di hilir
                    j = genome_seq.find(ra, s)
                    while j != -1 and j - s <= max_amp + len(rev):
                        rs = j + 12 - len(rev)
                        if rs >= 0:
                            rsub = genome_seq[rs:rs+len(rev)]
                            if len(rsub) == len(rev) and mismatches(rsub, rrev) <= max_mm and mismatches(rsub[-3:], rrev[-3:]) == 0:
                                amp = j + 12 - s
                                if min_amp <= amp <= max_amp:
                                    count += 1
                        j = genome_seq.find(ra, j + 1)
            i = genome_seq.find(fa, i + 1)
    return count

print("memuat genom untuk in-silico PCR ...")
genome_cache = {}
def load_genome(f):
    if f not in genome_cache:
        genome_cache[f] = dict(iter_fasta(f"{GENOMES}/{f}"))
    return genome_cache[f]

ins_rows = []
for d in designed:
    row = {"marker_id": d["marker_id"], "type": d["type"]}
    for f, (sp, grp) in sorted(G.items()):
        if d["type"] in ("M1-Ce", "M1-Cc", "REF"):
            # marker kopi: uji semua kopi + off-target kunci (niger, flavus[Pks mirip], cerevisiae)
            keep = grp in ("arabica", "canephora", "eugenioides") or "niger" in sp.lower() \
                   or "flavus" in sp.lower() or "cerevisiae" in sp.lower()
        else:
            # marker OTA: uji semua fungi + 2 genom kopi sebagai off-target
            keep = grp in ("mold_target", "mold_nontarget", "microflora") or "Cara_1.0" in f or "AUK" in f
        if not keep:
            continue
        n = insilico_pcr(load_genome(f), d["fwd"], d["rev"])
        row[f"{sp}|{grp}"] = n
    ins_rows.append(row)
    obs_grp = {k.split("|")[1] for k, v in row.items() if k not in ("marker_id", "type") and isinstance(v, int) and v > 0}
    tested = {k.split("|")[1] for k in row if k not in ("marker_id", "type")}
    ok_expect = all((e in obs_grp) or (e not in tested) for e in d["expect"])
    forbid_hit = any(g in obs_grp for g in d["forbid"])
    d["lolos"] = ok_expect and not forbid_hit
    print(d["marker_id"], d["type"], "obs:", obs_grp, "| expect:", d["expect"], "| lolos:", d["lolos"])

n_lolos = sum(1 for d in designed if d.get("lolos"))
print(f"primer lolos in-silico: {n_lolos}/{len(designed)}")

# ---------------- 4. simpan ----------------
def write_panel(path, items):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["marker_id", "type", "fwd", "rev", "probe", "amplicon_bp", "tm_fwd", "tm_rev", "tm_probe"])
        for d in items:
            w.writerow([d["marker_id"], d["type"], d["fwd"], d["rev"], d["probe"], d["amplicon"],
                        d["tm_f"], d["tm_r"], d["tm_p"]])
write_panel(f"{RES}/primers_v0.1.tsv", [d for d in designed if d.get("lolos")])
write_panel(f"{RES}/primers_all_debug.tsv", designed)
cols = ["marker_id", "type"] + sorted({k for r in ins_rows for k in r if k not in ("marker_id", "type")})
with open(f"{RES}/insilico_pcr.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for r in ins_rows:
        w.writerow([r.get(c, "") for c in cols])
print("tersimpan:", f"{RES}/primers_v0.1.tsv", "dan", f"{RES}/insilico_pcr.tsv")
