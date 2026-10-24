#!/usr/bin/env python3
"""Assay 5 — matriks validasi + simulasi campuran Robusta (doc Drey Tahap 5).
Input: results/insilico_pcr.tsv, results/primers_v0.1.tsv, kandidat lokus assay1.
Output: results/validation_matrix.tsv, results/validation_matrix.png,
        results/blend_simulation.tsv, results/blend_simulation.png, results/funnel.tsv
"""
import csv, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kopilib import *

RES = f"{BASE}/results"
OUT = f"{BASE}/output/assay1"
random.seed(42)

# ---------- matriks validasi ----------
EXPECT = {  # grup -> marker tipe yang diharapkan positif
    "arabica": {"M1-Ce", "M1-Cc", "REF"},
    "eugenioides": {"M1-Ce", "REF"},
    "canephora": {"M1-Cc", "REF"},
    "mold_target": {"OTA-B", "OTA-D"},
    "mold_nontarget": set(),
    "microflora": set(),
}
rows = list(csv.DictReader(open(f"{RES}/insilico_pcr.tsv"), delimiter="\t"))
genomes_cols = [c for c in rows[0] if c not in ("marker_id", "type")]
mtypes = sorted({r["type"] for r in rows})
mat, detail = [], []
from kopilib import MAP
ACC2GRP = {k: v[1] for k, v in MAP.items()}
for col in genomes_cols:
    sp, acc = col.rsplit("|", 1)
    acc_key = "_".join(acc.split("_")[:2]).rsplit(".", 1)[0]
    grp = ACC2GRP.get(acc_key, "?")
    exp = EXPECT.get(grp, set())
    row = [f"{sp} ({col.split('|')[1]})"]
    for mt in mtypes:
        obs = any(int(r[col] or 0) > 0 for r in rows if r["type"] == mt)
        expected = mt in exp
        if expected and obs: cell, ok = "OK+", True
        elif not expected and not obs: cell, ok = "OK-", True
        elif expected and not obs: cell, ok = "MISS", False
        else: cell, ok = "FP", False
        row.append(cell)
        detail.append((sp, mt, cell))
    mat.append(row)

with open(f"{RES}/validation_matrix.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genom"] + mtypes)
    w.writerows(mat)

n_ok = sum(1 for _, _, c in detail if c in ("OK+", "OK-"))
print(f"matriks validasi: {n_ok}/{len(detail)} sel sesuai harapan ({100*n_ok/len(detail):.0f}%)")
for sp, mt, c in detail:
    if c in ("MISS", "FP"):
        print("  ketidaksesuaian:", sp, mt, c)

# heatmap
import matplotlib
matplotlib.use("Agg")
from matplotlib import font_manager
font_manager.fontManager.addfont("/public/skills/skill_tools/font/fonts/NotoSansCJKsc-Regular.otf")
matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
import numpy as np

CMAP = {"OK+": 2, "OK-": 1, "MISS": -1, "FP": -2}
grid = np.array([[CMAP[c] for c in row[1:]] for row in mat])
fig, ax = plt.subplots(figsize=(1.2 + 1.4 * len(mtypes), 0.5 * len(mat) + 2))
colors = {-2: "#d62728", -1: "#ff9896", 1: "#d6e8d4", 2: "#2ca02c"}
for i in range(grid.shape[0]):
    for j in range(grid.shape[1]):
        ax.add_patch(plt.Rectangle((j, i), 1, 1, color=colors[grid[i, j]], ec="white"))
        ax.text(j + 0.5, i + 0.5, mat[i][j + 1], ha="center", va="center", fontsize=8)
ax.set_xlim(0, len(mtypes)); ax.set_ylim(len(mat), 0)
ax.set_xticks([j + 0.5 for j in range(len(mtypes))]); ax.set_xticklabels(mtypes)
ax.set_yticks([i + 0.5 for i in range(len(mat))])
ax.set_yticklabels([r[0] for r in mat], fontsize=8)
ax.set_title("Matriks validasi in-silico (OK+/OK- = sesuai harapan)", fontsize=10)
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
plt.tight_layout(); plt.savefig(f"{RES}/validation_matrix.png", dpi=150)
print("heatmap tersimpan")

# ---------- simulasi campuran ----------
# lokus marker pada genom referensi
def loci_from_fasta(path, prefix):
    loc = {}
    for name, seq in iter_fasta(path):
        if not name.startswith(prefix):
            continue
        chrom, se = name.split("|")[1].split(":")
        s, e = se.split("-")
        loc.setdefault(chrom, []).append((int(s), int(e)))
    return {k: merge_intervals(v) for k, v in loc.items()}

ce_loci = loci_from_fasta(f"{OUT}/assay1_candidates.fasta", "A1_")  # di eug ref
cc_loci = loci_from_fasta(f"{OUT}/candidates_cc_ref.fasta", "Cc_")  # di can ref
rf_loci = loci_from_fasta(f"{OUT}/candidates_cc_ref.fasta", "REF_")  # di eug ref
# batasi ke marker terpilih (yang masuk panel primer)
sel_ids = {r["marker_id"] for r in csv.DictReader(open(f"{RES}/primers_v0.1.tsv"), delimiter="\t")}

def load_loci_subset(fasta, prefix):
    loc = {}
    for name, seq in iter_fasta(fasta):
        mid = name.split("|")[0]
        if mid not in sel_ids:
            continue
        chrom, se = name.split("|")[1].split(":")
        s, e = se.split("-")
        loc.setdefault(chrom, []).append((int(s), int(e)))
    return loc

ce_sel = load_loci_subset(f"{OUT}/assay1_candidates.fasta", "A1_")
cc_sel = load_loci_subset(f"{OUT}/candidates_cc_ref.fasta", "Cc_")
rf_sel = load_loci_subset(f"{OUT}/candidates_cc_ref.fasta", "REF_")
print(f"lokus simulasi: Ce={sum(len(v) for v in ce_sel.values())}, Cc={sum(len(v) for v in cc_sel.values())}, REF={sum(len(v) for v in rf_sel.values())}")

def simulate(robusta_frac, cells=100000):
    """Model copy-number: arabika = tetraploid (Ce 2 kopi + Cc 2 kopi + REF 4 kopi per sel),
    robusta = diploid (Cc 2 kopi + REF 2 kopi). Sinyal ~ Poisson."""
    import math
    nara = cells * (1 - robusta_frac)
    nrob = cells * robusta_frac
    def pois(lam):
        # poisson sederhana via normal approx utk lam besar
        return max(0, int(random.gauss(lam, math.sqrt(max(lam, 1)))))
    ce = pois(nara * 2)
    cc = pois(nara * 2 + nrob * 2)
    ref = pois(nara * 4 + nrob * 2)
    return {"Ce": ce, "Cc": cc, "REF": ref}

sim_rows = []
for pct in [0, 2, 5, 10, 25, 50, 100]:
    h = simulate(pct / 100)
    ratio = h["Cc"] / h["Ce"] if h["Ce"] else float("nan")
    sim_rows.append((pct, h["Ce"], h["Cc"], h["REF"], ratio))
    print(f"robusta {pct:3}%: Ce={h['Ce']:5} Cc={h['Cc']:5} REF={h['REF']:5} rasio Cc/Ce={ratio:.3f}")

with open(f"{RES}/blend_simulation.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["robusta_pct", "frag_Ce", "frag_Cc", "frag_REF", "ratio_Cc_Ce"])
    w.writerows(sim_rows)

fig, ax = plt.subplots(figsize=(6, 4))
xs = [r[0] for r in sim_rows]; ys = [r[4] for r in sim_rows]
ax.plot(xs, ys, "o-", color="#6f4e37")
ax.set_xlabel("% Robusta dalam campuran")
ax.set_ylabel("Rasio sinyal Cc/Ce (simulasi, belum dikalibrasi)")
ax.set_title("Kerangka rasio marker: simulasi campuran in-silico")
ax.grid(alpha=0.3)
plt.tight_layout(); plt.savefig(f"{RES}/blend_simulation.png", dpi=150)
print("kurva simulasi tersimpan")

# ---------- funnel ----------
fun = [("assembly genom terunduh & tervalidasi", 25),
       ("kandidat lokus M1-Ce (alignment)", sum(1 for n, _ in iter_fasta(f"{OUT}/assay1_candidates.fasta"))),
       ("kandidat M1-Cc + REF", sum(1 for n, _ in iter_fasta(f"{OUT}/candidates_cc_ref.fasta"))),
       ("marker terpilih untuk primer", len(sel_ids)),
       ("pasangan primer lolos desain+insilico", len(rows))]
with open(f"{RES}/funnel.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["tahap", "jumlah"])
    w.writerows(fun)
print("funnel:", fun)
