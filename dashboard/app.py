import os, sys, time
import streamlit as st
import pandas as pd

st.set_page_config(page_title="KopiFied — Panel Marker DNA", page_icon="☕", layout="wide")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = f"{BASE}/results"
SAMPLES = f"{BASE}/dashboard/samples"
sys.path.insert(0, f"{BASE}/dashboard")
from insilico import count_amplicons

st.title("☕ KopiFied — Panel Marker DNA v0.1")
st.caption("In-silico proof of concept; persentase campuran dan kadar OTA belum dikalibrasi (TRL 3).")

PANEL = f"{RES}/primers_v0.1.tsv"
if not os.path.exists(PANEL):
    st.warning("Panel primer v0.1 belum tersedia — pipeline masih berjalan.")
    st.stop()
primers = pd.read_csv(PANEL, sep="\t")

# kurva simulasi + ambang
sim, THRESH = None, 1.6
if os.path.exists(f"{RES}/blend_simulation.tsv"):
    sim = pd.read_csv(f"{RES}/blend_simulation.tsv", sep="\t")
    r0 = float(sim.loc[sim.robusta_pct == 0, "ratio_Cc_Ce"].iloc[0])
    r25 = float(sim.loc[sim.robusta_pct == 25, "ratio_Cc_Ce"].iloc[0])
    THRESH = r0 + (r25 - r0) * 0.3

# label naratif sampel
LABELS = {
    "arabica_caturra.fa": ("🌱 Arabika Gayo — sampel pasar (seharusnya ASLI)", "arabika murni"),
    "mold_free.fa": ("🌱 Arabika bersih — tanpa kapang", "arabika murni, negatif OTA"),
    "blend_10pct.fa": ("⚠️ Campuran 10% Robusta — fraud halus", "campuran 10%"),
    "blend_25pct.fa": ("⚠️ Campuran 25% Robusta — kasus fraud tipikal", "campuran 25%"),
    "blend_50pct.fa": ("🚨 Campuran 50% Robusta — fraud berat", "campuran 50%"),
    "robusta.fa": ("🫘 Robusta murni — dijual sebagai 'arabika'", "robusta murni"),
    "arabica_with_carbonarius.fa": ("🍄 Arabika + kapang A. carbonarius (OTA+)", "arabika terkontaminasi kapang"),
    "mold_ota_positive.fa": ("🧫 Kontrol: ekstrak kapang murni (bukan kopi)", "kontrol kapang"),
}

st.sidebar.header("Sampel uji")
preset = [f for f in sorted(os.listdir(SAMPLES)) if f.endswith(".fa")]
disp = [f"{LABELS.get(f, (f,))[0]}" for f in preset]
idx = st.sidebar.selectbox("Pilih sampel", range(len(preset)), format_func=lambda i: disp[i])
if preset and preset[idx] in LABELS:
    st.sidebar.caption(LABELS[preset[idx]][1])
up = st.sidebar.file_uploader("atau unggah FASTA sendiri (maks 20 MB)", type=["fa", "fasta", "txt"])
run = st.sidebar.button("▶ Jalankan Panel", type="primary")

def read_fasta(src):
    seqs, buf = [], []
    fh = open(src) if isinstance(src, str) else src.read().decode().splitlines()
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if buf: seqs.append("".join(buf).upper())
            buf = []
        elif line:
            buf.append(line)
    if buf: seqs.append("".join(buf).upper())
    return seqs

def sim_figure(ratio_sample):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import font_manager
    _f = "/public/skills/skill_tools/font/fonts/NotoSansCJKsc-Regular.otf"
    if os.path.exists(_f):
        font_manager.fontManager.addfont(_f)
    matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK SC"]
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    xs = sim["robusta_pct"]; ys = sim["ratio_Cc_Ce"]
    ax.plot(xs, ys, "o-", color="#6f4e37", label="kurva simulasi (belum dikalibrasi)")
    ax.axhline(THRESH, color="#d62728", ls="--", lw=1.2, label=f"ambang flag ({THRESH:.2f})")
    if ratio_sample == ratio_sample:
        ax.scatter([0], [ratio_sample], s=140, color="#d62728", zorder=5, marker="*",
                   label=f"sampel Anda (rasio {ratio_sample:.2f})")
    ax.set_xlabel("% Robusta dalam campuran"); ax.set_ylabel("Rasio Cc/Ce")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    out = f"{RES}/_dash_sim.png"; fig.savefig(out, dpi=130); plt.close(fig)
    return out

if run:
    if up is not None and up.size > 20 * 1024 * 1024:
        st.error("File terlalu besar (maks 20 MB). Untuk genom utuh, gunakan pipeline offline.")
        st.stop()
    with st.spinner("Menjalankan panel 16 marker..."):
        t0 = time.time()
        src = up if up else f"{SAMPLES}/{preset[idx]}"
        label = getattr(up, "name", None) or preset[idx]
        seqs = read_fasta(src)
        hits = {p.marker_id: count_amplicons(seqs, p.fwd, p.rev, max_mm=2) for _, p in primers.iterrows()}
        dt = time.time() - t0

    ce = sum(hits.get(m, 0) for m in primers.marker_id[primers.type == "M1-Ce"])
    cc = sum(hits.get(m, 0) for m in primers.marker_id[primers.type == "M1-Cc"])
    ref_rows = primers[primers.type == "REF"]
    ref = sum(hits.get(m, 0) for m in ref_rows.marker_id)
    ref_n = sum(1 for m in ref_rows.marker_id if hits.get(m, 0) > 0)
    ota = {m: hits.get(m, 0) for m in primers.marker_id[primers.type.str.startswith("OTA")]}
    n_ce = max(1, int((primers.type == "M1-Ce").sum())); n_cc = max(1, int((primers.type == "M1-Cc").sum()))
    ratio = (cc / n_cc) / (ce / n_ce) if ce else float("nan")

    if ref == 0 and ce == 0 and cc == 0:
        call1, col1 = "INVALID — DNA kopi tidak terdeteksi", "red"
    elif ce == 0:
        call1, col1 = "FLAG: bukan Arabika (kemungkinan Robusta murni)", "red"
    elif ratio == ratio and ratio > THRESH:
        call1, col1 = f"FLAG: kemungkinan campuran Robusta (rasio {ratio:.2f} > ambang {THRESH:.2f})", "orange"
    else:
        call1, col1 = "ARABICA AUTHENTIC ✅", "green"

    ota_pos = [m for m, n in ota.items() if n > 0]
    if len(ota_pos) >= 2:
        call2, col2 = "FLAG: DNA kapang penghasil OTA terdeteksi — kirim lot untuk uji HPLC OTA", "red"
    elif len(ota_pos) == 1:
        call2, col2 = "CHECK: satu marker OTA positif — periksa ulang", "orange"
    else:
        call2, col2 = "Tidak ada DNA kapang penghasil OTA terdeteksi ✅", "green"

    st.success(f"Sampel: **{label}** — {len(seqs)} fragmen — dianalisis dalam {dt:.1f} detik")

    st.subheader("Modul 1 — Autentikasi Arabika")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("M1-Ce (bukti Arabika)", ce)
    c2.metric("M1-Cc (komponen Robusta)", cc)
    c3.metric("REF (kontrol kopi)", f"{ref} ({ref_n}/{len(ref_rows)} lokus)")
    c4.metric("Rasio Cc/Ce", f"{ratio:.2f}" if ratio == ratio else "—")
    st.markdown(f"### :{col1}[{call1}]")
    if sim is not None:
        st.image(sim_figure(ratio))
    st.caption("⚠️ Rasio dan ambang berasal dari simulasi copy-number yang BELUM dikalibrasi dengan sampel fisik — bukan klaim batas deteksi.")

    st.subheader("Modul 2 — Risiko Kapang OTA")
    cols = st.columns(max(1, len(ota)))
    for c, (m, n) in zip(cols, ota.items()):
        c.metric(m, n)
    st.markdown(f"### :{col2}[{call2}]")

    with st.expander("🔬 Rincian per marker (hit amplikon in-silico)"):
        det = pd.DataFrame([{"marker": m, "tipe": primers.set_index("marker_id").loc[m, "type"],
                             "amplikon": n} for m, n in hits.items()])
        st.dataframe(det, hide_index=True)

    html = f"""<html><body style="font-family:sans-serif">
    <h2>KopiFied — Laporan QA (v0.1)</h2>
    <p><b>Sampel:</b> {label} | <b>Fragmen:</b> {len(seqs)} | <b>Waktu analisis:</b> {dt:.1f}s</p>
    <h3>Modul 1 — Autentikasi Arabika</h3>
    <p>M1-Ce={ce}, M1-Cc={cc}, REF={ref} ({ref_n}/{len(ref_rows)} lokus), rasio Cc/Ce={'%.2f' % ratio if ratio == ratio else '—'}<br><b>{call1}</b></p>
    <h3>Modul 2 — Risiko Kapang OTA</h3>
    <p>{', '.join(f'{m}={n}' for m, n in ota.items())}<br><b>{call2}</b></p>
    <hr><p><i>In-silico proof of concept (TRL 3). Persentase campuran dan kadar OTA belum dikalibrasi; keputusan akhir memerlukan validasi lab.</i></p>
    </body></html>"""
    st.download_button("⬇ Unduh Laporan QA (HTML)", html, file_name=f"KopiFied_QA_{label.replace('.fa','')}_v0.1.html")

# ---------- tab bukti ----------
st.divider()
t1, t2, t3, t5, t4 = st.tabs(["Matriks Validasi", "Primer", "Genom", "Bukti SRA & Multiplex", "Metode & Batasan"])
with t1:
    p = f"{RES}/validation_matrix.png"
    if os.path.exists(p):
        st.image(p)
        st.caption("Hijau = sesuai harapan; merah = tidak sesuai. Sumber: results/validation_matrix.tsv")
with t2:
    st.dataframe(primers, hide_index=True)
with t3:
    g = f"{BASE}/config/genomes.tsv"
    if os.path.exists(g):
        st.dataframe(pd.read_csv(g, sep="\t"), hide_index=True)
with t5:
    st.markdown("**Validasi pada sampel independen (SRA SRR40981778, arabika lini mutan azide)**")
    sp = f"{BASE}/output/sra_check/sra_primer_sites.tsv"
    if os.path.exists(sp):
        st.dataframe(pd.read_csv(sp, sep="\t"), hide_index=True)
        st.caption("9/14 marker memiliki ≥1 read cocok primer penuh (≤1 mismatch) pada ~1,3× coverage; marker OTA nol di reads kopi (sesuai harapan).")
    st.markdown("**Kompatibilitas multiplex (primer3.calc_heterodimer)**")
    mp = f"{RES}/multiplex_check.tsv"
    if os.path.exists(mp):
        mux = pd.read_csv(mp, sep="\t")
        st.write(f"Hanya {len(mux)} dari 1.035 pasangan oligo (0,8%) berisiko dimer (ΔG < −9 kkal/mol) — semuanya interaksi probe–probe, dipisah ke reaksi berbeda.")
        st.dataframe(mux, hide_index=True)
with t4:
    r = f"{BASE}/README.md"
    if os.path.exists(r):
        st.markdown(open(r).read())
