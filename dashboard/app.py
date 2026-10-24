import os, sys, time
import streamlit as st
import pandas as pd

st.set_page_config(page_title="KopiFied — Panel Marker DNA", page_icon="☕", layout="wide",
                   initial_sidebar_state="expanded")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = f"{BASE}/results"
SAMPLES = f"{BASE}/dashboard/samples"
sys.path.insert(0, f"{BASE}/dashboard")
from insilico import count_amplicons

# ================= DESIGN SYSTEM ============
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Sora:wght@400;600;800&family=IBM+Plex+Mono:wght@400;600&display=swap');

:root {
  --coffee-900:#241610; --coffee-700:#3E2723; --coffee-500:#6f4e37; --coffee-300:#a98467;
  --cream:#f6f1e7; --gold:#c8a165; --green:#2e7d32; --red:#c62828; --amber:#ef6c00;
}
html, body, [class*="css"] { font-family:'Sora',sans-serif; }
.stApp { background: var(--cream); }

/* hero */
.hero { background: linear-gradient(120deg, var(--coffee-900) 0%, var(--coffee-700) 45%, var(--coffee-500) 100%);
  border-radius:18px; padding:34px 38px; color:#fff; margin-bottom:8px; position:relative; overflow:hidden;}
.hero::after{ content:""; position:absolute; right:-60px; top:-60px; width:260px; height:260px;
  background:radial-gradient(circle, rgba(200,161,101,.35), transparent 70%); }
.hero h1{ font-weight:800; font-size:34px; margin:0; letter-spacing:-.5px;}
.hero .sub{ opacity:.9; margin-top:6px; font-size:15px;}
.badges{ margin-top:14px; display:flex; gap:8px; flex-wrap:wrap;}
.badge{ background:rgba(255,255,255,.14); border:1px solid rgba(255,255,255,.25); padding:4px 12px;
  border-radius:999px; font-size:12px; font-weight:600;}
.badge.warn{ background:rgba(239,108,0,.25); border-color:rgba(239,108,0,.6);}

/* decision banner */
.verdict{ border-radius:14px; padding:20px 26px; margin:14px 0; color:#fff; display:flex; align-items:center; gap:18px;}
.verdict-ok{ background:linear-gradient(120deg,#1b5e20,#2e7d32);}
.verdict-flag{ background:linear-gradient(120deg,#b71c1c,#c62828);}
.verdict-warn{ background:linear-gradient(120deg,#e65100,#ef6c00);}
.verdict .big{ font-size:22px; font-weight:800; }
.verdict .small{ font-size:13px; opacity:.92;}
.verdict .icon{ font-size:34px;}

/* metric cards */
.mcard{ background:#fff; border:1px solid #e7ddc9; border-radius:14px; padding:16px 18px;
  box-shadow:0 2px 8px rgba(62,39,35,.06); text-align:center;}
.mcard .v{ font-size:26px; font-weight:800; color:var(--coffee-700); font-family:'IBM Plex Mono',monospace;}
.mcard .l{ font-size:11.5px; color:#8a7a5f; font-weight:600; text-transform:uppercase; letter-spacing:.5px; margin-top:4px;}
.mcard .s{ font-size:11px; color:#a89778; margin-top:2px;}

.section-title{ font-size:19px; font-weight:800; color:var(--coffee-700); margin:22px 0 4px;
  display:flex; align-items:center; gap:10px;}
.section-title .n{ background:var(--coffee-500); color:#fff; width:26px; height:26px; border-radius:8px;
  display:inline-flex; align-items:center; justify-content:center; font-size:13px;}
.note{ background:#fff8ec; border:1px solid #eeddb8; border-radius:10px; padding:10px 14px; font-size:12.5px; color:#7a5c2e;}

/* sidebar */
section[data-testid="stSidebar"]{ background:var(--coffee-900);}
section[data-testid="stSidebar"] *{ color:#f3ead9 !important;}
section[data-testid="stSidebar"] .stButton button{ background:var(--gold) !important; color:var(--coffee-900) !important;
  font-weight:800; border:none; border-radius:10px;}
div[data-testid="stTabs"] button{ font-weight:600;}
</style>
""", unsafe_allow_html=True)

# ================= HEADER =================
st.markdown("""
<div class="hero">
  <h1>☕ KopiFied</h1>
  <div class="sub">Panel Marker DNA untuk Autentikasi Kopi Arabika &amp; Deteksi Risiko Kapang OTA — Tim Kimpulita, GenoHack 2026</div>
  <div class="badges">
    <span class="badge">Panel v0.1 · 16 primer</span>
    <span class="badge">25 genom publik tervalidasi</span>
    <span class="badge">TRL 3</span>
    <span class="badge warn">In-silico PoC — belum dikalibrasi lab</span>
  </div>
</div>
""", unsafe_allow_html=True)

# ================= DATA =================
PANEL = f"{RES}/primers_v0.1.tsv"
if not os.path.exists(PANEL):
    st.warning("Panel primer belum tersedia."); st.stop()
primers = pd.read_csv(PANEL, sep="\t")

sim, THRESH = None, 1.6
if os.path.exists(f"{RES}/blend_simulation.tsv"):
    sim = pd.read_csv(f"{RES}/blend_simulation.tsv", sep="\t")
    r0 = float(sim.loc[sim.robusta_pct == 0, "ratio_Cc_Ce"].iloc[0])
    r25 = float(sim.loc[sim.robusta_pct == 25, "ratio_Cc_Ce"].iloc[0])
    THRESH = r0 + (r25 - r0) * 0.3

LABELS = {
    "arabica_caturra.fa": ("🌱 Arabika Gayo — sampel pasar", "arabika murni · harapan: AUTHENTIC"),
    "mold_free.fa": ("🌱 Arabika bersih — tanpa kapang", "arabika murni · harapan: AUTHENTIC + bersih"),
    "blend_10pct.fa": ("⚠️ Campuran 10% Robusta — fraud halus", "campuran 10% · harapan: FLAG"),
    "blend_25pct.fa": ("⚠️ Campuran 25% Robusta — fraud tipikal", "campuran 25% · harapan: FLAG"),
    "blend_50pct.fa": ("🚨 Campuran 50% Robusta — fraud berat", "campuran 50% · harapan: FLAG"),
    "robusta.fa": ("🫘 Robusta murni dijual sebagai 'arabika'", "robusta murni · harapan: FLAG bukan arabika"),
    "arabica_with_carbonarius.fa": ("🍄 Arabika + A. carbonarius (OTA+)", "kontaminasi kapang · harapan: AUTHENTIC + FLAG OTA"),
    "mold_ota_positive.fa": ("🧫 Kontrol: ekstrak kapang murni", "bukan kopi · harapan: INVALID + FLAG OTA"),
}

st.sidebar.markdown("### 🧪 Sampel Uji")
preset = [f for f in sorted(os.listdir(SAMPLES)) if f.endswith(".fa")]
disp = [LABELS.get(f, (f, ""))[0] for f in preset]
idx = st.sidebar.selectbox("Pilih sampel siap-pakai", range(len(preset)), format_func=lambda i: disp[i])
if preset and preset[idx] in LABELS:
    st.sidebar.caption(LABELS[preset[idx]][1])
up = st.sidebar.file_uploader("atau unggah FASTA sendiri (maks 20 MB)", type=["fa", "fasta", "txt"])
run = st.sidebar.button("▶  JALANKAN PANEL", use_container_width=True)
st.sidebar.divider()
st.sidebar.caption("Aturan keputusan mengikuti panel v0.1; rasio berasal dari simulasi copy-number (belum dikalibrasi).")

def read_fasta(src):
    seqs, buf = [], []
    fh = open(src) if isinstance(src, str) else src.read().decode().splitlines()
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if buf: seqs.append("".join(buf).upper())
            buf = []
        elif line: buf.append(line)
    if buf: seqs.append("".join(buf).upper())
    return seqs

def ratio_gauge(ratio):
    import plotly.graph_objects as go
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=ratio,
        number={"font": {"size": 34, "family": "IBM Plex Mono"}, "valueformat": ".2f"},
        title={"text": "Rasio Cc/Ce", "font": {"size": 14}},
        gauge={"axis": {"range": [0, 2.5]},
               "bar": {"color": "#6f4e37"},
               "steps": [{"range": [0, THRESH], "color": "#dff0d8"},
                         {"range": [THRESH, 2.5], "color": "#f8d7da"}],
               "threshold": {"line": {"color": "#c62828", "width": 3}, "value": THRESH}}))
    fig.update_layout(height=210, margin=dict(l=20, r=20, t=40, b=10), paper_bgcolor="rgba(0,0,0,0)")
    return fig

def sim_chart(ratio):
    import plotly.graph_objects as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sim["robusta_pct"], y=sim["ratio_Cc_Ce"], mode="lines+markers",
                             name="kurva simulasi", line=dict(color="#6f4e37", width=3),
                             marker=dict(size=7)))
    fig.add_hline(y=THRESH, line_dash="dash", line_color="#c62828",
                  annotation_text=f"ambang flag ({THRESH:.2f})", annotation_position="top left")
    if ratio == ratio:
        fig.add_trace(go.Scatter(x=[0], y=[ratio], mode="markers", name="sampel Anda",
                                 marker=dict(size=16, color="#c62828", symbol="star")))
    fig.update_layout(height=260, margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_title="% Robusta dalam campuran", yaxis_title="rasio Cc/Ce",
                      legend=dict(orientation="h", y=-0.3), paper_bgcolor="rgba(0,0,0,0)",
                      plot_bgcolor="rgba(0,0,0,0)")
    return fig

# ================= RUN =================
if run:
    if up is not None and up.size > 20 * 1024 * 1024:
        st.error("File terlalu besar (maks 20 MB). Untuk genom utuh, gunakan pipeline offline."); st.stop()
    with st.spinner("Menjalankan panel 16 marker…"):
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
        call1, cls1, ic1 = "INVALID — DNA kopi tidak terdeteksi", "verdict-flag", "⛔"
    elif ce == 0:
        call1, cls1, ic1 = "FLAG: bukan Arabika (kemungkinan Robusta murni)", "verdict-flag", "🚨"
    elif ratio == ratio and ratio > THRESH:
        call1, cls1, ic1 = f"FLAG: kemungkinan campuran Robusta", "verdict-warn", "⚠️"
    else:
        call1, cls1, ic1 = "ARABICA AUTHENTIC", "verdict-ok", "✅"

    ota_pos = [m for m, n in ota.items() if n > 0]
    if len(ota_pos) >= 2:
        call2, cls2, ic2 = "FLAG: DNA kapang penghasil OTA terdeteksi — kirim lot untuk uji HPLC", "verdict-flag", "🍄"
    elif len(ota_pos) == 1:
        call2, cls2, ic2 = "CHECK: satu marker OTA positif — periksa ulang", "verdict-warn", "🔎"
    else:
        call2, cls2, ic2 = "Tidak ada DNA kapang penghasil OTA", "verdict-ok", "✅"

    st.markdown(f"<div class='note'>🧾 Sampel <b>{label}</b> · {len(seqs):,} fragmen · dianalisis dalam {dt:.1f} detik · panel v0.1</div>",
                unsafe_allow_html=True)

    # ---- Modul 1 ----
    st.markdown("<div class='section-title'><span class='n'>1</span> Autentikasi Arabika</div>", unsafe_allow_html=True)
    c1, c2 = st.columns([1.15, 1])
    with c1:
        m = st.columns(4)
        m[0].markdown(f"<div class='mcard'><div class='v'>{ce}</div><div class='l'>M1-Ce</div><div class='s'>bukti arabika</div></div>", unsafe_allow_html=True)
        m[1].markdown(f"<div class='mcard'><div class='v'>{cc}</div><div class='l'>M1-Cc</div><div class='s'>komponen robusta</div></div>", unsafe_allow_html=True)
        m[2].markdown(f"<div class='mcard'><div class='v'>{ref}</div><div class='l'>REF</div><div class='s'>{ref_n}/{len(ref_rows)} lokus kontrol</div></div>", unsafe_allow_html=True)
        m[3].markdown(f"<div class='mcard'><div class='v'>{ratio:.2f}</div><div class='l'>rasio Cc/Ce</div><div class='s'>ambang {THRESH:.2f}</div></div>" if ratio == ratio else
                      f"<div class='mcard'><div class='v'>—</div><div class='l'>rasio Cc/Ce</div><div class='s'>n/a</div></div>", unsafe_allow_html=True)
        st.markdown(f"""<div class='verdict {cls1}'><div class='icon'>{ic1}</div>
        <div><div class='big'>{call1}</div><div class='small'>rasio Cc/Ce = {f'{ratio:.2f}' if ratio==ratio else '—'} · ambang {THRESH:.2f}</div></div></div>""",
                    unsafe_allow_html=True)
    with c2:
        if sim is not None:
            st.plotly_chart(ratio_gauge(ratio) if ratio == ratio else ratio_gauge(0), use_container_width=True)
            st.plotly_chart(sim_chart(ratio), use_container_width=True)
    st.markdown("<div class='note'>⚠️ Rasio &amp; ambang berasal dari simulasi copy-number yang <b>belum dikalibrasi</b> dengan sampel fisik — bukan klaim batas deteksi.</div>", unsafe_allow_html=True)

    # ---- Modul 2 ----
    st.markdown("<div class='section-title'><span class='n'>2</span> Risiko Kapang OTA</div>", unsafe_allow_html=True)
    cols = st.columns(max(1, len(ota)))
    for c, (mk, n) in zip(cols, ota.items()):
        col_color = "#c62828" if n > 0 else "#2e7d32"
        c.markdown(f"<div class='mcard'><div class='v' style='color:{col_color}'>{'●' if n>0 else '○'} {n}</div><div class='l'>{mk}</div></div>", unsafe_allow_html=True)
    st.markdown(f"""<div class='verdict {cls2}'><div class='icon'>{ic2}</div>
    <div><div class='big'>{call2}</div><div class='small'>deteksi gen biosintesis OTA (otaB/otaD) — skrining risiko, bukan pengukuran toksin</div></div></div>""",
                unsafe_allow_html=True)

    with st.expander("🔬 Rincian per marker"):
        det = pd.DataFrame([{"marker": mk, "tipe": primers.set_index("marker_id").loc[mk, "type"],
                             "amplikon": n} for mk, n in hits.items()])
        st.dataframe(det, hide_index=True, use_container_width=True)

    html = f"""<html><body style="font-family:sans-serif">
    <h2>KopiFied — Laporan QA (v0.1)</h2>
    <p><b>Sampel:</b> {label} | <b>Fragmen:</b> {len(seqs)} | <b>Waktu:</b> {dt:.1f}s</p>
    <h3>Modul 1 — Autentikasi Arabika</h3>
    <p>M1-Ce={ce}, M1-Cc={cc}, REF={ref} ({ref_n}/{len(ref_rows)} lokus), rasio={'%.2f' % ratio if ratio == ratio else '—'}<br><b>{call1}</b></p>
    <h3>Modul 2 — Risiko Kapang OTA</h3>
    <p>{', '.join(f'{m}={n}' for m, n in ota.items())}<br><b>{call2}</b></p>
    <hr><p><i>In-silico proof of concept (TRL 3). Belum dikalibrasi; keputusan akhir memerlukan validasi lab.</i></p>
    </body></html>"""
    st.download_button("⬇️ Unduh Laporan QA (HTML)", html,
                       file_name=f"KopiFied_QA_{label.replace('.fa','')}_v0.1.html")

# ================= TABS BUKTI =================
st.divider()
tabs = st.tabs(["🧪 Matriks Validasi", "🧬 Primer", "🗂️ Genom", "📡 Bukti SRA & Multiplex", "📖 Metode & Batasan"])
with tabs[0]:
    p = f"{RES}/validation_matrix.png"
    if os.path.exists(p):
        st.image(p, caption="110/125 sel sesuai harapan (88%). Sumber: results/validation_matrix.tsv")
with tabs[1]:
    st.dataframe(primers, hide_index=True, use_container_width=True)
with tabs[2]:
    g = f"{BASE}/config/genomes.tsv"
    if os.path.exists(g):
        st.dataframe(pd.read_csv(g, sep="\t"), hide_index=True, use_container_width=True)
with tabs[3]:
    st.markdown("**Validasi sampel independen (SRA SRR40981778 — arabika lini mutan azide)**")
    sp = f"{BASE}/output/sra_check/sra_primer_sites.tsv"
    if os.path.exists(sp):
        st.dataframe(pd.read_csv(sp, sep="\t"), hide_index=True, use_container_width=True)
        st.caption("9/14 marker dengan ≥1 read cocok primer penuh pada ~1,3× coverage; marker OTA nol di reads kopi (sesuai harapan).")
    st.markdown("**Kompatibilitas multiplex (primer3.calc_heterodimer)**")
    mp = f"{RES}/multiplex_check.tsv"
    if os.path.exists(mp):
        mux = pd.read_csv(mp, sep="\t")
        st.write(f"Hanya {len(mux)} dari 1.035 pasangan oligo (0,8%) berisiko dimer — semuanya interaksi probe–probe.")
        st.dataframe(mux, hide_index=True, use_container_width=True)
with tabs[4]:
    r = f"{BASE}/README.md"
    if os.path.exists(r):
        st.markdown(open(r).read())
