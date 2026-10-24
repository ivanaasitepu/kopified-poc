"""
KopiFied — Decision Support Dashboard
Panel Marker DNA untuk Autentikasi Kopi Arabika & Skrining Risiko Kapang OTA
Tim Kimpulita · GenoHack 2026
"""
import os, sys, time, hashlib, html as _html
from datetime import datetime
import streamlit as st
import pandas as pd

st.set_page_config(page_title="KopiFied · Panel Marker DNA", page_icon="☕", layout="wide",
                   initial_sidebar_state="expanded")

BASE = os.environ.get("KOPIFIED_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = f"{BASE}/results"
SAMPLES = f"{BASE}/dashboard/samples"
sys.path.insert(0, f"{BASE}/dashboard")
from insilico import count_amplicons

PANEL_VERSION = "v0.1"
LOGO_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static", "kopified_logo.png")
def logo_b64():
    import base64 as _b64
    if os.path.exists(LOGO_PATH):
        return "data:image/png;base64," + _b64.b64encode(open(LOGO_PATH, "rb").read()).decode()
    return ""
LOGO = logo_b64()
MAX_MM, MIN_AMP, MAX_AMP = 2, 50, 1000

# ---------- Palet (dipakai juga oleh grafik & laporan) ----------
C = dict(ink="#1C1917", ink2="#44403C", muted="#78716C", line="#E7E2D9", bg="#F7F5F0", paper="#FFFFFF",
         espresso="#2B1D16", roast="#5B3A29", copper="#A86B3C", gold="#C9A56A",
         ok="#2F6B4F", ok_bg="#EAF3EE", warn="#A15C00", warn_bg="#FDF3E2",
         flag="#9F2A2A", flag_bg="#F9E9E7", inv="#4B5563", inv_bg="#EEF0F3")

# ================= DESIGN SYSTEM =================
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&family=JetBrains+Mono:wght@500;700&display=swap');

:root{
  --ink:#1C1917; --ink-2:#44403C; --muted:#78716C; --line:#E7E2D9; --bg:#F7F5F0; --paper:#FFFFFF;
  --espresso:#2B1D16; --roast:#5B3A29; --copper:#A86B3C; --gold:#C9A56A;
  --ok:#2F6B4F; --ok-bg:#EAF3EE; --warn:#A15C00; --warn-bg:#FDF3E2;
  --flag:#9F2A2A; --flag-bg:#F9E9E7; --inv:#4B5563; --inv-bg:#EEF0F3;
}
html, body, .stApp, [class*="css"], .stMarkdown, p, li, label { font-family:'Inter',system-ui,sans-serif; color:var(--ink); }
.stApp{ background:var(--bg); }
.block-container{ padding-top:1.6rem; padding-bottom:3rem; max-width:1380px; }
#MainMenu, footer{ visibility:hidden; }
header[data-testid="stHeader"]{ background:transparent; }
.mono{ font-family:'JetBrains Mono',monospace; }

/* ---------- HERO ---------- */
.hero{ position:relative; overflow:hidden; border-radius:20px; padding:30px 36px;
  background: radial-gradient(1200px 400px at 100% 0%, rgba(201,165,106,.22), transparent 60%),
              linear-gradient(135deg, #1E140F 0%, #2B1D16 50%, #4A2F22 100%);
  color:#F8F3EA; display:grid; grid-template-columns: 1.5fr 1fr; gap:28px; align-items:center;
  box-shadow:0 10px 30px rgba(43,29,22,.18); }
.hero .eyebrow{ font-size:11.5px; letter-spacing:2.2px; text-transform:uppercase; color:var(--gold); font-weight:700; }
.hero h1{ font-family:'Source Serif 4',serif; font-weight:700; font-size:38px; margin:6px 0 6px; color:#FFF; letter-spacing:-.3px; }
.hero .sub{ font-size:14.5px; line-height:1.55; color:#E9DFCF; max-width:640px; }
.hero .chips{ display:flex; flex-wrap:wrap; gap:8px; margin-top:16px; }
.chip{ font-size:11.5px; font-weight:600; padding:5px 11px; border-radius:999px;
  background:rgba(255,255,255,.08); border:1px solid rgba(255,255,255,.18); color:#F3EADB; }
.chip.warn{ background:rgba(232,160,60,.16); border-color:rgba(232,160,60,.55); color:#FFD9A0; }
.hero-stats{ display:grid; grid-template-columns:1fr 1fr; gap:10px; }
.hs{ background:rgba(255,255,255,.06); border:1px solid rgba(255,255,255,.12); border-radius:14px; padding:14px 16px; }
.hs .v{ font-family:'JetBrains Mono',monospace; font-size:24px; font-weight:700; color:#FFF; }
.hs .l{ font-size:11px; letter-spacing:.8px; text-transform:uppercase; color:#CDBFA8; margin-top:2px; font-weight:600; }

/* ---------- KARTU ---------- */
.card{ background:var(--paper); border:1px solid var(--line); border-radius:16px; padding:20px 22px;
  box-shadow:0 1px 2px rgba(28,25,23,.04), 0 6px 18px rgba(28,25,23,.04); }
.card h4{ margin:0 0 4px; font-size:13px; color:var(--muted); letter-spacing:.8px; text-transform:uppercase; font-weight:700; }
.kpi{ background:var(--paper); border:1px solid var(--line); border-radius:14px; padding:14px 16px; height:100%; }
.kpi .l{ font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.8px; font-weight:700; }
.kpi .v{ font-family:'JetBrains Mono',monospace; font-size:28px; font-weight:700; color:var(--ink); margin-top:4px; line-height:1.1; }
.kpi .s{ font-size:12px; color:var(--ink-2); margin-top:4px; }
.kpi.accent{ border-top:3px solid var(--copper); }

/* ---------- STATUS ---------- */
.status{ border-radius:16px; padding:18px 20px; border:1px solid; display:flex; gap:16px; align-items:flex-start; height:100%; }
.status .tag{ font-size:10.5px; font-weight:800; letter-spacing:1.2px; padding:4px 9px; border-radius:6px; color:#fff; white-space:nowrap; margin-top:2px; }
.status .t{ font-size:17px; font-weight:700; line-height:1.3; }
.status .d{ font-size:13px; color:var(--ink-2); margin-top:4px; line-height:1.5; }
.status .m{ font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.8px; font-weight:700; margin-bottom:2px; }
.s-ok{ background:var(--ok-bg); border-color:#BFD9CA; } .s-ok .tag{ background:var(--ok); } .s-ok .t{ color:var(--ok); }
.s-warn{ background:var(--warn-bg); border-color:#F0D3A3; } .s-warn .tag{ background:var(--warn); } .s-warn .t{ color:var(--warn); }
.s-flag{ background:var(--flag-bg); border-color:#E8C1BC; } .s-flag .tag{ background:var(--flag); } .s-flag .t{ color:var(--flag); }
.s-invalid{ background:var(--inv-bg); border-color:#D1D5DB; } .s-invalid .tag{ background:var(--inv); } .s-invalid .t{ color:var(--inv); }

.overall{ border-radius:18px; padding:22px 26px; color:#fff; display:flex; justify-content:space-between; align-items:center; gap:20px; }
.overall .k{ font-size:11px; letter-spacing:2px; text-transform:uppercase; opacity:.85; font-weight:700; }
.overall .v{ font-family:'Source Serif 4',serif; font-size:30px; font-weight:700; margin-top:2px; }
.overall .d{ font-size:13.5px; opacity:.95; margin-top:4px; max-width:760px; line-height:1.5; }
.o-ok{ background:linear-gradient(120deg,#22513B,#2F6B4F); }
.o-warn{ background:linear-gradient(120deg,#7C4700,#A15C00); }
.o-flag{ background:linear-gradient(120deg,#7A1F1F,#9F2A2A); }
.o-invalid{ background:linear-gradient(120deg,#374151,#4B5563); }
.overall .meta{ text-align:right; font-family:'JetBrains Mono',monospace; font-size:12px; opacity:.9; line-height:1.7; white-space:nowrap; }

/* ---------- SECTION ---------- */
.sec{ display:flex; align-items:baseline; gap:12px; margin:30px 0 12px; padding-bottom:10px; border-bottom:1px solid var(--line); }
.sec .n{ font-family:'JetBrains Mono',monospace; font-size:12px; font-weight:700; color:var(--copper);
  border:1px solid var(--copper); border-radius:6px; padding:2px 7px; }
.sec .t{ font-family:'Source Serif 4',serif; font-size:23px; font-weight:700; color:var(--ink); }
.sec .s{ font-size:13px; color:var(--muted); margin-left:auto; }

.meta-strip{ display:flex; flex-wrap:wrap; gap:0; background:var(--paper); border:1px solid var(--line); border-radius:14px; overflow:hidden; margin:16px 0; }
.meta-strip > div{ padding:12px 18px; border-right:1px solid var(--line); flex:1 1 150px; }
.meta-strip > div:last-child{ border-right:none; }
.meta-strip .l{ font-size:10.5px; color:var(--muted); text-transform:uppercase; letter-spacing:.9px; font-weight:700; }
.meta-strip .v{ font-size:13.5px; color:var(--ink); font-weight:600; margin-top:3px; word-break:break-all; }

.note{ background:#FBF7EF; border:1px solid #EADFC9; border-left:3px solid var(--gold); border-radius:10px;
  padding:11px 15px; font-size:12.8px; color:#5C4A2E; line-height:1.55; margin-top:10px; }
.rec{ background:var(--paper); border:1px solid var(--line); border-radius:14px; padding:14px 18px; }
.rec li{ font-size:13.5px; color:var(--ink-2); margin:4px 0; line-height:1.5; }

.ota-chip{ background:var(--paper); border:1px solid var(--line); border-radius:14px; padding:14px 16px; display:flex; justify-content:space-between; align-items:center; }
.ota-chip .name{ font-family:'JetBrains Mono',monospace; font-size:13px; font-weight:700; color:var(--ink); }
.ota-chip .sub{ font-size:11.5px; color:var(--muted); margin-top:2px; }
.pill{ font-size:11px; font-weight:800; letter-spacing:.8px; padding:4px 10px; border-radius:999px; }
.pill.pos{ background:var(--flag-bg); color:var(--flag); border:1px solid #E8C1BC; }
.pill.neg{ background:var(--ok-bg); color:var(--ok); border:1px solid #BFD9CA; }

.step{ background:var(--paper); border:1px solid var(--line); border-radius:16px; padding:20px; height:100%; }
.step .n{ font-family:'JetBrains Mono',monospace; color:var(--copper); font-weight:700; font-size:13px; }
.step .t{ font-weight:700; font-size:16px; margin:6px 0 4px; color:var(--ink); }
.step .d{ font-size:13px; color:var(--ink-2); line-height:1.55; }

/* ---------- SIDEBAR ---------- */
section[data-testid="stSidebar"]{ background:#FBFAF7; border-right:1px solid var(--line); }
section[data-testid="stSidebar"] .brand{ font-family:'Source Serif 4',serif; font-size:22px; font-weight:700; color:var(--espresso); }
section[data-testid="stSidebar"] .brand-sub{ font-size:11.5px; color:var(--muted); letter-spacing:.6px; text-transform:uppercase; font-weight:600; margin-bottom:6px; }
section[data-testid="stSidebar"] h5{ font-size:11.5px; letter-spacing:1px; text-transform:uppercase; color:var(--copper); margin:14px 0 4px; font-weight:800; }
section[data-testid="stSidebar"] .stButton button{
  background:var(--espresso); color:#FFF; border:none; border-radius:12px; font-weight:700; letter-spacing:.4px; padding:.7rem 1rem; }
section[data-testid="stSidebar"] .stButton button:hover{ background:var(--roast); color:#FFF; }
section[data-testid="stSidebar"] .stButton button p{ color:#FFF !important; }
.rule{ font-size:12.3px; color:var(--ink-2); line-height:1.55; padding:8px 10px; border-radius:8px; background:#fff; border:1px solid var(--line); margin-bottom:6px; }
.rule b{ color:var(--ink); }

/* ---------- TABS & TOMBOL ---------- */
div[data-testid="stTabs"] button[role="tab"]{ font-weight:600; color:var(--ink-2); }
div[data-testid="stTabs"] button[role="tab"][aria-selected="true"]{ color:var(--espresso); }
.stDownloadButton button{ border-radius:10px; font-weight:700; border:1px solid var(--espresso); color:var(--espresso); background:#fff; }
.stDownloadButton button:hover{ background:var(--espresso); color:#fff; }
@media (max-width: 900px){ .hero{ grid-template-columns:1fr; } .overall{ flex-direction:column; align-items:flex-start; } .overall .meta{ text-align:left; } }
</style>
""", unsafe_allow_html=True)

# ================= DATA =================
PANEL = f"{RES}/primers_v0.1.tsv"
if not os.path.exists(PANEL):
    st.warning("Panel primer belum tersedia."); st.stop()
primers = pd.read_csv(PANEL, sep="\t")

sim, THRESH, R0 = None, 1.6, 1.0
if os.path.exists(f"{RES}/blend_simulation.tsv"):
    sim = pd.read_csv(f"{RES}/blend_simulation.tsv", sep="\t")
    R0 = float(sim.loc[sim.robusta_pct == 0, "ratio_Cc_Ce"].iloc[0])
    r25 = float(sim.loc[sim.robusta_pct == 25, "ratio_Cc_Ce"].iloc[0])
    THRESH = R0 + (r25 - R0) * 0.3
# A1_0008 hanya teramplifikasi di genom ET-39 (0 di Cara_1.0/Cara_2.4) — tetap dilaporkan, tapi tidak dipakai untuk rasio
RATIO_EXCLUDE = {"A1_0008"}
ARA_CAP = 95  # di atas ini selisih rasio terlalu kecil untuk dibedakan — tampilkan sebagai "≥95%"

LABELS = {
    "arabica_caturra.fa": ("Arabika Gayo — sampel pasar", "Arabika murni", "AUTHENTIC"),
    "mold_free.fa": ("Arabika bersih — tanpa kapang", "Arabika murni", "AUTHENTIC · bebas kapang"),
    "blend_10pct.fa": ("Campuran 10% Robusta", "Fraud halus", "FLAG campuran"),
    "blend_25pct.fa": ("Campuran 25% Robusta", "Fraud tipikal", "FLAG campuran"),
    "blend_50pct.fa": ("Campuran 50% Robusta", "Fraud berat", "FLAG campuran"),
    "robusta.fa": ("Robusta murni berlabel 'arabika'", "Substitusi total", "FLAG bukan arabika"),
    "arabica_with_carbonarius.fa": ("Arabika + A. carbonarius", "Kontaminasi kapang OTA+", "AUTHENTIC · FLAG OTA"),
    "mold_ota_positive.fa": ("Kontrol positif — ekstrak kapang", "Bukan matriks kopi", "INVALID · FLAG OTA"),
}
TYPE_DESC = {"M1-Ce": "Subgenom C. eugenioides (penanda Arabika)",
             "M1-Cc": "Subgenom C. canephora (penanda Robusta)",
             "REF": "Kontrol internal — lokus referensi kopi"}

def type_desc(t):
    return TYPE_DESC.get(t, "Gen biosintesis OTA (Aspergillus/Penicillium)" if str(t).startswith("OTA") else t)

def fmt_int(n):
    return f"{n:,}".replace(",", ".")

def esc(s):
    return _html.escape(str(s))

# ================= SIDEBAR =================
with st.sidebar:
    st.markdown((f"<img src='{LOGO}' style='max-width:150px; margin-bottom:4px;'>" if LOGO else "") +
                f"<div class='brand-sub'>Decision Support · Panel {PANEL_VERSION}</div>", unsafe_allow_html=True)
    st.markdown("<h5>Sampel uji</h5>", unsafe_allow_html=True)
    preset = [f for f in sorted(os.listdir(SAMPLES)) if f.endswith(".fa")] if os.path.isdir(SAMPLES) else []
    disp = [LABELS.get(f, (f,))[0] for f in preset]
    idx = st.selectbox("Sampel siap-pakai", range(len(preset)), format_func=lambda i: disp[i],
                       label_visibility="collapsed") if preset else None
    if idx is not None and preset[idx] in LABELS:
        _, cat, exp = LABELS[preset[idx]]
        st.caption(f"{cat} · hasil diharapkan: **{exp}**")
    up = st.file_uploader("Atau unggah file ekspor / FASTA (maks 20 MB)", type=["fa", "fasta", "txt"])
    run = st.button("Jalankan Analisis Panel", use_container_width=True)

    st.markdown("<h5>Aturan keputusan</h5>", unsafe_allow_html=True)
    st.markdown(f"""
    <div class='rule'><b>Invalid</b> — REF, M1-Ce, dan M1-Cc seluruhnya nol.</div>
    <div class='rule'><b>Bukan Arabika</b> — M1-Ce = 0 (penanda subgenom eugenioides absen).</div>
    <div class='rule'><b>Indikasi campuran</b> — rasio Cc/Ce &gt; <span class='mono'>{THRESH:.2f}</span>.</div>
    <div class='rule'><b>Risiko OTA</b> — ≥2 marker OTA positif = FLAG; 1 marker = periksa ulang.</div>
    """, unsafe_allow_html=True)
    st.caption("Ambang rasio diturunkan dari simulasi copy-number (r₀ + 0,3·(r₂₅ − r₀)); belum dikalibrasi lab. "
               f"Marker {', '.join(sorted(RATIO_EXCLUDE))} tidak dipakai dalam rasio (hanya teramplifikasi di sebagian genom Arabika).")

# ================= HELPERS =================
def read_fasta(src):
    seqs, buf = [], []
    lines = open(src) if isinstance(src, str) else src.getvalue().decode(errors="ignore").splitlines()
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if buf: seqs.append("".join(buf).upper())
            buf = []
        elif line: buf.append(line)
    if buf: seqs.append("".join(buf).upper())
    return seqs

def arabica_pct(ce, cc, ratio):
    """Estimasi % Arabika dari model copy-number: Ce ∝ fraksi Arabika, Cc konstan → r = r₀/(1−f), jadi %Arabika = 100·r₀/r."""
    if ce == 0:
        return 0.0 if cc > 0 else None
    if ratio != ratio or ratio <= 0:
        return None
    return max(0.0, min(100.0, 100 * R0 / ratio))

def fmt_pct(p):
    if p is None:
        return "—"
    return f"≥{ARA_CAP}%" if p >= ARA_CAP else f"≈{p:.0f}%"

def classify_auth(ce, cc, ref, ratio):
    if ref == 0 and ce == 0 and cc == 0:
        return dict(code="invalid", tag="INVALID", title="DNA kopi tidak terdeteksi",
                    desc="Kontrol internal dan kedua penanda subgenom bernilai nol — sampel tidak dapat dinilai.")
    if ce == 0:
        return dict(code="flag", tag="FLAG", title="Bukan Arabika — kemungkinan Robusta murni",
                    desc="Penanda subgenom C. eugenioides tidak terdeteksi, padahal Arabika (alotetraploid) selalu membawanya.")
    if ratio == ratio and ratio > THRESH:
        return dict(code="warn", tag="FLAG", title="Indikasi campuran Robusta",
                    desc=f"Rasio Cc/Ce {ratio:.2f} melampaui ambang {THRESH:.2f} — kelebihan sinyal subgenom canephora. "
                         f"Estimasi kandungan Arabika {fmt_pct(arabica_pct(ce, cc, ratio))} (indikatif).")
    return dict(code="ok", tag="PASS", title="Arabika autentik",
                desc=f"Penanda Arabika terdeteksi dan rasio Cc/Ce {ratio:.2f} berada di bawah ambang {THRESH:.2f}. "
                     f"Estimasi kandungan Arabika {fmt_pct(arabica_pct(ce, cc, ratio))} (indikatif).")

def classify_ota(ota):
    pos = [m for m, n in ota.items() if n > 0]
    if len(pos) >= 2:
        return dict(code="flag", tag="FLAG", title="DNA kapang penghasil OTA terdeteksi",
                    desc=f"{len(pos)} dari {len(ota)} marker gen biosintesis OTA positif. Kirim lot untuk kuantifikasi toksin (HPLC).")
    if len(pos) == 1:
        return dict(code="warn", tag="CHECK", title="Satu marker OTA positif",
                    desc=f"Hanya {pos[0]} yang positif — sinyal tunggal, perlu ekstraksi & amplifikasi ulang.")
    return dict(code="ok", tag="PASS", title="Tidak ada DNA kapang penghasil OTA",
                desc=f"Seluruh {len(ota)} marker OTA negatif.")

def overall_call(a, o):
    if a["code"] == "invalid":
        return dict(code="invalid", title="INVALID", desc="Sampel tidak memenuhi syarat analisis. Periksa ekstraksi DNA dan identitas matriks sampel.")
    if a["code"] == "ok" and o["code"] == "ok":
        return dict(code="ok", title="LOLOS", desc="Sampel konsisten dengan Arabika autentik dan tidak menunjukkan risiko kapang penghasil OTA.")
    if "flag" in (a["code"], o["code"]):
        return dict(code="flag", title="TIDAK LOLOS", desc="Ditemukan ketidaksesuaian yang memerlukan tindakan — lihat rekomendasi di bawah.")
    return dict(code="warn", title="PERLU TINDAK LANJUT", desc="Hasil belum konklusif — lakukan pengujian konfirmasi sebelum lot dirilis.")

def recommendations(a, o):
    r = []
    if a["code"] == "invalid":
        r.append("Ulangi ekstraksi DNA; pastikan sampel adalah matriks kopi (biji hijau/sangrai/bubuk).")
    elif a["code"] in ("flag", "warn"):
        r.append("Tahan lot; lakukan uji konfirmasi autentisitas (replikasi qPCR dan/atau analisis kimia 16-O-methylcafestol).")
        r.append("Telusuri rantai pasok lot terkait (pemasok, tanggal panen, gudang).")
    if o["code"] == "flag":
        r.append("Kirim lot untuk kuantifikasi okratoksin A dengan HPLC-FLD sebelum keputusan rilis.")
        r.append("Evaluasi kondisi pengeringan dan penyimpanan (kadar air, aw) pada lot terdampak.")
    elif o["code"] == "warn":
        r.append("Ulangi ekstraksi & amplifikasi marker OTA dalam duplikat untuk konfirmasi sinyal tunggal.")
    if not r:
        r.append("Tidak diperlukan tindakan tambahan; arsipkan laporan bersama catatan lot.")
    return r

PLOT_FONT = dict(family="Inter, sans-serif", color=C["ink"], size=12)

def ratio_gauge(ratio):
    import plotly.graph_objects as go
    top = max(2.5, (ratio if ratio == ratio else 0) * 1.15)
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=ratio if ratio == ratio else 0,
        number={"font": {"size": 40, "family": "JetBrains Mono", "color": C["ink"]}, "valueformat": ".2f"},
        title={"text": "Rasio Cc/Ce", "font": {"size": 13, "color": C["muted"]}},
        gauge={"axis": {"range": [0, top], "tickwidth": 1, "tickcolor": C["line"], "tickfont": {"color": C["muted"], "size": 11}},
               "bar": {"color": C["espresso"], "thickness": 0.28},
               "bgcolor": "white", "borderwidth": 0,
               "steps": [{"range": [0, THRESH], "color": "#DCEBE2"},
                         {"range": [THRESH, top], "color": "#F2D6D2"}],
               "threshold": {"line": {"color": C["flag"], "width": 3}, "thickness": 0.9, "value": THRESH}}))
    fig.update_layout(height=230, margin=dict(l=24, r=24, t=46, b=8), paper_bgcolor="rgba(0,0,0,0)", font=PLOT_FONT)
    return fig

def sim_position(ratio):
    """Posisi indikatif sampel pada kurva simulasi (interpolasi linear; hanya jika kurva monoton naik)."""
    if sim is None or ratio != ratio:
        return None
    s = sim.sort_values("robusta_pct")
    y, x = s["ratio_Cc_Ce"].tolist(), s["robusta_pct"].tolist()
    if any(b < a for a, b in zip(y, y[1:])) or ratio < y[0] or ratio > y[-1]:
        return None
    for i in range(len(y) - 1):
        if y[i] <= ratio <= y[i + 1]:
            return x[i] if y[i + 1] == y[i] else x[i] + (ratio - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    return None

def sim_chart(ratio):
    import plotly.graph_objects as go
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sim["robusta_pct"], y=sim["ratio_Cc_Ce"], mode="lines+markers",
                             name="Kurva simulasi", line=dict(color=C["roast"], width=2.6),
                             marker=dict(size=7, color="white", line=dict(color=C["roast"], width=2))))
    fig.add_hline(y=THRESH, line_dash="dash", line_color=C["flag"], line_width=1.5,
                  annotation_text=f"ambang {THRESH:.2f}", annotation_position="top left",
                  annotation_font=dict(color=C["flag"], size=11))
    if ratio == ratio:
        fig.add_hline(y=ratio, line_dash="dot", line_color=C["copper"], line_width=1.5)
        xp = sim_position(ratio)
        if xp is not None:
            fig.add_trace(go.Scatter(x=[xp], y=[ratio], mode="markers", name=f"Sampel (≈{xp:.0f}% indikatif)",
                                     marker=dict(size=15, color=C["copper"], symbol="diamond",
                                                 line=dict(color="white", width=2))))
        else:
            fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines", name=f"Rasio sampel {ratio:.2f}",
                                     line=dict(color=C["copper"], dash="dot")))
    fig.update_layout(height=270, margin=dict(l=8, r=8, t=24, b=8), font=PLOT_FONT,
                      xaxis=dict(title="% Robusta dalam campuran", gridcolor=C["line"], zeroline=False, ticksuffix="%"),
                      yaxis=dict(title="Rasio Cc/Ce", gridcolor=C["line"], zeroline=False),
                      legend=dict(orientation="h", y=-0.32, font=dict(size=11)),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="white")
    return fig

def marker_chart(det):
    import plotly.graph_objects as go
    col = {"M1-Ce": C["ok"], "M1-Cc": C["copper"], "REF": "#6B7280"}
    colors = [col.get(t, C["flag"]) for t in det["Tipe"]]
    fig = go.Figure(go.Bar(x=det["Amplikon"], y=det["Marker"], orientation="h", marker_color=colors,
                           text=det["Amplikon"], textposition="outside", textfont=dict(family="JetBrains Mono", size=11)))
    fig.update_layout(height=max(260, 26 * len(det) + 60), margin=dict(l=8, r=30, t=10, b=8), font=PLOT_FONT,
                      xaxis=dict(title="Jumlah amplikon in-silico", gridcolor=C["line"]),
                      yaxis=dict(autorange="reversed", tickfont=dict(family="JetBrains Mono", size=11)),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="white", showlegend=False)
    return fig

# ================= LAPORAN HTML =================
def svg_ratio_scale(ratio):
    W, H, pad = 680, 92, 30
    top = max(2.5, (ratio if ratio == ratio else 0) * 1.15)
    X = lambda v: pad + (W - 2 * pad) * v / top
    g = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" font-family="Arial,Helvetica,sans-serif">',
         f'<rect x="{pad}" y="34" width="{X(THRESH)-pad:.1f}" height="16" rx="3" fill="#DCEBE2"/>',
         f'<rect x="{X(THRESH):.1f}" y="34" width="{W-pad-X(THRESH):.1f}" height="16" rx="3" fill="#F2D6D2"/>',
         f'<line x1="{X(THRESH):.1f}" y1="26" x2="{X(THRESH):.1f}" y2="58" stroke="{C["flag"]}" stroke-width="2"/>',
         f'<text x="{X(THRESH):.1f}" y="20" font-size="10.5" fill="{C["flag"]}" text-anchor="middle" font-weight="700">ambang {THRESH:.2f}</text>',
         ]
    zone_labels = (
        f'<text x="{(pad+X(THRESH))/2:.1f}" y="46" font-size="10" fill="{C["ok"]}" text-anchor="middle" font-weight="700" stroke="#DCEBE2" stroke-width="4" paint-order="stroke">KONSISTEN ARABIKA</text>'
        f'<text x="{(X(THRESH)+W-pad)/2:.1f}" y="46" font-size="10" fill="{C["flag"]}" text-anchor="middle" font-weight="700" stroke="#F2D6D2" stroke-width="4" paint-order="stroke">INDIKASI CAMPURAN</text>')
    step = 0.5
    v = 0.0
    while v <= top + 1e-9:
        g.append(f'<line x1="{X(v):.1f}" y1="52" x2="{X(v):.1f}" y2="58" stroke="#A8A29E"/>'
                 f'<text x="{X(v):.1f}" y="72" font-size="10" fill="{C["muted"]}" text-anchor="middle">{v:.1f}</text>')
        v += step
    if ratio == ratio:
        x = X(ratio)
        g.append(f'<path d="M{x-7:.1f},{88} L{x+7:.1f},{88} L{x:.1f},{77} Z" fill="{C["espresso"]}"/>'
                 f'<line x1="{x:.1f}" y1="30" x2="{x:.1f}" y2="54" stroke="{C["espresso"]}" stroke-width="3"/>')
    g.append(zone_labels)
    g.append("</svg>")
    return "".join(g)

def svg_sim_curve(ratio):
    if sim is None:
        return ""
    s = sim.sort_values("robusta_pct")
    xs, ys = s["robusta_pct"].tolist(), s["ratio_Cc_Ce"].tolist()
    W, H, L, R, T, B = 680, 250, 54, 18, 16, 40
    xmax = max(xs) or 1
    ymax = max(max(ys), THRESH, ratio if ratio == ratio else 0) * 1.12
    X = lambda v: L + (W - L - R) * v / xmax
    Y = lambda v: T + (H - T - B) * (1 - v / ymax)
    g = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" font-family="Arial,Helvetica,sans-serif">']
    for k in range(6):
        yv = ymax * k / 5
        g.append(f'<line x1="{L}" y1="{Y(yv):.1f}" x2="{W-R}" y2="{Y(yv):.1f}" stroke="#EEE9E0"/>'
                 f'<text x="{L-8}" y="{Y(yv)+3:.1f}" font-size="10" fill="{C["muted"]}" text-anchor="end">{yv:.2f}</text>')
    for xv in xs:
        g.append(f'<text x="{X(xv):.1f}" y="{H-B+16}" font-size="10" fill="{C["muted"]}" text-anchor="middle">{xv:g}%</text>')
    g.append(f'<line x1="{L}" y1="{H-B}" x2="{W-R}" y2="{H-B}" stroke="#A8A29E"/>')
    g.append(f'<text x="{(L+W-R)/2}" y="{H-4}" font-size="10.5" fill="{C["ink2"]}" text-anchor="middle">% Robusta dalam campuran (simulasi)</text>')
    g.append(f'<text x="12" y="{(T+H-B)/2}" font-size="10.5" fill="{C["ink2"]}" text-anchor="middle" transform="rotate(-90 12 {(T+H-B)/2})">Rasio Cc/Ce</text>')
    g.append(f'<line x1="{L}" y1="{Y(THRESH):.1f}" x2="{W-R}" y2="{Y(THRESH):.1f}" stroke="{C["flag"]}" stroke-dasharray="6 4" stroke-width="1.5"/>'
             f'<text x="{W-R-4}" y="{Y(THRESH)-5:.1f}" font-size="10" fill="{C["flag"]}" text-anchor="end" font-weight="700">ambang {THRESH:.2f}</text>')
    pts = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in zip(xs, ys))
    g.append(f'<polyline points="{pts}" fill="none" stroke="{C["roast"]}" stroke-width="2.4"/>')
    g += [f'<circle cx="{X(a):.1f}" cy="{Y(b):.1f}" r="3.6" fill="#fff" stroke="{C["roast"]}" stroke-width="2"/>' for a, b in zip(xs, ys)]
    if ratio == ratio:
        g.append(f'<line x1="{L}" y1="{Y(ratio):.1f}" x2="{W-R}" y2="{Y(ratio):.1f}" stroke="{C["copper"]}" stroke-dasharray="2 3" stroke-width="1.5"/>')
        xp = sim_position(ratio)
        if xp is not None:
            g.append(f'<rect x="{X(xp)-6:.1f}" y="{Y(ratio)-6:.1f}" width="12" height="12" fill="{C["copper"]}" stroke="#fff" stroke-width="2" transform="rotate(45 {X(xp):.1f} {Y(ratio):.1f})"/>')
        g.append(f'<text x="{L+6}" y="{Y(ratio)-6:.1f}" font-size="10" fill="{C["copper"]}" font-weight="700">sampel {ratio:.2f}</text>')
    g.append("</svg>")
    return "".join(g)

def build_report(r):
    a, o, ov = r["auth"], r["ota_call"], r["overall"]
    col = {"ok": C["ok"], "warn": C["warn"], "flag": C["flag"], "invalid": C["inv"]}
    bg = {"ok": C["ok_bg"], "warn": C["warn_bg"], "flag": C["flag_bg"], "invalid": C["inv_bg"]}
    ratio = r["ratio"]
    rs = f"{ratio:.2f}" if ratio == ratio else "—"
    xp = sim_position(ratio)
    rows = "".join(
        f"<tr><td class='mono'>{esc(d['Marker'])}</td><td>{esc(d['Tipe'])}</td><td>{esc(d['Fungsi'])}</td>"
        f"<td class='num mono'>{d['Amplikon']}</td>"
        f"<td><span class='pill {'pos' if d['Amplikon'] > 0 else 'neg'}'>{'Terdeteksi' if d['Amplikon'] > 0 else 'Tidak terdeteksi'}</span></td></tr>"
        for d in r["det"].to_dict("records"))
    ota_rows = "".join(
        f"<tr><td class='mono'>{esc(m)}</td><td class='num mono'>{n}</td>"
        f"<td><span class='pill {'pos' if n > 0 else 'neg'}'>{'Positif' if n > 0 else 'Negatif'}</span></td></tr>"
        for m, n in r["ota"].items())
    recs = "".join(f"<li>{esc(x)}</li>" for x in r["recs"])
    qc = ("Memadai" if r["ref_n"] == r["ref_total"] and r["ref_total"] else
          "Parsial" if r["ref_n"] > 0 else "Tidak terdeteksi")

    def box(title, c):
        return (f"<div class='res' style='border-color:{col[c['code']]};background:{bg[c['code']]}'>"
                f"<div class='res-k'>{title}</div>"
                f"<div class='res-row'><span class='tag' style='background:{col[c['code']]}'>{c['tag']}</span>"
                f"<span class='res-t' style='color:{col[c['code']]}'>{esc(c['title'])}</span></div>"
                f"<div class='res-d'>{esc(c['desc'])}</div></div>")

    return f"""<!doctype html>
<html lang="id"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Laporan QA {esc(r['report_id'])} — KopiFied</title>

<style>
  @page {{ size: A4; margin: 14mm 13mm 16mm; }}
  * {{ box-sizing: border-box; -webkit-print-color-adjust: exact; print-color-adjust: exact; }}
  body {{ margin:0; background:#EDE9E2; font-family:Arial, Helvetica, sans-serif; color:{C['ink']}; font-size:12.5px; line-height:1.5; }}
  .page {{ max-width: 860px; margin: 28px auto; background:#fff; box-shadow:0 8px 30px rgba(0,0,0,.08); }}
  .mono {{ font-family:Consolas, 'Courier New', monospace; }}
  header {{ background:linear-gradient(135deg,#1E140F,#2B1D16 55%,#4A2F22); color:#F8F3EA; padding:26px 34px; display:flex; justify-content:space-between; gap:20px; }}
  .brand {{ font-family:Arial, Helvetica, sans-serif; font-size:26px; font-weight:700; color:#fff; }}
  .brand small {{ display:block; font-family:Arial, Helvetica, sans-serif; font-size:10.5px; letter-spacing:2px; text-transform:uppercase; color:{C['gold']}; font-weight:700; margin-top:2px; }}
  h1.doc {{ font-family:Arial, Helvetica, sans-serif; font-size:19px; margin:12px 0 0; color:#fff; font-weight:600; }}
  .hmeta {{ text-align:right; font-size:11px; color:#DCCFBB; line-height:1.8; }}
  .hmeta b {{ color:#fff; font-family:Consolas, 'Courier New', monospace; }}
  main {{ padding: 26px 34px 10px; }}
  .verdict {{ display:flex; justify-content:space-between; align-items:center; border-radius:12px; padding:18px 22px; color:#fff; background:{col[ov['code']]}; }}
  .verdict .k {{ font-size:10px; letter-spacing:2px; text-transform:uppercase; opacity:.85; font-weight:700; }}
  .verdict .v {{ font-family:Arial, Helvetica, sans-serif; font-size:26px; font-weight:700; }}
  .verdict .d {{ font-size:12px; opacity:.95; max-width:520px; }}
  .stamp {{ border:2px solid rgba(255,255,255,.8); border-radius:8px; padding:6px 12px; font-weight:800; letter-spacing:2px; font-size:13px; transform:rotate(-4deg); white-space:nowrap; }}
  h2 {{ font-family:Arial, Helvetica, sans-serif; font-size:16.5px; margin:26px 0 10px; padding-bottom:6px; border-bottom:1.5px solid {C['line']}; display:flex; gap:10px; align-items:baseline; }}
  h2 .n {{ font-family:Consolas, 'Courier New', monospace; font-size:11px; color:{C['copper']}; border:1px solid {C['copper']}; border-radius:5px; padding:1px 6px; }}
  table {{ width:100%; border-collapse:collapse; }}
  .kv td {{ padding:6px 10px; border-bottom:1px solid #F1EDE6; vertical-align:top; }}
  .kv td:first-child {{ width:34%; color:{C['muted']}; font-weight:600; font-size:11.5px; }}
  .data th {{ text-align:left; font-size:10.5px; letter-spacing:.8px; text-transform:uppercase; color:{C['muted']}; background:#FAF8F4; padding:8px 10px; border-bottom:1.5px solid {C['line']}; }}
  .data td {{ padding:7px 10px; border-bottom:1px solid #F1EDE6; }}
  .data tr:nth-child(even) td {{ background:#FDFCFA; }}
  .num {{ text-align:right; }}
  .grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:12px; margin-top:12px; }}
  .res {{ border:1px solid; border-left-width:4px; border-radius:10px; padding:12px 14px; }}
  .res-k {{ font-size:10px; letter-spacing:1.4px; text-transform:uppercase; color:{C['muted']}; font-weight:800; }}
  .res-row {{ display:flex; gap:8px; align-items:center; margin-top:4px; }}
  .tag {{ color:#fff; font-size:9.5px; font-weight:800; letter-spacing:1px; padding:3px 7px; border-radius:4px; }}
  .res-t {{ font-weight:700; font-size:13.5px; }}
  .res-d {{ font-size:11.5px; color:{C['ink2']}; margin-top:4px; }}
  .kpis {{ display:grid; grid-template-columns:repeat(5,1fr); gap:10px; margin:6px 0 12px; }}
  .kpi {{ border:1px solid {C['line']}; border-radius:10px; padding:10px 12px; }}
  .kpi .l {{ font-size:9.5px; letter-spacing:1px; text-transform:uppercase; color:{C['muted']}; font-weight:800; }}
  .kpi .v {{ font-family:Consolas, 'Courier New', monospace; font-size:20px; font-weight:700; }}
  .kpi .s {{ font-size:10.5px; color:{C['ink2']}; }}
  .fig {{ border:1px solid {C['line']}; border-radius:10px; padding:12px 14px; margin-top:10px; }}
  .cap {{ font-size:10.5px; color:{C['muted']}; margin-top:6px; }}
  .cap b {{ color:{C['ink2']}; }}
  .pill {{ font-size:10px; font-weight:700; padding:2px 8px; border-radius:999px; white-space:nowrap; }}
  .pill.pos {{ background:{C['flag_bg']}; color:{C['flag']}; }}
  .pill.neg {{ background:{C['ok_bg']}; color:{C['ok']}; }}
  ul.recs {{ margin:6px 0 0; padding-left:18px; }} ul.recs li {{ margin:3px 0; }}
  .rules li {{ margin:3px 0; }}
  .disc {{ background:#FBF7EF; border:1px solid #EADFC9; border-left:3px solid {C['gold']}; border-radius:8px; padding:10px 14px; font-size:11.2px; color:#5C4A2E; }}
  .sign {{ display:grid; grid-template-columns:1fr 1fr 1fr; gap:22px; margin-top:16px; }}
  .sign div {{ border-top:1px solid {C['ink2']}; padding-top:6px; font-size:10.5px; color:{C['muted']}; margin-top:44px; }}
  footer {{ padding:14px 34px 22px; font-size:10px; color:{C['muted']}; display:flex; justify-content:space-between; border-top:1px solid {C['line']}; margin-top:18px; }}
  .toolbar {{ max-width:860px; margin:16px auto -12px; text-align:right; }}
  .toolbar button {{ background:#2B1D16; color:#fff; border:none; border-radius:8px; padding:8px 16px; font-weight:700; cursor:pointer; font-family:Arial, Helvetica, sans-serif; }}
  section {{ break-inside: avoid; }}
  @media print {{ body {{ background:#fff; }} .page {{ margin:0; box-shadow:none; max-width:none; }} .toolbar {{ display:none; }} }}
</style></head>
<body>
<div class="toolbar"><button onclick="window.print()">Cetak / Simpan PDF</button></div>
<div class="page">
<header>
  <div><div class="brand">KopiFied<small>Panel Marker DNA · Tim Kimpulita</small></div>
  <h1 class="doc">Laporan Analisis Autentikasi Arabika &amp; Skrining Risiko Kapang OTA</h1></div>
  <div class="hmeta">No. laporan <b>{esc(r['report_id'])}</b><br>Tanggal <b>{esc(r['ts'])}</b><br>Panel <b>{PANEL_VERSION}</b> · TRL 3</div>
</header>
<main>
  <div class="verdict"><div><div class="k">Kesimpulan keseluruhan</div><div class="v">{esc(ov['title'])}</div><div class="d">{esc(ov['desc'])}</div></div>
  <div class="stamp">{esc(ov['title'])}</div></div>
  <div class="grid2">{box('Modul 1 · Autentikasi Arabika', a)}{box('Modul 2 · Risiko kapang OTA', o)}</div>

  <section><h2><span class="n">01</span>Informasi sampel &amp; parameter analisis</h2>
  <table class="kv">
    <tr><td>Identitas sampel</td><td class="mono">{esc(r['label'])}</td></tr>
    <tr><td>Deskripsi</td><td>{esc(r['desc'])}</td></tr>
    <tr><td>Jumlah fragmen / sekuens</td><td class="mono">{fmt_int(r['n_seq'])} ({fmt_int(r['n_bp'])} bp total)</td></tr>
    <tr><td>Sidik jari data (SHA-256)</td><td class="mono">{esc(r['sha'])}</td></tr>
    <tr><td>Parameter in-silico PCR</td><td>Anchor 3′ 12-mer · maks. {MAX_MM} mismatch per primer · 3 basa ujung 3′ wajib cocok · amplikon {MIN_AMP}–{MAX_AMP} bp</td></tr>
    <tr><td>Waktu analisis</td><td class="mono">{r['dt']:.2f} detik</td></tr>
    {f"<tr><td>Estimasi indikatif komposisi</td><td>≈ <b>{100-xp:.0f}% Arabika</b> / {xp:.0f}% Robusta <i>(interpolasi kurva simulasi; belum dikalibrasi)</i></td></tr>" if xp is not None else ""}
  </table></section>

  <section><h2><span class="n">02</span>Modul 1 — Autentikasi Arabika</h2>
  <div class="kpis">
    <div class="kpi"><div class="l">Estimasi Arabika</div><div class="v">{fmt_pct(r.get('ara_pct'))}</div><div class="s">100·r₀/rasio · indikatif</div></div>
    <div class="kpi"><div class="l">M1-Ce</div><div class="v">{r['ce']}</div><div class="s">bukti subgenom Arabika</div></div>
    <div class="kpi"><div class="l">M1-Cc</div><div class="v">{r['cc']}</div><div class="s">sinyal subgenom canephora</div></div>
    <div class="kpi"><div class="l">REF</div><div class="v">{r['ref']}</div><div class="s">{r['ref_n']}/{r['ref_total']} lokus kontrol</div></div>
    <div class="kpi"><div class="l">Rasio Cc/Ce</div><div class="v">{rs}</div><div class="s">ambang {THRESH:.2f}</div></div>
  </div>
  <div class="fig">{svg_ratio_scale(ratio)}
  <div class="cap"><b>Gambar 1.</b> Posisi rasio Cc/Ce sampel (penanda ▲) terhadap ambang keputusan. Rasio dinormalisasi per jumlah marker: (ΣCc/{r['n_cc']}) ÷ (ΣCe/{r['n_ce']}){'; ' + ', '.join(sorted(RATIO_EXCLUDE)) + ' tidak dihitung (tidak konsisten antar genom Arabika)' if RATIO_EXCLUDE else ''}.</div></div>
  {f'<div class="fig">{svg_sim_curve(ratio)}<div class="cap"><b>Gambar 2.</b> Kurva simulasi copy-number rasio Cc/Ce terhadap proporsi Robusta.' + (f' Posisi indikatif sampel ≈ <b>{xp:.0f}%</b> Robusta (interpolasi; bukan estimasi kuantitatif).' if xp is not None else '') + '</div></div>' if sim is not None else ''}
  </section>

  <section><h2><span class="n">03</span>Modul 2 — Skrining risiko kapang OTA</h2>
  <table class="data"><thead><tr><th>Marker gen biosintesis OTA</th><th class="num">Amplikon</th><th>Status</th></tr></thead><tbody>{ota_rows}</tbody></table>
  <div class="cap">Deteksi DNA gen biosintesis okratoksin A merupakan skrining <b>risiko</b>, bukan pengukuran kadar toksin.</div></section>

  <section><h2><span class="n">04</span>Kontrol kualitas</h2>
  <table class="kv">
    <tr><td>Kontrol internal (REF)</td><td>{r['ref_n']} dari {r['ref_total']} lokus terdeteksi — <b>{qc}</b></td></tr>
    <tr><td>Validitas matriks kopi</td><td>{'Tidak valid — DNA kopi tidak terdeteksi' if a['code']=='invalid' else 'Valid — DNA kopi terdeteksi'}</td></tr>
  </table></section>

  <section><h2><span class="n">05</span>Rincian hasil per marker</h2>
  <table class="data"><thead><tr><th>Marker</th><th>Tipe</th><th>Fungsi</th><th class="num">Amplikon</th><th>Status</th></tr></thead><tbody>{rows}</tbody></table></section>

  <section><h2><span class="n">06</span>Rekomendasi tindak lanjut</h2><ul class="recs">{recs}</ul></section>

  <section><h2><span class="n">07</span>Kriteria keputusan</h2>
  <ul class="rules">
    <li><b>INVALID</b> — REF, M1-Ce, dan M1-Cc seluruhnya nol.</li>
    <li><b>Bukan Arabika</b> — M1-Ce = 0.</li>
    <li><b>Indikasi campuran Robusta</b> — rasio Cc/Ce &gt; {THRESH:.2f} (r₀ + 0,3·(r₂₅ − r₀) dari kurva simulasi). Marker {', '.join(sorted(RATIO_EXCLUDE))} tidak dipakai dalam rasio.</li>
    <li><b>Risiko OTA</b> — ≥2 marker OTA positif = FLAG; tepat 1 marker = CHECK (periksa ulang).</li>
    <li><b>Estimasi % Arabika</b> — 100 · r₀ / rasio Cc/Ce (r₀ = {R0:.2f}, rasio Arabika murni pada kurva simulasi); nilai ≥{ARA_CAP}% dilaporkan sebagai "≥{ARA_CAP}%". Indikatif, per genom (bukan % berat).</li>
  </ul></section>

  <section><h2><span class="n">08</span>Batasan</h2>
  <div class="disc">Hasil ini merupakan <b>in-silico proof of concept (TRL 3)</b>. Rasio, ambang, dan estimasi % Arabika berasal dari simulasi copy-number yang <b>belum dikalibrasi</b> dengan sampel fisik dan tidak merupakan klaim batas deteksi. Keputusan rilis/penolakan lot memerlukan validasi laboratorium (qPCR terkalibrasi, HPLC untuk OTA).</div>
  <div class="sign"><div>Dianalisis oleh</div><div>Ditinjau oleh</div><div>Tanggal &amp; paraf</div></div></section>
</main>
<footer><span>KopiFied {PANEL_VERSION} · Tim Kimpulita · GenoHack 2026</span><span class="mono">{esc(r['report_id'])}</span></footer>
</div></body></html>"""

# ================= HERO =================
n_primers = len(primers)
st.markdown(f"""
<div class="hero">
  <div>
    <div class="eyebrow">Decision Support · DNA-based Coffee QA</div>
    <h1>KopiFied</h1>
    <div class="sub">Panel marker DNA untuk autentikasi kopi Arabika dan skrining risiko kapang penghasil okratoksin A (OTA) —
    menerjemahkan data ekspor amplifikasi menjadi keputusan QA yang dapat ditelusuri.</div>
    <div class="chips">
      <span class="chip">Panel {PANEL_VERSION}</span><span class="chip">Tim Kimpulita · GenoHack 2026</span>
      <span class="chip warn">In-silico PoC — belum dikalibrasi lab</span>
    </div>
  </div>
  <div class="hero-stats"><div class="hs" style="grid-column:1/3; background:#fff; display:flex; align-items:center; justify-content:center; padding:10px 14px;"><img src='{LOGO}' style='max-width:100%; height:56px; object-fit:contain;'></div>
    <div class="hs"><div class="v">{n_primers}</div><div class="l">Primer dalam panel</div></div>
    <div class="hs"><div class="v">25</div><div class="l">Genom publik tervalidasi</div></div>
    <div class="hs"><div class="v">2</div><div class="l">Modul uji sekaligus</div></div>
    <div class="hs"><div class="v">&lt;1 dtk</div><div class="l">Waktu analisis per sampel</div></div>
  </div>
</div>
""", unsafe_allow_html=True)

# ================= RUN =================
if run:
    if up is not None and up.size > 20 * 1024 * 1024:
        st.error("File terlalu besar (maks 20 MB). Untuk genom utuh, gunakan pipeline offline."); st.stop()
    if up is None and idx is None:
        st.error("Tidak ada sampel yang dipilih."); st.stop()
    with st.spinner(f"Menjalankan panel {n_primers} marker…"):
        t0 = time.time()
        src = up if up else f"{SAMPLES}/{preset[idx]}"
        label = getattr(up, "name", None) or preset[idx]
        seqs = read_fasta(src)
        hits = {p.marker_id: count_amplicons(seqs, p.fwd, p.rev, max_mm=MAX_MM) for _, p in primers.iterrows()}
        dt = time.time() - t0

    ce_ids = [m for m in primers.marker_id[primers.type == "M1-Ce"] if m not in RATIO_EXCLUDE]
    cc_ids = [m for m in primers.marker_id[primers.type == "M1-Cc"] if m not in RATIO_EXCLUDE]
    ce = sum(hits.get(m, 0) for m in ce_ids)
    cc = sum(hits.get(m, 0) for m in cc_ids)
    ref_rows = primers[primers.type == "REF"]
    ref = sum(hits.get(m, 0) for m in ref_rows.marker_id)
    ref_n = sum(1 for m in ref_rows.marker_id if hits.get(m, 0) > 0)
    ota = {m: hits.get(m, 0) for m in primers.marker_id[primers.type.str.startswith("OTA")]}
    n_ce = max(1, len(ce_ids)); n_cc = max(1, len(cc_ids))
    ratio = (cc / n_cc) / (ce / n_ce) if ce else float("nan")

    auth = classify_auth(ce, cc, ref, ratio)
    ota_call = classify_ota(ota)
    tmap = primers.set_index("marker_id")["type"]
    det = pd.DataFrame([{"Marker": mk, "Tipe": tmap.loc[mk], "Fungsi": type_desc(tmap.loc[mk]), "Amplikon": int(n)}
                        for mk, n in hits.items()])
    h = hashlib.sha256()
    for s in seqs: h.update(s.encode())
    now = datetime.now()
    st.session_state["res"] = dict(
        label=label, desc=LABELS.get(label, ("Unggahan pengguna", "—", "—"))[0] if not up else "Unggahan pengguna",
        n_seq=len(seqs), n_bp=sum(map(len, seqs)), dt=dt, sha=h.hexdigest()[:32],
        ce=ce, cc=cc, ref=ref, ref_n=ref_n, ref_total=len(ref_rows), ota=ota, ratio=ratio, ara_pct=arabica_pct(ce, cc, ratio),
        n_ce=n_ce, n_cc=n_cc, auth=auth, ota_call=ota_call, overall=overall_call(auth, ota_call),
        recs=recommendations(auth, ota_call), det=det,
        ts=now.strftime("%d-%m-%Y %H:%M:%S"), report_id=f"KPF-{now:%Y%m%d-%H%M%S}")

res = st.session_state.get("res")

if res is None:
    # ---------- Empty state ----------
    st.markdown("<div class='sec'><span class='n'>00</span><span class='t'>Mulai analisis</span>"
                "<span class='s'>Pilih sampel di panel kiri lalu jalankan</span></div>", unsafe_allow_html=True)
    s1, s2, s3 = st.columns(3)
    for col, n, t, d in [
        (s1, "01", "Muat data sampel", "Pilih sampel demo atau unggah file ekspor/FASTA hasil sekuensing atau amplifikasi."),
        (s2, "02", "Jalankan panel 16 marker", "In-silico PCR menghitung amplikon untuk penanda subgenom, kontrol internal, dan gen OTA."),
        (s3, "03", "Tinjau keputusan & laporan", "Dapatkan status LOLOS / TIDAK LOLOS beserta laporan QA siap cetak yang dapat ditelusuri.")]:
        col.markdown(f"<div class='step'><div class='n'>{n}</div><div class='t'>{t}</div><div class='d'>{d}</div></div>",
                     unsafe_allow_html=True)
else:
    r = res
    a, o, ov = r["auth"], r["ota_call"], r["overall"]
    ratio = r["ratio"]

    # ---------- Ringkasan ----------
    st.markdown(f"""
    <div class='meta-strip'>
      <div><div class='l'>Sampel</div><div class='v mono'>{esc(r['label'])}</div></div>
      <div><div class='l'>Deskripsi</div><div class='v'>{esc(r['desc'])}</div></div>
      <div><div class='l'>Fragmen</div><div class='v mono'>{fmt_int(r['n_seq'])}</div></div>
      <div><div class='l'>Waktu analisis</div><div class='v mono'>{r['dt']:.2f} s</div></div>
      <div><div class='l'>No. laporan</div><div class='v mono'>{r['report_id']}</div></div>
    </div>""", unsafe_allow_html=True)

    st.markdown(f"""
    <div class='overall o-{ov['code']}'>
      <div><div class='k'>Kesimpulan keseluruhan</div><div class='v'>{ov['title']}</div><div class='d'>{ov['desc']}</div></div>
      <div class='meta'>Panel {PANEL_VERSION}<br>{r['ts']}<br>SHA {r['sha'][:12]}…</div>
    </div>""", unsafe_allow_html=True)

    st.write("")
    c1, c2 = st.columns(2)
    for col, mod, c in [(c1, "Modul 1 · Autentikasi Arabika", a), (c2, "Modul 2 · Risiko kapang OTA", o)]:
        col.markdown(f"<div class='status s-{c['code']}'><span class='tag'>{c['tag']}</span><div>"
                     f"<div class='m'>{mod}</div><div class='t'>{c['title']}</div><div class='d'>{c['desc']}</div></div></div>",
                     unsafe_allow_html=True)

    # ---------- Modul 1 ----------
    st.markdown("<div class='sec'><span class='n'>01</span><span class='t'>Autentikasi Arabika</span>"
                "<span class='s'>Rasio penanda subgenom C. canephora / C. eugenioides</span></div>", unsafe_allow_html=True)
    k = st.columns(5)
    for col, l, v, s in [(k[0], "Estimasi Arabika", fmt_pct(r.get("ara_pct")), "100·r₀/rasio · indikatif"),
                         (k[1], "M1-Ce", r["ce"], "bukti subgenom Arabika"),
                         (k[2], "M1-Cc", r["cc"], "sinyal subgenom canephora"),
                         (k[3], "REF", r["ref"], f"{r['ref_n']}/{r['ref_total']} lokus kontrol"),
                         (k[4], "Rasio Cc/Ce", f"{ratio:.2f}" if ratio == ratio else "—", f"ambang {THRESH:.2f}")]:
        col.markdown(f"<div class='kpi accent'><div class='l'>{l}</div><div class='v'>{v}</div><div class='s'>{s}</div></div>",
                     unsafe_allow_html=True)
    if sim is not None:
        g1, g2 = st.columns([1, 1.5])
        with g1:
            st.plotly_chart(ratio_gauge(ratio), use_container_width=True, config={"displayModeBar": False})
        with g2:
            st.plotly_chart(sim_chart(ratio), use_container_width=True, config={"displayModeBar": False})
    st.markdown("<div class='note'>Rasio, ambang, dan estimasi % Arabika berasal dari simulasi copy-number yang <b>belum dikalibrasi</b> "
                "dengan sampel fisik — proporsi dihitung per genom (bukan % berat), bersifat indikatif, bukan klaim batas deteksi "
                f"atau kuantifikasi. Di atas {ARA_CAP}% selisih rasio terlalu kecil untuk dibedakan.</div>",
                unsafe_allow_html=True)

    # ---------- Modul 2 ----------
    st.markdown("<div class='sec'><span class='n'>02</span><span class='t'>Risiko kapang OTA</span>"
                "<span class='s'>Deteksi gen biosintesis okratoksin A</span></div>", unsafe_allow_html=True)
    if r["ota"]:
        oc = st.columns(len(r["ota"]))
        for col, (mk, n) in zip(oc, r["ota"].items()):
            pos = n > 0
            col.markdown(f"<div class='ota-chip'><div><div class='name'>{esc(mk)}</div><div class='sub'>{n} amplikon</div></div>"
                         f"<span class='pill {'pos' if pos else 'neg'}'>{'POSITIF' if pos else 'NEGATIF'}</span></div>",
                         unsafe_allow_html=True)
    st.markdown("<div class='note'>Skrining <b>risiko</b> berbasis DNA kapang — bukan pengukuran kadar toksin. "
                "Hasil positif harus dikonfirmasi dengan kuantifikasi OTA (HPLC).</div>", unsafe_allow_html=True)

    # ---------- Rekomendasi & rincian ----------
    st.markdown("<div class='sec'><span class='n'>03</span><span class='t'>Rekomendasi &amp; rincian marker</span></div>",
                unsafe_allow_html=True)
    d1, d2 = st.columns([1, 1.3])
    with d1:
        st.markdown("<div class='rec'><b>Tindak lanjut yang disarankan</b><ul>" +
                    "".join(f"<li>{esc(x)}</li>" for x in r["recs"]) + "</ul></div>", unsafe_allow_html=True)
        st.write("")
        st.download_button("Unduh Laporan QA (HTML · siap cetak PDF)", build_report(r),
                           file_name=f"KopiFied_QA_{r['label'].rsplit('.', 1)[0]}_{r['report_id']}.html",
                           mime="text/html", use_container_width=True)
        st.download_button("Unduh data per marker (CSV)", r["det"].to_csv(index=False),
                           file_name=f"KopiFied_markers_{r['report_id']}.csv", mime="text/csv",
                           use_container_width=True)
    with d2:
        st.plotly_chart(marker_chart(r["det"]), use_container_width=True, config={"displayModeBar": False})
    with st.expander("Tabel rincian per marker"):
        st.dataframe(r["det"], hide_index=True, use_container_width=True,
                     column_config={"Amplikon": st.column_config.ProgressColumn(
                         "Amplikon", format="%d", min_value=0, max_value=max(1, int(r["det"]["Amplikon"].max())))})

# ================= TABS BUKTI =================
st.markdown("<div class='sec'><span class='n'>◆</span><span class='t'>Bukti ilmiah &amp; dokumentasi panel</span></div>",
            unsafe_allow_html=True)
tabs = st.tabs(["Matriks validasi", "Primer", "Genom referensi", "Bukti SRA & multiplex", "Metode & batasan"])
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
    st.markdown("**Validasi sampel independen** — SRA SRR40981778 (Arabika, lini mutan azide)")
    sp = f"{BASE}/output/sra_check/sra_primer_sites.tsv"
    if os.path.exists(sp):
        st.dataframe(pd.read_csv(sp, sep="\t"), hide_index=True, use_container_width=True)
        st.caption("9/14 marker dengan ≥1 read cocok primer penuh pada ~1,3× coverage; marker OTA nol di reads kopi (sesuai harapan).")
    st.markdown("**Kompatibilitas multiplex** — `primer3.calc_heterodimer`")
    mp = f"{RES}/multiplex_check.tsv"
    if os.path.exists(mp):
        mux = pd.read_csv(mp, sep="\t")
        st.write(f"Hanya {len(mux)} dari 1.035 pasangan oligo (0,8%) berisiko dimer — semuanya interaksi probe–probe.")
        st.dataframe(mux, hide_index=True, use_container_width=True)
with tabs[4]:
    rd = f"{BASE}/README.md"
    if os.path.exists(rd):
        st.markdown(open(rd).read())
