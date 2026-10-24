#!/usr/bin/env python3
"""Buat laporan HTML mandiri KopiFied PoC v0.1 (gambar base64-embedded)."""
import base64, csv, html, os

BASE = "/work/ivanasitepu/genohack_Kopicert/build"
RES = f"{BASE}/results"

def img64(path):
    return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()

def tsv_table(path, max_rows=40, cls=""):
    rows = list(csv.reader(open(path), delimiter="\t"))
    h = "".join(f"<th>{html.escape(c)}</th>" for c in rows[0])
    body = ""
    for r in rows[1:max_rows+1]:
        body += "<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in r) + "</tr>"
    more = f"<p class='note'>… {len(rows)-1-max_rows} baris lainnya tidak ditampilkan</p>" if len(rows)-1 > max_rows else ""
    return f"<table class='{cls}'><thead><tr>{h}</tr></thead><tbody>{body}</tbody></table>{more}"

vm_png = img64(f"{RES}/validation_matrix.png")
bl_png = img64(f"{RES}/blend_simulation.png")
primers_tbl = tsv_table(f"{RES}/primers_v0.1.tsv")
funnel_tbl = tsv_table(f"{RES}/funnel.tsv")
ota_tbl = tsv_table(f"{RES}/ota_gene_presence.tsv")
genomes_tbl = tsv_table(f"{RES}/genome_qc.tsv", max_rows=30)
blend_tbl = tsv_table(f"{RES}/blend_simulation.tsv")

demo_rows = [
    ("arabica_caturra.fa", "ARABICA AUTHENTIC (rasio 1,00)", "Bersih", "ok"),
    ("robusta.fa", "FLAG: bukan Arabika", "Bersih", "warn"),
    ("blend_10pct.fa", "FLAG: campuran Robusta (rasio 1,11)", "Bersih", "warn"),
    ("blend_25pct.fa", "FLAG: campuran Robusta (rasio 1,25)", "Bersih", "warn"),
    ("blend_50pct.fa", "FLAG: campuran Robusta (rasio 2,00)", "Bersih", "warn"),
    ("arabica_with_carbonarius.fa", "ARABICA AUTHENTIC", "FLAG: kapang OTA terdeteksi", "mix"),
    ("mold_free.fa", "ARABICA AUTHENTIC", "Bersih", "ok"),
    ("mold_ota_positive.fa", "INVALID (bukan kopi)", "FLAG: kapang OTA terdeteksi", "mix"),
]
demo_tbl = "<table><thead><tr><th>Sampel demo</th><th>Modul 1 — Autentikasi</th><th>Modul 2 — Risiko OTA</th></tr></thead><tbody>"
for s, m1, m2, k in demo_rows:
    demo_tbl += f"<tr><td>{s}</td><td>{html.escape(m1)}</td><td>{html.escape(m2)}</td></tr>"
demo_tbl += "</tbody></table>"

page = f"""<!DOCTYPE html>
<html lang="id"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>KopiFied — PoC Panel Marker DNA v0.1</title>
<style>
 body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; background: #faf7f2; color: #2b2118; }}
 .wrap {{ max-width: 980px; margin: 0 auto; padding: 24px; }}
 header {{ background: linear-gradient(135deg, #4b2e1e, #6f4e37); color: #fff; padding: 36px 24px; border-radius: 12px; }}
 header h1 {{ margin: 0 0 6px; font-size: 28px; }}
 header p {{ margin: 4px 0; opacity: .92; }}
 .badge {{ display:inline-block; background:#c8a165; color:#2b2118; border-radius: 999px; padding: 3px 12px; font-size: 12px; font-weight: 600; margin-right: 6px; }}
 section {{ background: #fff; border-radius: 12px; padding: 20px 24px; margin-top: 18px; box-shadow: 0 1px 4px rgba(0,0,0,.06); }}
 h2 {{ color: #4b2e1e; border-bottom: 2px solid #e8dcc8; padding-bottom: 6px; font-size: 20px; }}
 table {{ border-collapse: collapse; width: 100%; font-size: 12.5px; }}
 th {{ background: #6f4e37; color: #fff; text-align: left; padding: 6px 8px; }}
 td {{ border-bottom: 1px solid #eee; padding: 5px 8px; }}
 tr:nth-child(even) td {{ background: #faf5ee; }}
 img {{ max-width: 100%; border-radius: 8px; border: 1px solid #e5dccb; }}
 .note {{ color: #7a6a58; font-size: 12px; }}
 .claim {{ background: #eef6ee; border-left: 4px solid #2ca02c; padding: 10px 14px; border-radius: 6px; }}
 .limit {{ background: #fdf0ee; border-left: 4px solid #d62728; padding: 10px 14px; border-radius: 6px; }}
 code {{ background: #f3ece2; padding: 1px 5px; border-radius: 4px; font-size: 12px; }}
 ol li, ul li {{ margin: 4px 0; }}
 footer {{ text-align: center; color: #9a8a76; font-size: 12px; padding: 24px; }}
</style></head><body><div class="wrap">

<header>
  <h1>☕ KopiFied — Panel Marker DNA v0.1</h1>
  <p>Autentikasi kopi Arabika &amp; deteksi risiko kapang penghasil ochratoxin A — proof-of-concept in silico</p>
  <p><span class="badge">GenoHack 2026</span><span class="badge">Tim Kimpulita</span><span class="badge">TRL 3</span><span class="badge">25 genom publik</span></p>
</header>

<section><h2>1. Apa yang kami klaim (TRL 3)</h2>
<div class="claim"><ol>
<li>Marker <b>M1-Ce</b> mengamplifikasi semua assembly <i>C. arabica</i> (3) dan <i>C. eugenioides</i> (2) publik, dan <b>tidak</b> mengamplifikasi <i>C. canephora</i> — bukti keaslian Arabika.</li>
<li>Marker <b>M1-Cc</b> mengamplifikasi arabika + canephora tetapi tidak eugenioides — pasangan kerangka rasio untuk deteksi campuran Robusta.</li>
<li>Marker <b>REF</b> (3 lokus redundan) mengamplifikasi ketiga spesies kopi — kontrol internal DNA kopi.</li>
<li>Marker <b>OTA</b> (gen otaB/otaD) mengamplifikasi kapang penghasil OTA (westerdijkiae, ochraceus, carbonarius) dan tidak mengamplifikasi kapang non-penghasil, kopi, atau mikroflora fermentasi.</li>
<li>16 pasangan primer lolos desain Primer3 + in-silico PCR (≤2 mismatch, 0 mismatch 3 basa ujung 3′), termasuk 3 marker REF redundan.</li>
</ol></div>
<div class="limit" style="margin-top:10px"><b>Yang TIDAK kami klaim:</b> persentase Robusta dalam campuran, kadar OTA, sensitivitas pada biji sangrai, dan batas deteksi. Semua menunggu validasi lab (TRL 4). Kurva simulasi berlabel “belum dikalibrasi”.</div>
</section>

<section><h2>2. Corong penemuan marker</h2>{funnel_tbl}
<p class="note">25 assembly genom publik (6 kopi + 16 kapang + 2 mikroflora + Penicillium) → kandidat dari alignment genom-ke-genom (minimap2) → verifikasi BLASTN → desain primer → in-silico PCR.</p></section>

<section><h2>3. Matriks validasi in-silico</h2>
<p>110/125 sel sesuai harapan (88%). Hijau = sesuai (OK+/OK−); merah muda/merah = ketidaksesuaian yang terdokumentasi (lihat bagian 7).</p>
<img src="{vm_png}" alt="Matriks validasi"></section>

<section><h2>4. Kerangka rasio marker (simulasi campuran)</h2>
<p>Model copy-number: arabika alotetraploid (Ce 2 kopi, Cc 2 kopi, REF 4 kopi per sel), robusta diploid (Cc 2 kopi, REF 2 kopi). Rasio Cc/Ce naik monoton terhadap persentase Robusta — inilah dasar estimasi campuran setelah kalibrasi lab.</p>
<img src="{bl_png}" alt="Kurva simulasi">{blend_tbl}</section>

<section><h2>5. Panel primer v0.1</h2>{primers_tbl}
<p class="note">Amplicon 70–130 bp (siap untuk DNA terfragmentasi); probe TaqMan bila tersedia. Detail per genom: <code>results/insilico_pcr.tsv</code>.</p></section>

<section><h2>5b. Kompatibilitas multiplex (heterodimer check)</h2>
<p>Cek silang semua 34 oligo (primer + probe) dengan <code>primer3.calc_heterodimer</code>: hanya <b>8 dari 1.035 pasangan (0,8%)</b> yang berisiko dimer (ΔG &lt; −9 kkal/mol), semuanya interaksi probe–probe yang dapat dipisah ke reaksi multiplex berbeda.</p>
<table><thead><tr><th>Oligo 1</th><th>Oligo 2</th><th>ΔG (kkal/mol)</th></tr></thead><tbody>
<tr><td>A1_0397:P</td><td>OTA-D-carb:P</td><td>−12,8</td></tr>
<tr><td>A1_0397:P</td><td>REF_0008:P</td><td>−11,4</td></tr>
<tr><td>REF_0008:P</td><td>REF_0013:P</td><td>−11,1</td></tr>
<tr><td>OTA-D-carb:P</td><td>REF_0015:P</td><td>−10,6</td></tr>
<tr><td>A1_0397:P</td><td>A1_0170:P</td><td>−10,3</td></tr>
</tbody></table>
<p class="note">Lengkap: <code>results/multiplex_check.tsv</code></p></section>
<section><h2>6. Kehadiran gen OTA per genom kapang (BLASTN, task blastn)</h2>{ota_tbl}
<p class="note">Produsen OTA kopi (westerdijkiae, ochraceus, carbonarius, steynii) membawa klaster; non-penghasil (niger, flavus, oryzae, tubingensis) bersih. <i>P. verrucosum</i> membawa homolog klaster (benar secara biologis — OTA+ di serealia) tetapi tidak diamplifikasi panel v0.1.</p></section>

<section><h2>7. Uji demo end-to-end (8 sampel FASTA dari genom asli)</h2>{demo_tbl}
<p class="note">Sampel = fragmen 150 bp dari genom publik + amplikon marker; keputusan dihitung oleh logika yang sama dengan dashboard. Waktu analisis &lt; 1 detik per sampel.</p></section>

<section><h2>8. Keputusan desain &amp; keterbatasan (dokumentasi jujur)</h2>
<ol>
<li>Pendekatan alignment (minimap2) dipakai sebagai rencana cadangan KMC karena batas RAM kontainer; logika biologis identik.</li>
<li>Primer OTA belum menangkap <i>A. steynii</i> dan sebagian <i>A. ochraceus</i> — klaster mereka divergen dari westerdijkiae; produksi butuh primer degenerat / set per klad.</li>
<li>Panel diperkuat 3 marker REF redundan setelah REF_0008 terbukti lemah di satu assembly arabika.</li>
<li>Marker OTA menarget <b>gen biosintesis</b>, bukan spesies — lebih dekat ke risiko toksin, tetapi bukan pengukuran toksin.</li>
<li>Risiko diketahui: introgresi Timor Hybrid pada kultivar Catimor/Sarchimor dapat menaikkan sinyal M1-Cc; mitigasi: kalibrasi per asal geografis.</li>
</ol></section>

<section><h2>9. Kurasi data &amp; reproduksibilitas</h2>
<p>25 genom publik NCBI terverifikasi integritas (checksum) — ringkasan QC:</p>{genomes_tbl}
<p class="note">Pipeline dapat dijalankan ulang dengan <code>bash run_all.sh</code> (environment.yml terlampir). Tools: minimap2 2.31, BLAST 2.17, seqkit 2.13, primer3-py 2.4.0, Python 3.12.</p></section>

<section><h2>10. Validasi pada sampel independen (SRA)</h2>
<p>Situs primer dicek pada reads resekuensing arabika dari studi lain (SRA <code>SRR40981778</code>, lini mutan azide — materi biologis berbeda dari ketiga assembly referensi). Pada cakupan ~1,3× (8 juta reads): 9/14 marker memiliki ≥1 read yang cocok dengan primer penuh (≤1 mismatch); marker OTA tidak ditemukan sama sekali pada reads kopi (sesuai harapan). Ini bukti pendukung konservasi situs primer pada sampel independen; cakupan penuh memerlukan reads lebih dalam.</p>
</section>

<footer>KopiFied PoC v0.1 — Tim Kimpulita — GenoHack 2026 · ITB Hackathon Series · dibuat 24 Okt 2026<br>
In-silico proof of concept; bukan pengganti uji laboratorium.</footer>
</div></body></html>"""

out = "/work/ivanasitepu/genohack_Kopicert/KopiFied_PoC_report_v0.1.html"
open(out, "w").write(page)
print(out, os.path.getsize(out) // 1024, "KB")
