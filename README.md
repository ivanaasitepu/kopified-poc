# KopiFied PoC — Panel Marker DNA untuk Autentikasi Kopi Arabika & Deteksi Risiko Kapang OTA

Proof-of-concept in silico untuk GenoHack 2026 (Tim Kimpulita).

## Klaim (TRL 3)

Panel marker DNA yang:

1. mengenai setiap assembly *C. arabica* dan *C. eugenioides* publik dan tidak mengenai *C. canephora* (marker M1-Ce);
2. mengenai komponen subgenom canephora pada arabika (M1-Cc) sebagai pasangan kerangka rasio;
3. mengenai genom kapang penghasil OTA (gen otaA/otaB/otaC/otaR/otaD) dan tidak mengenai kapang non-penghasil, kopi, atau mikroflora umum;
4. primer dan probe yang lolos desain Primer3, in-silico PCR, dan BLAST off-target;
5. kerangka rasio marker yang *bisa* dikalibrasi untuk persentase Robusta di kemudian hari.

## Yang TIDAK kami klaim

Persentase Robusta dalam campuran, kadar OTA, sensitivitas pada biji sangrai, dan batas deteksi. Semua itu menunggu validasi lab (TRL 4). Kurva simulasi campuran berlabel "belum dikalibrasi".

## Cara menjalankan (3 perintah)

```bash
conda env create -f environment.yml && conda activate kopified-poc
bash run_all.sh
streamlit run dashboard/app.py
```

## Struktur

- `config/genomes.tsv` — 25 genom publik terkurasi (aksesi, peran, checksum md5)
- `results/genome_qc.tsv` — QC per assembly
- `output/assay1/` — alignment + kandidat + panel M1-Ce
- `results/ota_gene_presence.tsv` — kehadiran gen OTA per genom kapang
- `results/` — marker panel, primer, matriks validasi, simulasi

## Keputusan desain (dokumentasi jujur)

1. **Pendekatan alignment, bukan k-mer.** Rencana utama KMC k=31; kami memakai rencana cadangan alignment (minimap2 asm5) karena keterbatasan RAM kontainer. Logika biologis identik.
2. **_P. verrucosum_ membawa homolog klaster OTA** (benar secara biologis — OTA+ di serealia), tetapi primer panel v0.1 tidak mengamplifikasinya (terlalu divergen dari klaster kopi). Keputusan PoC: panel fokus pada produsen OTA di kopi (seksi Circumdati + carbonarius); deteksi P. verrucosum masuk roadmap.
3. **_A. ochraceus_ terdeteksi parsial** (identitas tinggi, cakupan rendah) — klaster OTA-nya divergen dari westerdijkiae; primer produksi perlu query multi-spesies atau primer degenerat. Dicatat sebagai pekerjaan lanjutan.
4. **Marker gen, bukan spesies, untuk OTA** — dalam satu spesies (A. niger) ada strain OTA+ dan OTA−; marker gen lebih dekat ke risiko toksin, tetapi tetap bukan pengukuran toksin.
5. **M1-cp (kloroplas)** tidak dikerjakan di PoC ini — kandidat pengembangan berikutnya (sensitif untuk DNA terdegradasi, tetapi tidak bisa dipakai untuk rasio).

## Apa yang belum kami ketahui

- Batas deteksi pada sampel nyata dan pada biji sangrai (DNA terfragmentasi) — perlu validasi lab.
- Efek introgresi Timor Hybrid pada kultivar Catimor/Sarchimor: sinyal M1-Cc bisa lebih tinggi dari arabika murni tipikal — mitigasi: kalibrasi per asal geografis.
- Stabilitas rasio Ce/Cc antar kultivar arabika Indonesia (Gayo, Lintong, Kintamani) — menunggu data resekuensing.

## Versi

v0.1-hackathon — 25 Oktober 2026
