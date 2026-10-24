#!/bin/bash
# run_all.sh — KopiFied PoC: jalankan ulang seluruh pipeline (Tahap 1-5)
# Penggunaan: bash run_all.sh
set -e
cd "$(dirname "$0")"
export PATH=/public/tools/blast/2.17.0:/public/tools/seqkit/2.13.0:$PATH
mkdir -p logs
echo "=== KopiFied PoC pipeline — $(date) ==="

echo "[1/6] Unduh genom (skip jika sudah ada .DONE)"
bash download_genomes.sh  2>&1 | tee logs/01_download.log
bash download_extra.sh    2>&1 | tee -a logs/01_download.log

echo "[2/6] Assay 1: marker M1-Ce (autentikasi Arabica)"
python3 pipeline/assay1_markers.py 2>&1 | tee logs/02_assay1.log

echo "[3/6] Assay 1b: marker M1-Cc + REF (kerangka rasio)"
python3 pipeline/assay1b_cc_ref.py 2>&1 | tee logs/03_assay1b.log

echo "[4/6] Assay 2: panel marker kapang OTA"
python3 pipeline/assay2_markers.py 2>&1 | tee logs/04_assay2.log
python3 pipeline/assay2d_genes_from_cluster.py 2>&1 | tee logs/05_assay2d.log

echo "[5/6] Primer + in-silico PCR (Tahap 3)"
python3 pipeline/assay3_primers.py 2>&1 | tee logs/06_primers.log

echo "[6/6] Matriks validasi + simulasi campuran (Tahap 5)"
python3 pipeline/assay5_validate.py 2>&1 | tee logs/07_validate.log

echo "=== selesai ==="
minimap2 --version; blastn -version | head -1; seqkit version; python3 -c "import primer3; print('primer3-py', primer3.__version__)"
