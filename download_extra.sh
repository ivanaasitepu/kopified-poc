#!/bin/bash
# Unduh non-target tambahan sesuai doc Drey (Tahap 1 pelengkap)
BASE=/work/ivanasitepu/genohack_Kopicert/build/data/genomes
cd $BASE
dl() {
  local acc="$1" path="$2" label="$3"
  local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${path}/${acc}_${label}/${acc}_${label}_genomic.fna.gz"
  local f="${acc}_${label}_genomic.fna.gz"
  [ -f "${f}.DONE" ] && { echo "SKIP $acc"; return 0; }
  for i in $(seq 1 60); do
    curl -s --proxy-insecure -k -L --connect-timeout 60 --speed-time 120 --speed-limit 500 -C - "$url" -o "$f" || true
    if gzip -t "$f" 2>/dev/null; then touch "${f}.DONE"; echo "OK $acc $(du -sh $f | cut -f1)"; return 0; fi
    sleep 5
  done
  echo "FAIL $acc"
}
dl GCA_986281605.1 "GCA/986/281/605" "DAH1045FM"          # A. flavus NRRL 3357
dl GCA_000184455.4 "GCA/000/184/455" "ASM184455v4"        # A. oryzae RIB40
dl GCA_988135455.1 "GCA/988/135/455" "TATVAM_FUN7"        # A. tubingensis
dl GCF_028828655.1 "GCF/028/828/655" "ASM2882865v1"       # P. verrucosum
dl GCA_051861825.1 "GCA/051/861/825" "ASM5186182v1"       # S. cerevisiae
dl GCF_061378395.1 "GCF/061/378/395" "ASM6137839v1"       # L. plantarum
echo "DONE EXTRA"
