#!/bin/bash
# KopiFied Tahap 0 — download genomes with retry+resume loops
BASE=/work/ivanasitepu/genohack_Kopicert/build/data/genomes
mkdir -p $BASE
cd $BASE

dl() {
  local acc="$1" path="$2" label="$3"
  local url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${path}/${acc}_${label}/${acc}_${label}_genomic.fna.gz"
  local f="${acc}_${label}_genomic.fna.gz"
  [ -f "${f}.DONE" ] && { echo "SKIP $acc (done)"; return 0; }
  for i in $(seq 1 60); do
    curl -s --proxy-insecure -k -L --connect-timeout 60 --speed-time 120 --speed-limit 500 -C - "$url" -o "$f" || true
    if gzip -t "$f" 2>/dev/null; then touch "${f}.DONE"; echo "OK $acc after $i attempt(s), $(du -sh $f | cut -f1)"; return 0; fi
    sleep 5
  done
  echo "FAIL $acc"
  return 1
}

# fungi first (small, quick wins)
dl GCA_001990825.1 "GCA/001/990/825" "Aspca3"                 # A. carbonarius
dl GCF_002849105.1 "GCF/002/849/105" "Aspste1"                # A. steynii
dl GCA_001307345.1 "GCA/001/307/345" "ASM130734v1"            # A. westerdijkiae
dl GCA_026259865.1 "GCA/026/259/865" "ASM2625986v1"           # A. westerdijkiae 2
dl GCA_019976475.1 "GCA/019/976/475" "ASM1997647v1"           # A. westerdijkiae 3
dl GCA_004849945.1 "GCA/004/849/945" "AoFc"                   # A. ochraceus
dl GCA_025594125.1 "GCA/025/594/125" "ASM2559412v1"           # A. ochraceus 2
dl GCA_005784425.1 "GCA/005/784/425" "ASM578442v1"            # A. ochraceus 3
dl GCA_026108165.1 "GCA/026/108/165" "ASM2610816v1"           # A. ochraceus 4
dl GCA_037127345.1 "GCA/037/127/345" "ASM3712734v1"           # A. ochraceus 5
dl GCA_051529875.1 "GCA/051/529/875" "UCR_AnigDFA-N_1.0"      # A. niger (non-target, chr)
dl GCA_047651775.1 "GCA/047/651/775" "ASM4765177v1"           # A. niger 2 (non-target)

# coffee (larger)
dl GCF_003713205.1 "GCF/003/713/205" "Ceug_1.0"               # C. eugenioides (JHU ref)
dl GCA_036785875.1 "GCA/036/785/875" "ASM3678587v1"           # C. eugenioides 2
dl GCA_900059795.1 "GCA/900/059/795" "AUK_PRJEB4211_v1"       # C. canephora (Genoscope)
dl GCA_036785865.1 "GCA/036/785/865" "ASM3678586v1"           # C. canephora 2
dl GCF_003713225.1 "GCF/003/713/225" "Cara_1.0"               # C. arabica (JHU ref)
dl GCF_036785885.1 "GCF/036/785/885" "Coffea_Arabica_ET-39_HiFi"  # C. arabica 2 (HiFi)
dl GCA_030873655.1 "GCA/030/873/655" "IGA_Cara_2.4"           # C. arabica 3

echo "=== SUMMARY ==="
ls -la $BASE/*.fna.gz 2>/dev/null
echo "DONE ALL"
