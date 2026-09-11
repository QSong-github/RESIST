#!/usr/bin/env bash
# =============================================================================
# RESIST | Module D - Immunogenomic features of resistance
# D2_prepare_cellranger_fastqs.sh
# -----------------------------------------------------------------------------
# Purpose  : Reshape SRA-style FASTQ triplets into the file names Cell Ranger
#            expects, either by symlinking one run per sample or by concatenating
#            several runs that belong to the same library.
# Inputs   : A directory of per-run FASTQs named <SRR>_{1,2,3}.fastq.gz
# Outputs  : cellranger_fastqs/<name>_S1_L001_{R1,R2,I1}_001.fastq.gz
# Usage    : bash D2_prepare_cellranger_fastqs.sh link   <src_dir> <out_dir> <SRR> [SRR ...]
#            bash D2_prepare_cellranger_fastqs.sh merge  <src_dir> <out_dir> <sample> <SRR> [SRR ...]
# Origin   : APA_analysis/Untitled-1.sh
# Notes    : Read 2 of the SRA triplet carries the cell barcode and UMI and becomes
#            R1; read 3 carries the cDNA and becomes R2. Verify this against your
#            own submission before running.
# =============================================================================

set -euo pipefail

usage() {
  cat >&2 <<'USAGE'
usage:
  D2_prepare_cellranger_fastqs.sh link  <src_dir> <out_dir> <SRR> [SRR ...]
      One run per sample. Symlinks <SRR>_2 -> R1 and <SRR>_3 -> R2.

  D2_prepare_cellranger_fastqs.sh merge <src_dir> <out_dir> <sample> <SRR> [SRR ...]
      Several runs of one library. Concatenates reads 1, 2 and 3 of every run
      into <sample>_S1_L001_{R1,R2,I1}_001.fastq.gz.
USAGE
  exit 2
}

[[ $# -ge 4 ]] || usage
MODE="$1"; SRC="$2"; OUT="$3"; shift 3
mkdir -p "$OUT"

case "$MODE" in
  link)
    for srr in "$@"; do
      ln -sf "${SRC}/${srr}/${srr}_2.fastq.gz" "${OUT}/${srr}_S1_L001_R1_001.fastq.gz"
      ln -sf "${SRC}/${srr}/${srr}_3.fastq.gz" "${OUT}/${srr}_S1_L001_R2_001.fastq.gz"
      echo "[D2] linked ${srr}"
    done
    ;;
  merge)
    [[ $# -ge 2 ]] || usage
    SAMPLE="$1"; shift
    for idx in 1 2 3; do
      case "$idx" in 1) tag=R1 ;; 2) tag=R2 ;; 3) tag=I1 ;; esac
      dest="${OUT}/${SAMPLE}_S1_L001_${tag}_001.fastq.gz"
      : > "$dest"
      for srr in "$@"; do
        cat "${SRC}/${srr}/${srr}_${idx}.fastq.gz" >> "$dest"
      done
      echo "[D2] wrote $dest"
    done
    ;;
  *) usage ;;
esac
