#!/usr/bin/env bash
# =============================================================================
# RESIST | Module D - Immunogenomic features of resistance
# D4_run_optitype.sh
# -----------------------------------------------------------------------------
# Purpose  : Call MHC class I genotypes (HLA-A, -B, -C) for one sample with OptiType.
#            The alleles it reports define which peptides can be presented and are
#            therefore required input to the neoantigen prediction in D5.
# Inputs   : One FASTQ (single-end) or two FASTQs (paired-end) of reads that survived
#            HLA-region fishing
# Outputs  : <out_dir>/<timestamp>_result.tsv - predicted class I alleles
# Usage    : bash D4_run_optitype.sh --pipeline /path/to/OptiTypePipeline.py \
#            --type rna --out <out_dir> <fastq1> [fastq2]
# Origin   : README.md, PART II Step 1 (command documented, not previously scripted)
# Notes    : Requires an OptiType installation: https://github.com/FRED-2/OptiType
# =============================================================================

set -euo pipefail

PIPELINE=""
SEQTYPE="rna"
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pipeline) PIPELINE="$2"; shift 2 ;;
    --type)     SEQTYPE="$2";  shift 2 ;;
    --out)      OUTDIR="$2";   shift 2 ;;
    -h|--help)
      echo "usage: $0 --pipeline <OptiTypePipeline.py> [--type rna|dna] --out <dir> <fq1> [fq2]" >&2
      exit 0 ;;
    *) break ;;
  esac
done

if [[ -z "$PIPELINE" || -z "$OUTDIR" || $# -lt 1 ]]; then
  echo "usage: $0 --pipeline <OptiTypePipeline.py> [--type rna|dna] --out <dir> <fq1> [fq2]" >&2
  exit 2
fi
case "$SEQTYPE" in rna|dna) ;; *) echo "error: --type must be rna or dna" >&2; exit 2 ;; esac

mkdir -p "$OUTDIR"
python "$PIPELINE" -i "$@" "--${SEQTYPE}" --outdir "$OUTDIR"

echo "[D4] OptiType results in $OUTDIR"
