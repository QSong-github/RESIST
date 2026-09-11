#!/usr/bin/env bash
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C1_call_variants_cellsnp.sh
# -----------------------------------------------------------------------------
# Purpose  : Genotype every cell at a panel of reference SNP positions with
#            cellSNP-lite, the first step of the resistance-associated variant
#            analysis. Run once per condition (sensitive, resistant) so that the
#            Fisher test in C3 can contrast the two.
# Inputs   : <bam>       position-sorted, indexed alignment for one condition
#            <barcodes>  one cell barcode per line
#            <ref_snps>  reference SNP VCF
# Outputs  : <outdir>/cellSNP.cells.vcf.gz - per-cell AD/DP, consumed by C2
# Usage    : bash C1_call_variants_cellsnp.sh <bam> <barcodes> <ref_snps> <outdir> [threads]
# Origin   : README.md, PART I Step 1 (command documented, not previously scripted)
# Notes    : Requires cellSNP-lite on PATH: https://github.com/single-cell-genetics/cellsnp-lite
# =============================================================================

set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "usage: $0 <bam> <barcodes.txt> <reference_snps.vcf> <outdir> [threads]" >&2
  exit 2
fi

BAM="$1"
BARCODES="$2"
REF_SNPS="$3"
OUTDIR="$4"
THREADS="${5:-8}"

command -v cellsnp-lite >/dev/null 2>&1 || {
  echo "error: cellsnp-lite not found on PATH" >&2; exit 127; }

mkdir -p "$OUTDIR"

# Thresholds follow the RESIST protocol: variants are retained at a minor
# allele frequency of at least 0.1 and a total depth of at least 20 reads,
# which keeps per-cell allele-frequency estimates interpretable while
# discarding positions supported by a handful of UMIs.
cellsnp-lite \
    -s "$BAM" \
    -b "$BARCODES" \
    -O "$OUTDIR" \
    -R "$REF_SNPS" \
    --cellTAG CB \
    --UMItag UB \
    --gzip \
    --genotype \
    --minMAF 0.1 \
    --minCOUNT 20 \
    -p "$THREADS"

echo "[C1] wrote ${OUTDIR}/cellSNP.cells.vcf.gz"
