#!/usr/bin/env bash
# =============================================================================
# RESIST | Module D input preparation
# prepare_cell_annotations.sh
# -----------------------------------------------------------------------------
# Purpose  : Join a Cell Ranger graph-based clustering table with its t-SNE
#            projection and prefix each barcode with the run accession, producing the
#            per-sample cell annotation table used when assembling multi-run objects.
# Inputs   : <cellranger_outs>/analysis/clustering/gene_expression_graphclust/clusters.csv
#            <cellranger_outs>/analysis/tsne/gene_expression_2_components/projection.csv
# Outputs  : <out_dir>/<accession>_annots.csv - cell_id, cluster, tsne_1, tsne_2
# Usage    : bash scripts/D/prepare_cell_annotations.sh <cellranger_outs> <accession> <out_dir>
# Origin   : APA_analysis/Untitled-1.sh (four repeated join commands)
# Notes    : The original wrote every sample to the same `annots.csv`, so only the
#            last one survived; the output name now carries the accession.
# =============================================================================

set -euo pipefail

[[ $# -eq 3 ]] || { echo "usage: $0 <cellranger_outs> <accession> <out_dir>" >&2; exit 2; }
OUTS="$1"; ACC="$2"; DEST="$3"
mkdir -p "$DEST"

CLUST="${OUTS}/analysis/clustering/gene_expression_graphclust/clusters.csv"
TSNE="${OUTS}/analysis/tsne/gene_expression_2_components/projection.csv"
for f in "$CLUST" "$TSNE"; do
  [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
done

join -t',' "$CLUST" "$TSNE" \
  | sed -E "s/([ACGT]{16})-1/${ACC}_count_\1/" \
  | sed "1 s/^.*$/cell_id,cluster,tsne_1,tsne_2/" \
  > "${DEST}/${ACC}_annots.csv"

echo "[utils] wrote ${DEST}/${ACC}_annots.csv"
