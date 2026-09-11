#!/usr/bin/env bash
# =============================================================================
# RESIST | scripts/utils/rename_results.sh
# -----------------------------------------------------------------------------
# Purpose : Strip a dataset-specific prefix from result file names, so that
#           files written as <GSE>_seurat_afterAnno.RDS_<plot>.pdf become
#           <plot>.pdf for upload to the RESIST browser.
# Usage   : bash scripts/utils/rename_results.sh <results_dir> <prefix>
#           e.g. bash scripts/utils/rename_results.sh data/results/B GSE104987
# Origin  : scripts/rename.sh
# Notes   : Existing targets are never overwritten; conflicts are reported.
# =============================================================================
set -euo pipefail

[[ $# -eq 2 ]] || { echo "usage: $0 <results_dir> <prefix>" >&2; exit 2; }
DIR="$1"; PREFIX="$2"
cd "$DIR"

shopt -s nullglob
for f in "${PREFIX}"_*; do
  new="${f#${PREFIX}_seurat_afterAnno.RDS_}"
  new="${new#${PREFIX}_deg.sig_}"
  new="${new#${PREFIX}_}"
  [[ "$new" != "$f" ]] || continue
  if [[ -e "$new" ]]; then
    echo "SKIP (target exists): $f -> $new"
    continue
  fi
  mv -- "$f" "$new"
done
