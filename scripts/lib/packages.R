# =============================================================================
# RESIST | scripts/lib/packages.R
# -----------------------------------------------------------------------------
# Packages required by every module entry point. Analysis-specific packages
# (CellChat, GSVA, clusterProfiler, fgsea, cmapR, ...) stay declared at the top
# of the script that needs them, so that each entry point remains readable as a
# self-contained method description.
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(Seurat)
  library(patchwork)
  library(stringr)
  library(ggsci)
  library(R.utils)
  library(scales)
  library(tibble)
})
