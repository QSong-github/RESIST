#!/usr/bin/env Rscript
# =============================================================================
# RESIST | config/setup/install_r_packages.R
# -----------------------------------------------------------------------------
# Installs the R packages that are not available from conda, or that must be
# taken from GitHub. Run once after creating the conda environment.
#   Rscript config/setup/install_r_packages.R
# =============================================================================

cran <- c("tidyverse", "Seurat", "yaml", "jsonlite", "RColorBrewer", "ggsci", "ggrepel",
          "cowplot", "patchwork", "pheatmap", "viridis", "reshape2", "scales",
          "broom", "msigdbr", "treemap", "treemapify", "packcircles", "R.utils", "devtools")
bioc <- c("SingleCellExperiment", "GenomicRanges", "rtracklayer",
          "clusterProfiler", "GSVA", "GSEABase", "fgsea", "biomaRt",
          "org.Hs.eg.db", "org.Mm.eg.db", "cmapR", "cytolib", "flowCore", "ComplexHeatmap")
gh   <- c(CellChat = "jinworks/CellChat")

need <- function(p) !requireNamespace(p, quietly = TRUE)

miss <- Filter(need, cran)
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

if (need("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
miss <- Filter(need, bioc)
if (length(miss)) BiocManager::install(miss, ask = FALSE, update = FALSE)

for (p in names(gh)) if (need(p)) devtools::install_github(gh[[p]])

missing <- Filter(need, unique(c(cran, bioc, names(gh))))
if (length(missing)) stop("Packages still missing: ", paste(missing, collapse = ", "))
cat("R dependencies satisfied.\n")
