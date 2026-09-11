#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module D - Immunogenomic features of resistance
# D3_apa_quantification_and_plots.R
# -----------------------------------------------------------------------------
# Purpose  : Combine 10x gene expression with scUTRquant poly(A)-site
#            quantification to measure 3' UTR relative expression (RE) per
#            cell, then compare RE between sensitive and resistant states
#            overall and within each cell type. Produces the violin and ECDF
#            panels of the APA section of Module D.
# Inputs   : config/apa_config.yaml   - GTF, annotation object, cohort id
#            config/apa_samples.csv   - sample_id, group, tenx_dir, txs_rds
# Outputs  : data/results/D/APA/<cohort_id>/
#              RE_outputs/          - RE matrices and differential tables
#              violin_by_celltype/  - per-cell-type RE violins
#              *.pdf / *.png        - cohort-level violin and ECDF figures
# Usage    : Rscript scripts/D/D3_apa_quantification_and_plots.R
# Origin   : APA_analysis/APA_processing_graph.R
# Notes    : `annot_info` is consumed by map_from_annot_to_sce() and
#            diagnose_barcode_overlap() but is never assigned in the inherited
#            code; the annotation object named in the configuration has to be
#            loaded into it before this script will run end to end. See
#            docs/repository-guide.md, Known issues.
# =============================================================================

## --- RESIST bootstrap --------------------------------------------------------
## Locates the repository root (RESIST_HOME, the invoking script, or an upward
## search for config/config.yaml) and loads the shared configuration, palettes
## and I/O helpers. Identical in every RESIST entry point.
.resist_home <- local({
  h <- Sys.getenv("RESIST_HOME", unset = NA_character_)
  if (!is.na(h) && nzchar(h)) return(normalizePath(h, mustWork = TRUE))
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else normalizePath(getwd())
  while (!file.exists(file.path(d, "config", "config.yaml")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "config", "config.yaml")))
    stop("cannot locate the RESIST repository root; set RESIST_HOME.", call. = FALSE)
  d
})
source(file.path(.resist_home, "scripts", "lib", "init.R"))
## -----------------------------------------------------------------------------

# Load all required packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(GenomicRanges)
  library(rtracklayer)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(readr)
  library(tidyr)
  library(rlang)
})

## -----------------------------------------------------------------------------
## 1. CONFIGURATION & PATHS
## -----------------------------------------------------------------------------
# Cohort configuration lives in `config/apa_config.yaml`; the sample registry
# lives in the CSV that file names. Nothing below this block needs editing to
# process a new cohort.
APA_CFG <- yaml::read_yaml(Sys.getenv(
  "RESIST_APA_CONFIG",
  unset = file.path(RESIST_HOME, "config", "apa_config.yaml")))

.apa_path <- function(p) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) p else file.path(RESIST_HOME, p)
}

# gff_path: custom UTROME GTF used for PAS coordinate mapping
gff_path <- .apa_path(APA_CFG$gff_path)

# annotation_seurat_path: reference Seurat object carrying cell-type metadata
annotation_seurat_path <- .apa_path(APA_CFG$annotation_seurat_path)

# Output directory structure, under the Module D results tree
OUTDIR <- file.path(RESULTS_D, "APA", APA_CFG$cohort_id)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "violin_by_celltype"), showWarnings = FALSE)
dir.create(file.path(OUTDIR, "RE_outputs"), showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 2. SAMPLE DEFINITION
## -----------------------------------------------------------------------------
# One row per sample, mapping the 10x count matrix to its scUTRquant TXS
# object. Required columns: sample_id, group, tenx_dir, txs_rds.
sample_sheet <- .apa_path(APA_CFG$sample_sheet)
sample_tbl   <- readr::read_csv(sample_sheet, show_col_types = FALSE)
stopifnot(all(c("sample_id", "group", "tenx_dir", "txs_rds") %in% names(sample_tbl)))

## -----------------------------------------------------------------------------
## 3. CORE COMPUTATIONAL FUNCTIONS
## -----------------------------------------------------------------------------

#' Aggregate Sparse Matrix by Groups
#' @param mat Sparse matrix (dgCMatrix)
#' @param group Vector of group identifiers
aggregate_sparse_rowsum <- function(mat, group) {
  stopifnot(nrow(mat) == length(group))
  ug <- unique(group)
  G  <- sparseMatrix(i = match(group, ug), j = seq_len(nrow(mat)),
                     x = 1, dims = c(length(ug), nrow(mat)),
                     dimnames = list(ug, NULL))
  G %*% mat
}

#' Normalize 10x Barcodes
norm_cells_10x <- function(ids, sample_prefix) {
  raw_bc <- sub("^.*_", "", ids)
  bc     <- sub("-\\d+$", "", raw_bc)
  paste0(sample_prefix, "_", bc)
}

#' Normalize TXS Barcodes
norm_cells_txs <- function(ids, sample_prefix) {
  bc <- sub("^.*_", "", ids)
  bc <- sub("-\\d+$", "", bc)
  paste0(sample_prefix, "_", bc)
}

#' Calculate Relative Expression (RE) of PAS usage
#' @param sce SingleCellExperiment object from scUTRquant
#' @param gff_path GTF file for PAS distance calculation
#' @return List containing RE matrix, mean RE per cell, and gene-level counts
compute_RE_from_txs_sce <- function(sce, gff_path) {
  cnt       <- assay(sce, "counts")
  gene_name <- as.vector(rowData(sce)$gene_name)
  
  rr <- rowRanges(sce)
  if (inherits(rr, "GRangesList")) {
    pas_str <- vapply(rr, function(gr) as.character(strand(gr)[1]), character(1))
    pas_pos <- ifelse(
      pas_str == "+",
      vapply(rr, function(gr) max(end(gr)),   integer(1)),
      vapply(rr, function(gr) min(start(gr)), integer(1))
    )
  } else if (inherits(rr, "GRanges")) {
    pas_str <- as.character(strand(rr))
    pas_pos <- ifelse(pas_str == "+", end(rr), start(rr))
  } else stop("Invalid rowRanges format.")
  
  keep <- pas_str %in% c("+", "-")
  cnt <- cnt[keep, ]; gene_name <- gene_name[keep]; pas_str <- pas_str[keep]; pas_pos <- pas_pos[keep]
  
  prox_flag <- logical(nrow(cnt))
  idx_plus  <- which(pas_str == "+")
  idx_minus <- which(pas_str == "-")
  if (length(idx_plus)) {
    min_by_gn_plus <- tapply(pas_pos[idx_plus],  gene_name[idx_plus],  min)
    prox_flag[idx_plus]  <- pas_pos[idx_plus]  == min_by_gn_plus[gene_name[idx_plus]]
  }
  if (length(idx_minus)) {
    max_by_gn_minus <- tapply(pas_pos[idx_minus], gene_name[idx_minus], max)
    prox_flag[idx_minus] <- pas_pos[idx_minus] == max_by_gn_minus[gene_name[idx_minus]]
  }
  prox_rows   <- which(prox_flag)
  prox_gnames <- gene_name[prox_rows]
  prox_pos    <- pas_pos[prox_rows]
  
  utrome <- rtracklayer::import(gff_path)
  if ("gene_name" %in% colnames(mcols(utrome))) {
    gname_u <- as.character(mcols(utrome)$gene_name)
  } else {
    attr_txt <- as.character(mcols(utrome)$attribute)
    gname_u  <- sub('.*gene_name "([^"]+)".*', '\\1', attr_txt)
  }
  str_u <- as.character(strand(utrome))
  stop_side <- ifelse(str_u == "+", start(utrome), end(utrome))
  last_side <- ifelse(str_u == "+", end(utrome), start(utrome))
  names(stop_side) <- names(last_side) <- gname_u
  
  ok <- prox_gnames %in% names(stop_side) & prox_gnames %in% names(last_side)
  prox_rows   <- prox_rows[ok]
  prox_gnames <- prox_gnames[ok]
  prox_pos    <- prox_pos[ok]
  
  cUTR_len <- abs(prox_pos - stop_side[prox_gnames]); cUTR_len[cUTR_len < 1] <- 1
  aUTR_len <- abs(last_side[prox_gnames] - prox_pos);  aUTR_len[aUTR_len < 1] <- 1
  
  prox_cnt  <- cnt[prox_rows, , drop = FALSE]
  gene_sum  <- aggregate_sparse_rowsum(cnt, group = gene_name)
  gene_sum  <- gene_sum[prox_gnames, , drop = FALSE]
  other_cnt <- gene_sum - prox_cnt
  other_cnt@x[other_cnt@x < 0] <- 0
  
  RD_cUTR <- prox_cnt / cUTR_len
  RD_aUTR <- other_cnt / aUTR_len
  eps <- 1e-8
  RE <- log2((RD_aUTR + eps) / (RD_cUTR + eps))
  rownames(RE) <- prox_gnames
  
  RE_mean <- as.numeric(Matrix::colMeans(RE, na.rm = TRUE))
  names(RE_mean) <- colnames(cnt)
  
  list(RE = RE, RE_mean = RE_mean,
       gene_counts_all_genes = aggregate_sparse_rowsum(cnt, group = gene_name))
}

## -----------------------------------------------------------------------------
## 4. METADATA MAPPING & UTILITIES
## -----------------------------------------------------------------------------

ensure_writable_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(normalizePath(path))
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 1e-4) "****" else if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "ns"
}

load_annotation_seurat <- function(path, object_name = NULL) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  seus <- nm[sapply(nm, function(x) inherits(e[[x]], "Seurat"))]
  if (!is.null(object_name)) {
    if (!object_name %in% seus) stop("Named annotation object is not a Seurat object: ", object_name)
    return(e[[object_name]])
  }
  if (length(seus) != 1L) stop("Annotation RData must contain exactly one Seurat object; save the intended object as RDS.")
  return(e[[seus[1]]])
}

get_annot_meta <- function(seu_annot,
                           ct_col_candidates = c("sctype_classification","sctype_class","cell_type","sctype","celltype","CellType","annotation","annot","predicted.id"),
                           cond_col_candidates = c("condition","Condition","group","Group","status","Status")) {
  df <- seu_annot@meta.data |> tibble::rownames_to_column("barcode_full")
  suf_dash_re <- "-\\d+(?:_\\d+)*$"
  df$barcode_pure_dash <- ifelse(grepl(suf_dash_re, df$barcode_full), sub(suf_dash_re, "", df$barcode_full), df$barcode_full)
  df$suffix_dash       <- ifelse(grepl(suf_dash_re, df$barcode_full), sub(paste0(".*(", suf_dash_re, ")$"), "\\1", df$barcode_full), "")
  df$barcode_pure_dot <- ifelse(grepl("\\.", df$barcode_full), sub("\\..*$", "", df$barcode_full), df$barcode_full)
  df$suffix_dot       <- ifelse(grepl("\\.", df$barcode_full), sub("^[^.]*", "", df$barcode_full), "")
  
  ct_col   <- intersect(ct_col_candidates,   names(df))
  cond_col <- intersect(cond_col_candidates, names(df))
  
  df$cond_std <- NA_character_
  if (length(cond_col)) {
    low <- tolower(as.character(df[[cond_col[1]]]))
    df$cond_std <- dplyr::case_when(
      low %in% c("resistant","drug","treated","res","r","resistance")      ~ "Resistant",
      low %in% c("sensitive","control","untreated","sen","s","sensitivity") ~ "Sensitive",
      TRUE ~ NA_character_
    )
  }
  list(meta = df, ct_col = ct_col[1], cond_col = if (length(cond_col)) cond_col[1] else NA_character_)
}

diagnose_barcode_overlap <- function(bc_sce, annot_meta, sid) {
  meta <- annot_meta
  ndot  <- sum(bc_sce %in% meta$barcode_pure_dot)
  ndash <- sum(bc_sce %in% meta$barcode_pure_dash)
  meta_sen <- meta[!is.na(meta$cond_std) & meta$cond_std == "Sensitive", ]
  meta_res <- meta[!is.na(meta$cond_std) & meta$cond_std == "Resistant", ]
  
  data.frame(
    sample = sid, n_cells = length(bc_sce),
    match_dot_total = ndot, match_dash_total = ndash,
    match_dot_Sensitive = sum(bc_sce %in% meta_sen$barcode_pure_dot),
    match_dot_Resistant = sum(bc_sce %in% meta_res$barcode_pure_dot),
    rate_dot  = sprintf("%.1f%%", 100*ndot/length(bc_sce))
  )
}

map_from_annot_to_sce <- function(sce, sid, annot_info) {
  key_sce <- colnames(sce)
  bc_sce  <- sub("^.*_", "", key_sce)
  meta    <- annot_info$meta
  ct_col  <- annot_info$ct_col
  
  meta_use_dot  <- meta[!duplicated(meta$barcode_pure_dot), ]
  meta_use_dash <- meta[!duplicated(meta$barcode_pure_dash), ]
  
  idx_dot  <- match(bc_sce, meta_use_dot$barcode_pure_dot)
  idx_dash <- match(bc_sce, meta_use_dash$barcode_pure_dash)
  
  if (sum(!is.na(idx_dot)) >= sum(!is.na(idx_dash)) && sum(!is.na(idx_dot)) > 0) {
    idx <- idx_dot; meta_chosen <- meta_use_dot; used_mode <- "pure(dot)"
  } else {
    idx <- idx_dash; meta_chosen <- meta_use_dash; used_mode <- "pure(dash)"
  }
  
  tibble::tibble(
    cell = key_sce, used_suffix = used_mode,
    seurat_barcode_full = meta_chosen$barcode_full[idx],
    sctype_class = meta_chosen[[ct_col]][idx],
    condition = meta_chosen$cond_std[idx],
    sample = sid, sce_barcode = bc_sce
  )
}

save_and_show <- function(p, out_prefix, outdir, layout = "auto", n_panels = NULL, ncol = NULL, w_percol = 4, h_perrow = 2.75, dpi = 300) {
  outdir <- ensure_writable_dir(outdir)
  if (layout == "single") { w <- 4; h <- 2.75 } 
  else if (layout == "two") { w <- 8; h <- 2.75 } 
  else if (layout == "2x2") { w <- 8; h <- 5.5 } 
  else {
    ncol <- if (is.null(ncol)) 2 else ncol
    w <- w_percol * ncol
    h <- h_perrow * ceiling(n_panels / ncol)
  }
  ggsave(file.path(outdir, paste0(out_prefix, ".pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, paste0(out_prefix, ".png")), p, width = w, height = h, dpi = dpi)
}
## =============================================================================
## SECTION 5: MAIN EXECUTION FLOW
## =============================================================================
set.seed(20250814)

annot_info <- get_annot_meta(load_annotation_seurat(annotation_seurat_path))
if (is.na(annot_info$ct_col)) stop("Annotation object has no recognised cell-type column.")

# Handle duplicate sample IDs by creating unique identifiers
if (any(duplicated(sample_tbl$sample_id))) {
  sample_tbl$sample_uid <- make.unique(sample_tbl$sample_id)
} else {
  sample_tbl$sample_uid <- sample_tbl$sample_id
}

# Initialize data containers
seu_list     <- list()
all_RE_tbl   <- list()
all_map_rows <- list()
diag_rows    <- list()

# Iterate through sample registry
for (i in seq_len(nrow(sample_tbl))) {
  sid_raw  <- sample_tbl$sample_id[i]
  sid      <- sample_tbl$sample_uid[i]
  tenx_dir <- .apa_path(sample_tbl$tenx_dir[i])
  txs_rds  <- .apa_path(sample_tbl$txs_rds[i])
  
  ## 5.1 Load Gene Expression Matrix (10x Genomics)
  # Standard Read10X call and barcode normalization
  mtx10x <- Read10X(data.dir = tenx_dir)
  colnames(mtx10x) <- norm_cells_10x(colnames(mtx10x), sid)
  
  # Deduplicate barcodes after normalization by aggregating counts
  if (any(duplicated(colnames(mtx10x)))) {
    mtx10x <- t(rowsum(t(mtx10x), group = colnames(mtx10x)))
  }
  
  ## 5.2 Load scUTRquant TXS results (SingleCellExperiment)
  # Ensure barcodes are synchronized with the expression matrix
  sce <- readRDS(txs_rds)
  txs_ids <- if (!is.null(colnames(sce))) colnames(sce) else as.character(colData(sce)$cell_id)
  colnames(sce) <- norm_cells_txs(txs_ids, sid)
  if ("cell_id" %in% colnames(colData(sce))) colData(sce)$cell_id <- colnames(sce)
  
  # Run barcode overlap diagnostics against reference annotation
  bc_sce <- sub("^.*_", "", colnames(sce))
  diag_rows[[sid]] <- diagnose_barcode_overlap(bc_sce, annot_info$meta, sid = sid)
  
  ## 5.3 Quantitative PAS Analysis
  # Calculate Relative Expression (RE) scores per gene and cell-wise mean RE
  txs_out <- compute_RE_from_txs_sce(sce, gff_path = gff_path)
  RE_mat  <- txs_out$RE
  RE_mean <- txs_out$RE_mean
  
  ## 5.4 Feature & Cell Synchronization
  # Intersection of barcodes present in both 10x and scUTRquant outputs
  common_cells <- intersect(colnames(mtx10x), colnames(RE_mat))
  if (length(common_cells) == 0) {
    stop(sprintf("Fatal Error: No overlapping cells found for sample [%s]", sid))
  }
  
  mtx10x  <- mtx10x[, common_cells, drop = FALSE]
  RE_mat  <- RE_mat[, common_cells, drop = FALSE]
  RE_mean <- RE_mean[common_cells]
  
  ## 5.5 Seurat Object Initialization
  # Standard creation with minimum cell/feature filtering
  seu <- CreateSeuratObject(counts = mtx10x, project = sid, min.cells = 3, min.features = 200)
  seu$sample  <- sid_raw
  seu$RE_mean <- as.numeric(RE_mean[colnames(seu)])
  
  ## 5.6 Metadata & Annotation Integration
  # Mapping pre-defined annotations (cell types, conditions) from reference Seurat
  map_df <- map_from_annot_to_sce(sce, sid, annot_info)
  all_map_rows[[sid]] <- map_df
  
  # Standardize cell type nomenclature (Handle NA/Unknowns)
  bads <- c("Unknown", "unknown", "UNK", "Unk", "UNKN", "Other", "")
  map_df$sctype_class[map_df$sctype_class %in% bads] <- NA
  
  # Assign metadata attributes to Seurat object
  v_ct   <- setNames(map_df$sctype_class, map_df$cell)
  v_cond <- setNames(map_df$condition,    map_df$cell)
  
  seu$cell_type_mapped <- unname(v_ct[colnames(seu)])
  seu$condition_mapped <- unname(v_cond[colnames(seu)])
  
  # Group definition: Prioritize Sensitive/Resistant status
  seu$group <- seu$condition_mapped
  seu$group[is.na(seu$group)] <- "Unknown"
  seu$group <- factor(seu$group, levels = c("Sensitive", "Resistant", "Unknown"))
  
  ## 5.7 Data Persistence
  # Export processed RE matrices and summary stats for downstream analysis
  saveRDS(RE_mat, file.path(OUTDIR, "RE_outputs", paste0(sid, "_RE_gene_by_cell.rds")))
  write.csv(
    data.frame(cell    = names(RE_mean), 
               RE_mean = as.numeric(RE_mean),
               sample  = sid_raw, 
               group   = as.character(seu$group[match(names(RE_mean), colnames(seu))])),
    file.path(OUTDIR, "RE_outputs", paste0(sid, "_RE_mean.csv")),
    row.names = FALSE
  )
  
  seu_list[[sid]]   <- seu
  all_RE_tbl[[sid]] <- data.frame(
    cell    = names(RE_mean), 
    RE_mean = as.numeric(RE_mean),
    sample  = sid_raw,
    group   = as.character(seu$group[match(names(RE_mean), colnames(seu))])
  )
}
## =============================================================================
## 6. BARCODE MAPPING & DIAGNOSTICS
## =============================================================================

# Export comprehensive mapping table for reproducibility
map_out_all <- dplyr::bind_rows(all_map_rows)
readr::write_csv(map_out_all, file.path(OUTDIR, "sce_barcode_to_seurat_mapping.csv"))

# Export diagnostic summary of barcode overlap rates
diag_tbl <- dplyr::bind_rows(diag_rows)
readr::write_csv(diag_tbl, file.path(OUTDIR, "mapping_diagnostics.csv"))

## =============================================================================
## 7. DATA INTEGRATION & PRE-PROCESSING
## =============================================================================

# Merge individual Seurat objects into a unified study object
seu_all <- Reduce(function(a, b) merge(a, b), seu_list)
RE_df   <- dplyr::bind_rows(all_RE_tbl)

# Synchronize RE scores into the merged meta.data
if (!"RE_mean" %in% colnames(seu_all@meta.data)) seu_all$RE_mean <- NA_real_
mt <- match(colnames(seu_all), RE_df$cell)
ok <- !is.na(mt)
seu_all$RE_mean[ok] <- RE_df$RE_mean[mt[ok]]

# Filter for cells with valid annotations and target conditions
label_col  <- "cell_type_mapped"
keep_cells <- !is.na(seu_all@meta.data[[label_col]]) & 
              seu_all$group %in% c("Sensitive", "Resistant")

seu_all <- subset(seu_all, cells = colnames(seu_all)[keep_cells])
seu_all$group <- factor(seu_all$group, levels = c("Sensitive", "Resistant"))

# Standard Seurat normalization and dimensionality reduction pipeline
seu_all <- NormalizeData(seu_all, verbose = FALSE)
seu_all <- FindVariableFeatures(seu_all, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu_all <- ScaleData(seu_all, features = rownames(seu_all), verbose = FALSE)
seu_all <- RunPCA(seu_all, features = VariableFeatures(seu_all), verbose = FALSE)
seu_all <- FindNeighbors(seu_all, dims = 1:20, verbose = FALSE)
seu_all <- FindClusters(seu_all, resolution = 0.8, verbose = FALSE)

## =============================================================================
## 8. RE ASSAY CONSTRUCTION
## =============================================================================

# Collect intermediate RE matrices from disk
re_files <- list.files(file.path(OUTDIR, "RE_outputs"), pattern = "_RE_gene_by_cell\\.rds$", full.names = TRUE)
re_list  <- lapply(re_files, readRDS)

# Subset to common features and synchronized cells
all_cells    <- colnames(seu_all)
re_list      <- lapply(re_list, function(M) M[, intersect(colnames(M), all_cells), drop = FALSE])
common_genes <- Reduce(intersect, lapply(re_list, rownames))
re_list      <- lapply(re_list, function(M) methods::as(M[common_genes, , drop = FALSE], "CsparseMatrix"))

# Build global RE matrix
RE_mat_tmp <- do.call(cbind, re_list)
RE_mat     <- RE_mat_tmp[, all_cells, drop = FALSE]

# Initialize RE assay within Seurat object
seu_all[["RE"]] <- CreateAssayObject(data = RE_mat)

## =============================================================================
## 9. CUMULATIVE DISTRIBUTION ANALYSIS (ECDF)
## =============================================================================

# Aesthetic mapping for conditions
cond_cols    <- c(Resistant = '#E4C66F', Sensitive = '#5B9BD5')
linetype_map <- c(Resistant = "solid",   Sensitive = "22")

# Prepare long-format data for ECDF plotting
meta_ecdf <- seu_all@meta.data %>%
  dplyr::transmute(
    cell      = rownames(.),
    cell_type = .data[[label_col]],
    group     = factor(group, levels = c("Sensitive", "Resistant")),
    RE_mean   = RE_mean
  ) %>% dplyr::filter(is.finite(RE_mean))

# Generate and save ECDF plots by cell type

if (nrow(meta_ecdf)) {
  p_ecdf <- ggplot(meta_ecdf, aes(x = RE_mean, colour = group, linetype = group)) +
    stat_ecdf(linewidth = 1.0) +
    facet_wrap(~ cell_type, scales = "free_y") +
    scale_color_manual(values = cond_cols) +
    scale_linetype_manual(values = linetype_map) +
    labs(x = "Preference of proximal PAS (RE_mean = log2(aUTR/cUTR))", y = "Cumulative ratio") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "top")
  
  save_and_show(p_ecdf, "RE_mean_ECDF_by_CellType", file.path(OUTDIR, "ECDF_by_celltype"), 
                n_panels = length(unique(meta_ecdf$cell_type)))
}

## =============================================================================
## 10. DIFFERENTIAL APA ANALYSIS & VIOLIN PLOTS
## =============================================================================

celltypes <- sort(na.omit(unique(seu_all@meta.data[[label_col]])))
DE_results <- list()

for (ct in celltypes) {
  # Subset to specific cell type
  obj_ct <- subset(seu_all, cells = rownames(seu_all@meta.data)[seu_all@meta.data[[label_col]] == ct])
  
  # Differential testing using Wilcoxon Rank Sum test on RE scores
  # FindMarkers used for Resistant vs Sensitive comparison
  de <- tryCatch({
    FindMarkers(obj_ct, ident.1 = "Resistant", ident.2 = "Sensitive", 
                assay = "RE", test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
  }, error = function(e) NULL)
  
  if (is.null(de) || nrow(de) == 0) next
  de$gene      <- rownames(de)
  de$cell_type <- ct
  DE_results[[ct]] <- de

  # Visualize Top 4 differentially used PAS
  top_genes <- de %>% dplyr::slice_max(order_by = abs(avg_log2FC), n = 4) %>% pull(gene)
  
  
  p_vln <- VlnPlot(obj_ct, features = top_genes, group.by = "group", cols = cond_cols, pt.size = 0.1) & 
           theme_bw() & labs(y = "RE Score")
  
  save_and_show(p_vln, paste0(gsub("[^A-Za-z0-9]+", "_", ct), "_RE_top_markers"), 
                file.path(OUTDIR, "violin_by_celltype"), layout = "2x2")
}

# Export consolidated DE summary
all_de <- dplyr::bind_rows(DE_results)
readr::write_csv(all_de, file.path(OUTDIR, "RE_Differential_Analysis_Results.csv"))