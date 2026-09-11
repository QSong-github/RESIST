#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C7_rbp_target_enrichment.R
# -----------------------------------------------------------------------------
# Purpose  : Test whether the genes that change with resistance are enriched among the
#            target sets of individual RNA-binding proteins, implicating
#            post-transcriptional regulators in the resistant phenotype.
# Inputs   : data/results/B/*_deg.csv (from B1)
#            data/reference/All_RBP_TargetGene/ - one target list per RBP
# Outputs  : data/results/C/*_deg_rbp_enrichment.csv
# Usage    : Rscript scripts/C/C7_rbp_target_enrichment.R
# Origin   : scripts/pipeline8_RBP1_enrichment_replicate.R
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

library(dplyr)
library(readr)
library(tibble)

out_dir <- RESULTS_C
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

rbp_files <- list.files(resist_ref("rbp_target_dir"), pattern = "\\.txt$", full.names = TRUE) # Get list of RBP files
deg_files <- list.files(RESULTS_B, pattern = "deg\\.csv$", full.names = TRUE) # Get list of DEG files

# these are mouse 
mouse_dat<- c('GSE201765_deg',
              'GSE213183_deg',
              'GSE261898_post_deg',
              'GSE261898_pre_deg')



#'-----------------
#' human
#'-----------------

# Loop through each DEG file
for (deg_path in deg_files) { # remove the [17:18] to run all files
  # Load DEG data
  deg_all <- read.csv(deg_path)
  deg_celltypes <- unique(deg_all$cell_type)
  deg_file_short <- tools::file_path_sans_ext(basename(deg_path))
  
  if (!"Tumor cells" %in% deg_celltypes) next
  if (deg_file_short %in% mouse_dat) next
  
  deg <- deg_all %>% filter(cell_type == "Tumor cells")
  rownames(deg) <- NULL
  deg <- deg %>% column_to_rownames(var = "gene")
  
  # Define universe
  all_genes <- rownames(deg)
  N <- length(all_genes)
  
  # Define DEGs
  up <- deg %>% filter(p_val_adj < 0.05, avg_log2FC >= 1)
  dn <- deg %>% filter(p_val_adj < 0.05, avg_log2FC <= -1)
  up_genes <- rownames(up)
  dn_genes <- rownames(dn)
  
  # Result table
  result_list <- list()
  
  for (rbp_path in rbp_files) {
    # Get RBP name from filename
    rbp_name <- tools::file_path_sans_ext(basename(rbp_path))
    
    # Read RBP targets
    rbp_data <- tryCatch(
      read.delim(rbp_path, skip = 3, header = TRUE),
      error = function(e) return(NULL)
    )
    if (is.null(rbp_data) || !"geneName" %in% colnames(rbp_data)) next
    
    targets <- unique(as.character(rbp_data$geneName))
    M <- sum(all_genes %in% targets) #' M: total number of ACIN1 target genes in the universe
    
    # Hypergeometric for UP genes
    K_up <- length(up_genes)
    x_up <- sum(up_genes %in% targets)
    p_up <- phyper(q = x_up - 1, m = M, n = N - M, k = K_up, lower.tail = FALSE)
    
    # Hypergeometric for DOWN genes
    K_dn <- length(dn_genes)
    x_dn <- sum(dn_genes %in% targets)
    p_dn <- phyper(q = x_dn - 1, m = M, n = N - M, k = K_dn, lower.tail = FALSE)
    
    result_list[[rbp_name]] <- data.frame(
      RBP = rbp_name,
      Targets = M, #' M: total number of ACIN1 target genes in the universe
      Up_DEG = K_up, #' number of up-regulated genes
      Up_overlap = x_up, #' number of overlapping genes between targets and up-regulated DEGs
      Up_pval = p_up,
      
      Down_DEG = K_dn, #' number of dn-regulated genes
      Down_overlap = x_dn, #' number of overlapping genes between targets and dn-regulated DEGs
      Down_pval = p_dn,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  results_df <- bind_rows(result_list)
  results_df$FDR_up<- p.adjust(results_df$Up_pval, method = "holm")
  results_df$FDR_dn<- p.adjust(results_df$Down_pval, method = "holm")
  
  # Output name based on input DEG file
  write.csv(results_df, file = file.path(RESULTS_C, paste0(deg_file_short, "_rbp_enrichment.csv")), row.names = FALSE)
}



