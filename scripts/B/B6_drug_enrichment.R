#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B6_drug_enrichment.R
# -----------------------------------------------------------------------------
# Purpose  : Query the resistance signature of each dataset against LINCS L1000 Level-4
#            profiles to rank compounds whose transcriptional response opposes the
#            resistant state - the candidate resistance-reversal agents surfaced by
#            the RESIST browser.
# Inputs   : data/results/B/*.sig.csv (from B1)
#            data/reference/GSE70138_LINCS_Level4.rds
# Outputs  : data/results/B/*_drug_enrichment.csv
# Usage    : Rscript scripts/B/B6_drug_enrichment.R
# Origin   : scripts/pipeline9_drug_enrichment_replicate.R
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

#'--------------------
#' set up
#' --------------------

library(tidyverse)
library(devtools)
#install_github("RGLab/cytolib")
library(cytolib)
#BiocManager::install("flowCore")
library(flowCore)
#BiocManager::install("cmapR")
library(cmapR)

#gct_file<- './GSE70138_Broad_LINCS_Level4_ZSPCINF_mlr12k_n78980x22268_2015-06-30.gct'
#gct_data<- parse_gctx(gct_file)

#saveRDS(gct_data, file = "GSE70138_LINCS_Level4.rds") ## save to rds for future use, it is too large, 13G and need very large RAM 50G and larger later to analyze near 70G or larger 
gct_data<- readRDS(resist_ref("lincs_level4_rds"))
# Access matrix
expr_matrix <- gct_data@mat

# Access gene info
gene_metadata <- gct_data@rdesc

# Access sample info
sample_metadata <- gct_data@cdesc


# Get all unique compound names (excluding DMSO)
drug_names <- unique(sample_metadata$SM_Name)
drug_names <- drug_names[drug_names != "DMSO"]




#'--------------------------------------
#' Define the KS-based enrichment functions
#'--------------------------------------
library(fgsea)

compute_es <- function(deg_set, drug_profile, type) {
  r <- length(drug_profile)  # Total genes in profile (~12,328)
  s <- length(deg_set)       # Number of DEGs
  
  # Get ranks of DEGs in drug profile (ascending: 1 = most downregulated)
  v_j <- rank(drug_profile)[deg_set] %>% sort()
  
  # Compute KS statistics
  j <- 1:s
  a <- max(j/s - v_j/r)
  b <- max(v_j/r - (j-1)/s)
  
  if (type == "up") {
    return(ifelse(a > b, a, -b))  # ES_up
  } else {
    return(ifelse(a > b, a, -b))  # ES_down
  }
}



#'--------------------------------------------------------------------
#' loop through all drugs for each dataset to get enrichment score
#' -------------------------------------------------------------------

outdir <- RESULTS_B
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

#' my deg, List all CSV files in the folder
deg_files <- list.files(RESULTS_B, pattern = '\\.sig\\.csv$', full.names = TRUE)

# Loop through each file
for (file in deg_files) {
  cat("Reading:", file, "\n")  # optional: print filename
  file_name <- basename(file)

  my_deg <- read.csv(file)
  cell_types<- my_deg$cell_type
  if (!'Tumor cells' %in% cell_types) next
  
  my_deg<- my_deg %>% dplyr::filter(cell_type=='Tumor cells')
  up_deg<- my_deg %>% dplyr::filter(avg_log2FC>0) %>% pull(gene)
  dn_deg<- my_deg %>% dplyr::filter(avg_log2FC<0) %>% pull(gene)



  # Initialize results data frame
  results <- data.frame(
    pert_id = character(),
    es_up = numeric(),
    es_down = numeric(),
    es_final = numeric(),
    pval = numeric(),
    fdr = numeric(),
    stringsAsFactors = FALSE
  )



  for (drug in drug_names) {
    filtered_samples <- sample_metadata$SM_Name %in% c('DMSO', drug) &
    sample_metadata$SM_Pert_Type %in% c("trt_cp", "ctl_vehicle") & sample_metadata$SM_Time == '24'

    # Get the corresponding column names (sample IDs)
    selected_sample_ids <- rownames(sample_metadata)[filtered_samples]

    # Subset the expression matrix
    expr_subset <- expr_matrix[, selected_sample_ids]

    # subset sample metadata
    sample_metadata_subset <- sample_metadata[filtered_samples, ]
    identical(rownames(expr_subset), rownames(gene_metadata))
    rownames(expr_subset)<- gene_metadata$pr_gene_symbol



    # Calculate mean z-score per gene for the drug
    gene_zscore_mean <- rowMeans(expr_subset, na.rm = TRUE)
    # Assign gene symbols as names
    gene_symbols <- rownames(expr_subset)
    names(gene_zscore_mean) <- gene_symbols
    # Remove "-666" and NA gene names
    valid_genes <- !is.na(names(gene_zscore_mean)) & names(gene_zscore_mean) != "-666"
    gene_zscore_mean <- gene_zscore_mean[valid_genes]
    # Sort by z-score ascending (most downregulated first)
    gene_zscore_sorted <- sort(gene_zscore_mean, decreasing = FALSE)
    # Keep only the first occurrence of each gene symbol
    drug_profile <- gene_zscore_sorted[!duplicated(names(gene_zscore_sorted))]


    if (length(intersect(names(drug_profile), up_deg)) < 10 |
      length(intersect(names(drug_profile), dn_deg)) < 10 |
      length(intersect(names(drug_profile), union(up_deg, dn_deg))) >= 2000
      ) next


    # Compute ES_up and ES_down
    es_up <- compute_es(intersect(names(drug_profile), up_deg), drug_profile, "up")
    es_down <- compute_es(intersect(names(drug_profile), dn_deg), drug_profile, "down")



    # Combine scores
    if (sign(es_up) != sign(es_down)) {
      es_final <- es_up - es_down
    } else {
      es_final <- 0
    }


    # permutation test - shuffle DEG sets
    n_perm <- 1000
    null_es <- numeric(n_perm)
    n_up <- length(up_deg)
    n_down <- length(dn_deg)
    all_genes <- names(drug_profile)

    for(p in 1:n_perm) {
      # Random gene sets
      rand_up <- sample(all_genes, n_up)
      rand_down <- sample(all_genes, n_down)
    
      # Compute null ES 
      es_up_rand <- compute_es(intersect(names(drug_profile), rand_up), drug_profile, "up")
      es_down_rand <- compute_es(intersect(names(drug_profile), rand_down), drug_profile, "down")
    
      if(sign(es_up_rand) != sign(es_down_rand)) {
        null_es[p] <- es_up_rand - es_down_rand
      } else {
        null_es[p] <- 0
    }
  }

    # Empirical p-value
    if (es_final > 0) {
      pval <- mean(null_es >= es_final)
    } else if (es_final < 0) {
      pval <- mean(null_es <= es_final)
    } else {
      pval <- 1
    }


    # Add result for this drug
    results <- rbind(results, data.frame(
      pert_id = drug,
      es_up = es_up,
      es_down = es_down,
      es_final = es_final,
      pval = pval,
      fdr = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate FDR (after all drugs)
  results$fdr <- p.adjust(results$pval, method = "fdr")

  # filter sig
  results<- results %>% dplyr::filter(fdr<0.05) %>% arrange (-es_final)

  # save results for each dataset's deg
  path_sig <- file.path(outdir, paste0(file_name, "_drug_enrichment.csv"))
  write.csv(results, file = path_sig, row.names = FALSE)
  
  print(paste0(file_name, ' is done!'))

}