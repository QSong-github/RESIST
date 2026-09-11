#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C9_mirna_target_enrichment.R
# -----------------------------------------------------------------------------
# Purpose  : Test the resistance-associated gene sets for over-representation of miRDB
#            predicted targets, nominating microRNAs that plausibly drive the
#            transcriptional change.
# Inputs   : data/results/B/*_deg.csv (from B1)
#            data/reference/miRDB_v6.0_prediction_result.txt
#            data/reference/mart_{h,m}.rds, data/reference/all_mappings_{h,m}.csv
# Outputs  : data/results/C/*_mirna_enrichment.csv
# Usage    : Rscript scripts/C/C9_mirna_target_enrichment.R
# Origin   : scripts/pipeline7_miRNA_replicate.R
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

library(biomaRt)
library(ggplot2)
library(dplyr)
library(tibble)
miRDB <- read.delim(resist_ref("mirdb"), header = F)
colnames(miRDB)[1]<- 'miRNA'
colnames(miRDB)[2]<- 'TargetID'
colnames(miRDB)[3]<- 'Score'

#mart_h <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#mart_m <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

mart_h<- readRDS(resist_ref("biomart_human"))
mart_m<- readRDS(resist_ref("biomart_mouse"))
all_mappings_h<- read.csv(resist_ref("gene_map_human"))
all_mappings_m<- read.csv( resist_ref("gene_map_mouse"))

deg_dir <- RESULTS_B

out_dir <- RESULTS_C
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

deg_files <- list.files(deg_dir, pattern = "deg\\.sig\\.csv$", full.names = TRUE)


# these are mouse 
mouse_dat<- c('GSE201765_deg.sig',
'GSE213183_deg.sig',
'GSE261898_post_deg.sig',
'GSE261898_pre_deg.sig')



#'-----------------------
#' human
#'---------------------

for (file in deg_files) { # remove the [17:18] to run all files
  tryCatch({
    
    deg_all <- read.csv(file)
    base_name <- tools::file_path_sans_ext(basename(file))
    
    if (base_name %in% mouse_dat) next
    
    
    cell_types <- unique(deg_all$cell_type)
    if (!"Tumor cells" %in% cell_types) next
    
    deg <- deg_all %>% filter(cell_type == "Tumor cells")
    up <- deg %>% filter(p_val_adj < 0.05, avg_log2FC >= 1) %>%  
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = 10)
    dn <- deg %>% filter(p_val_adj < 0.05, avg_log2FC <= -1) %>%
      arrange(desc(abs(avg_log2FC))) %>%
      slice_head(n = 10)
    
    up.gene <- up$gene
    dn.gene <- dn$gene
    
    
    # UP------------------------------------------------------------------------------------------------------
    # Map Gene Symbols to RefSeq IDs
    mapped_up <- subset(all_mappings_h, hgnc_symbol %in% up.gene)
    mapped_up<- mapped_up[!is.na(mapped_up$refseq_mrna) & mapped_up$refseq_mrna != "", ]
    mapped_up<- mapped_up[!duplicated(mapped_up$hgnc_symbol), ]
    
    # Query miRDB for Each Transcript
    filtered_up <- miRDB %>%
      filter(TargetID %in% mapped_up$refseq_mrna & Score > 80) %>%
      left_join(mapped_up, by = c("TargetID" = "refseq_mrna")) %>% dplyr::select(hgnc_symbol, GeneBank_accession = TargetID, miRNA, Score) %>%
      arrange(hgnc_symbol)
    
    # Save output
    path1 <- file.path(out_dir, paste0(base_name, "_up_miRNA.csv"))
    write.csv(filtered_up, file = path1, row.names = FALSE)
    
    # DN------------------------------------------------------------------------------------------------------
    mapped_dn <- subset(all_mappings_h, hgnc_symbol %in% dn.gene)
    
    if (nrow(mapped_dn) > 0) {
      mapped_dn <- mapped_dn[!duplicated(mapped_dn$hgnc_symbol), ]
      mapped_dn<- mapped_dn[!is.na(mapped_dn$refseq_mrna) & mapped_dn$refseq_mrna != "", ]
      filtered_dn <- miRDB %>%
        filter(TargetID %in% mapped_dn$refseq_mrna & Score > 80) %>%
        left_join(mapped_dn, by = c("TargetID" = "refseq_mrna")) %>% dplyr::select(hgnc_symbol, GeneBank_accession = TargetID, miRNA, Score) %>%
        arrange(hgnc_symbol)
      
      path2 <- file.path(out_dir, paste0(base_name, "_dn_miRNA.csv"))
      write.csv(filtered_dn, file = path2, row.names = FALSE)
    }
    message(paste0(base_name, " - miRNA saved."))
  }, error = function(e) {
    message(paste0(file, " - skipped due to error: ", conditionMessage(e)))
  })
}



