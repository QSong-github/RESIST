#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B2_pathway_enrichment.R
# -----------------------------------------------------------------------------
# Purpose  : Map the differentially expressed genes of each dataset and cell type onto
#            GO and pathway gene sets, separating programmes that are induced in the
#            resistant state from those that are lost.
# Inputs   : data/results/B/*_deg.csv (from B1)
# Outputs  : data/results/B/*_GO_enrichment_up.csv
#            data/results/B/*_GO_enrichment_down.csv
# Usage    : Rscript scripts/B/B2_pathway_enrichment.R
# Origin   : scripts/pipeline3_enrichment_replicate.R
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

library(ggplot2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)

deg_dir <- RESULTS_B
out_dir <- RESULTS_B
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
deg_files <- list.files(deg_dir, pattern = "deg\\.sig\\.csv$", full.names = TRUE)


#'----------------------
#' GO & KEGG
#' ----------------------
for (file in deg_files) { ## remove the [17:18] to run all files
  deg_all <- read.csv(file)
  base_name <- tools::file_path_sans_ext(basename(file))
  cell_types <- unique(deg_all$cell_type)
  
  if (!"Tumor cells" %in% cell_types) next
  
  deg <- deg_all %>% filter(cell_type == "Tumor cells")
  rownames(deg) <- NULL
  deg <- deg %>% column_to_rownames(var = "gene")
  
  up <- deg %>% filter(p_val_adj < 0.05, avg_log2FC >= 1)
  dn <- deg %>% filter(p_val_adj < 0.05, avg_log2FC <= -1)
  
  
#  first_gene <- rownames(dn)[1]
#  if (grepl("^[A-Z0-9]+$", first_gene)) {
    org <- org.Hs.eg.db
    orga<- 'hsa'
 # } else {
 #   org <- org.Mm.eg.db
 #   orga<- 'mmu'
 # }
    
  
  dn.gene<- bitr(rownames(dn), fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org)
  up.gene<- bitr(rownames(up), fromType = "SYMBOL", 
                 toType = "ENTREZID",
                 OrgDb = org)
  
  #------------ GO ------------#
  up.go <- enrichGO(
    gene           = up.gene$ENTREZID,
    OrgDb          = org,
    ont            = "BP",
    pAdjustMethod  = "BH"
  )
  dn.go <- enrichGO(
    gene           = dn.gene$ENTREZID,
    OrgDb          = org,
    ont            = "BP",
    pAdjustMethod  = "BH"
  )
  
  up.go.res<- up.go@result
  dn.go.res<- dn.go@result
  
  
  # Output
  path_up <- file.path(out_dir, paste0(base_name, "_GO_enrichment_up.csv"))
  write.csv(up.go.res, path_up)
  
  path_dn <- file.path(out_dir, paste0(base_name, "_GO_enrichment_dn.csv"))
  write.csv(dn.go.res, path_dn)
  
  #------------ KEGG ------------#
  
  up.kegg <- enrichKEGG(gene = up.gene$ENTREZID, organism = orga)
  dn.kegg <- enrichKEGG(gene = dn.gene$ENTREZID, organism = orga)
  
  up.kegg.res<- up.kegg@result
  dn.kegg.res<- dn.kegg@result

  # Output
  path_up <- file.path(out_dir, paste0(base_name, "_KEGG_enrichment_up.csv"))
  write.csv(up.kegg.res, path_up)
  
  path_dn <- file.path(out_dir, paste0(base_name, "_KEGG_enrichment_dn.csv"))
  write.csv(dn.kegg.res, path_dn)
  
  message(paste0(base_name, " enrichment files saved."))
}


#'----------------------
#' Hallmark
#' ----------------------
for (file in deg_files) { ## remove the [17:18] to run all files
  deg_all <- read.csv(file)
  base_name <- tools::file_path_sans_ext(basename(file))
  cell_types <- unique(deg_all$cell_type)
  
  if (!"Tumor cells" %in% cell_types) next
  
  deg <- deg_all %>% filter(cell_type == "Tumor cells")
  rownames(deg) <- NULL
  deg <- deg %>% column_to_rownames(var = "gene")
  
  up <- deg %>% filter(p_val_adj < 0.05, avg_log2FC >= 1)
  dn <- deg %>% filter(p_val_adj < 0.05, avg_log2FC <= -1)
  
  
 # first_gene <- rownames(dn)[1]
 # if (grepl("^[A-Z0-9]+$", first_gene)) {
    spe = "Homo sapiens"
    org <- org.Hs.eg.db
 # } else {
 #   spe = "Mus musculus"
  #  org <- org.Mm.eg.db
 # }
  
  
  dn.gene<- bitr(rownames(dn), fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org)
  up.gene<- bitr(rownames(up), fromType = "SYMBOL", 
                 toType = "ENTREZID",
                 OrgDb = org)
  
  #------------ hallmark ------------#
  
  m_t2g <- msigdbr(species = spe, category = 'H') %>% 
    dplyr::select(gs_name, entrez_gene)
  
  up.em<- enricher(up.gene$ENTREZID,TERM2GENE = m_t2g)
  dn.em<- enricher(dn.gene$ENTREZID,TERM2GENE = m_t2g)
  
  up.em.res<- up.em@result
  dn.em.res<- dn.em@result
  
  
  # Output
  path_up <- file.path(out_dir, paste0(base_name, "_Hallmark_enrichment_up.csv"))
  write.csv(up.em.res, path_up)
  
  path_dn <- file.path(out_dir, paste0(base_name, "_Hallmark_enrichment_dn.csv"))
  write.csv(dn.em.res, path_dn)
  
  message(paste0(base_name, " enrichment files saved."))
  
}
