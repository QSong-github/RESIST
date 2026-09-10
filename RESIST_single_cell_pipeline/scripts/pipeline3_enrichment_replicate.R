library(ggplot2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
source('/blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/color.R')

deg_dir <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/results"
out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
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
