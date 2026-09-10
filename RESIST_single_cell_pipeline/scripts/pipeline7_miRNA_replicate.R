library(biomaRt)
library(ggplot2)
library(dplyr)
library(tibble)
miRDB <- read.delim("/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/miRDB_v6.0_prediction_result.txt", header = F)
colnames(miRDB)[1]<- 'miRNA'
colnames(miRDB)[2]<- 'TargetID'
colnames(miRDB)[3]<- 'Score'

#mart_h <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#mart_m <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

mart_h<- readRDS('/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/mart_h.rds')
mart_m<- readRDS('/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/mart_m.rds')
all_mappings_h<- read.csv('/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/all_mappings_h.csv')
all_mappings_m<- read.csv( '/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/all_mappings_m.csv')

deg_dir <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/results"

out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
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



