# these are the datasets that have pre and post samples
ids_timppint<- c('GSE120575_anti_PD1_seu.rds',
                 'GSE120575_anti_CTLA4_PD1_seu.rds',
                 'GSE123813_bcc_pembrolizumab_seu.rds',
                 'GSE123813_scc_pembrolizumab_seu.rds',
                 'GSE123813_scc_cemiplimab_seu.rds',
                 'GSE261898_seu.rds')


# these are datasets with harmony
ids_harmony<- c('GSE223779_crizotinib_seu.rds',
                'GSE223779_alectinib_seu.rds',
                'GSE233766_1uPLX_seu.rds')


source('/blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/color.R')  

#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

outdir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

rds_files <- list.files('/blue/qsong1/sen.guo/resist_02/resist_test_share/data', pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)


# Loop through each file
for (file_path in rds_files) { ## remove the [12:13] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_seu.*", "", file_name)
  
  if (file_name %in% ids_timppint || file_name %in% ids_harmony) {
    next
  }
  
  seu <- readRDS(file_path)
  # DEG between resistant and sensitive in each cell type
  seu@meta.data$celltype[seu@meta.data$celltype == "Malignant cells"] <- "Tumor cells" ### the input rds data use "Malignant cells"
  colnames(seu@meta.data)[colnames(seu@meta.data) == "celltype"] <- "cell_type" ## change the column name to "cell_type" for consistency

  cell_types <- unique(seu$cell_type) 
  deg_list <- list()
   for (ct in cell_types) {
    subset_seu <- subset(seu, subset = cell_type == ct)
    if (length(unique(subset_seu$condition)) < 2) next
    cond_counts <- table(subset_seu$condition)
    if (any(cond_counts < 3)) next
  
    Idents(subset_seu) <- "condition"
    
    deg <- FindMarkers(subset_seu, ident.1 = "resistant", ident.2 = "sensitive")
    
    deg_list[[ct]] <- deg
   }
  
  deg_combined <- do.call(rbind, lapply(names(deg_list), function(ct) {
    df <- deg_list[[ct]]
    df$cell_type <- ct
    df$gene <- rownames(df)
    return(df)
  }))
  
  path <- file.path(outdir, paste0(dataset_id, "_deg.csv"))
  write.csv(deg_combined, file = path, row.names = FALSE)
  
  deg_combined_sig<- deg_combined %>% filter(p_val_adj<0.05, abs(avg_log2FC) >= 1)
  path_sig <- file.path(outdir, paste0(dataset_id, "_deg.sig.csv"))
  write.csv(deg_combined_sig, file = path_sig, row.names = FALSE)
  
  print(paste0(dataset_id, ' is done!'))
}
  
 
#'-------------------------------------------------------------------------
#' part 4: volcano
#'-------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tibble)

deg_dir <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/results"
out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
deg_files <- list.files(deg_dir, pattern = "deg\\.csv$", full.names = TRUE)



for (file in deg_files) { ## remove the [17:18] to run all files
  deg_all <- read.csv(file)
  
  base_name <- tools::file_path_sans_ext(basename(file))
  cell_types <- unique(deg_all$cell_type)
  
  if (!"Tumor cells" %in% cell_types) next
  
  deg <- deg_all %>% filter(cell_type == "Tumor cells")
  rownames(deg) <- NULL
  
  deg <- deg %>% column_to_rownames(var = "gene")
  
  deg.sig <- deg %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 1)
  
  # Select top 10 up and down genes
  top_up <- deg.sig %>% filter(avg_log2FC > 0) %>% arrange(-avg_log2FC) %>% slice_head(n = 10)
  top_down <- deg.sig %>% filter(avg_log2FC < 0) %>% arrange(avg_log2FC) %>% slice_head(n = 10)
  top_genes <- bind_rows(top_up, top_down)
  
  # Volcano base
  p <- ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(shape = 19, color = 'grey80', size = 1) +
    geom_point(data = deg.sig, aes(x = avg_log2FC, y = -log10(p_val_adj)), color = "red", size = 1) +
    labs(x = "log2FoldChange", y = "-log10(padj)", title = "Tumor cells") +
    theme_classic() +
    theme(
      legend.position = 'none',
      plot.margin = margin(20, 20, 20, 20),
      axis.text = element_text(size = 10, face = 'bold', color = 'black'),
      axis.title = element_text(size = 11, face = 'bold'),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  
  # Labeled plot
  p_labeled <- p +
    geom_text(data = top_genes, aes(label = rownames(top_genes)),
              check_overlap = TRUE, fontface = 'bold.italic', hjust = 1, vjust = -0.5, size = 3)
  
  # Output paths
  pdf_path <- file.path(out_dir, paste0(base_name, "_Tumor_cells_volcano.pdf"))
  png_path <- file.path(out_dir, paste0(base_name, "_Tumor_cells_volcano.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 5.5)
  print(p_labeled)
  dev.off()
  
  # Save PNG
  png(png_path, width = 8, height = 5.5, units = "in", res = 300)
  print(p_labeled)
  dev.off()
  
  message(paste0(base_name, " - Tumor cells volcano saved."))
}


#'-------------------------------------------------------------------------