#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B1_differential_expression.R
# -----------------------------------------------------------------------------
# Purpose  : Test every cell type for genes that change between the drug-sensitive
#            and drug-resistant state. This is the upstream step of Module B and of
#            the miRNA and RBP analyses in Module C - they all consume the DEG tables
#            written here.
# Inputs   : paths.data/*.rds or *.RDS - annotated Seurat objects
# Outputs  : data/results/B/<dataset>_deg.csv - full statistics
#            data/results/B/<dataset>_deg.sig.csv - significant subset
#            data/results/B/<dataset>_deg_Tumor_cells_volcano.{pdf,png}
# Usage    : Rscript scripts/B/B1_differential_expression.R
# Origin   : scripts/pipeline2_DEG_replicate.R
# Notes    : Harmony-integrated datasets are held out of the default contrast;
#            see `datasets.harmony_integrated` in config/config.yaml.
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

# these are the datasets that have pre and post samples
ids_timppint <- resist_datasets("paired_pre_post_seu")
# these are datasets with harmony
ids_harmony <- resist_datasets("harmony_integrated")

#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

outdir <- RESULTS_B

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

rds_files <- list.files(PATH_DATA, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)


# Loop through each file
for (file_path in rds_files) { ## remove the [12:13] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_seu.*", "", file_name)
  
  if (file_name %in% ids_timppint || file_name %in% ids_harmony) {
    next
  }
  
  seu <- resist_read_seurat(file_path)
  # DEG between resistant and sensitive in each cell type

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

deg_dir <- RESULTS_B
out_dir <- RESULTS_B

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