#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B5_emt_score.R
# -----------------------------------------------------------------------------
# Purpose  : Score every tumour cell against the MSigDB Hallmark EMT signature and
#            compare the distribution between sensitive and resistant states, one of
#            the recurrent transcriptional programmes of acquired resistance.
# Inputs   : data/*.rds - annotated Seurat objects
#            data/reference/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt
# Outputs  : data/results/B/<input_filename>_EMT_box.{pdf,png}
# Usage    : Rscript scripts/B/B5_emt_score.R
# Origin   : scripts/pipeline5_EMT_replicate.R
# Notes    : This entry point inherited a different hold-out list from the other
#            Module A/B scripts (`paired_pre_post_anno`). The discrepancy is
#            preserved deliberately - see docs/repository-guide.md, Known issues.
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

library(GSVA)
library(GSEABase)
library(tidyverse)
library(Seurat)

gmt_human <- resist_ref("emt_gmt_human")
gene_sets_human <- getGmt(gmt_human)
gene_sets_human<- geneIds(gene_sets_human)



# these are the datasets that have pre and post samples
ids_timppint <- resist_datasets("paired_pre_post_anno")
#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

out_dir <- RESULTS_B
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
#rds_files <- list.files('../DRM_data', full.names = TRUE)
rds_files <- list.files(PATH_DATA, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)

# Loop through each file
for (file_path in rds_files) { ## remove the [10:33] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_.*", "", file_name)
  
  if (file_name %in% ids_timppint) {
    next
  }
  
  seu <- resist_read_seurat(file_path)
  celltypes <- unique(seu$cell_type)
  if (!"Tumor cells" %in% celltypes) next
  
  malignant_cells <- subset(seu, subset = cell_type == "Tumor cells")
  
  #first_gene <- rownames(malignant_cells)[5]
 # if (grepl("^[A-Z]{2}", first_gene)) {
    geneset <- gene_sets_human
 # } else {
  #  geneset <- gene_sets_mouse
 # }
  
  expr_mat <- as.matrix(GetAssayData(malignant_cells, slot = "data"))  # genes x cells
  
  gsvaPar<- gsvaParam(exprData = expr_mat, geneset) 
  gsva.res <- gsva(gsvaPar, verbose=FALSE)
  emt_scores <- as.numeric(gsva.res[, colnames(malignant_cells)])

  
  #  Add EMT score to metadata
  malignant_cells$EMT <- emt_scores
  
  
  # draw box plot
  df<- malignant_cells@meta.data
  df$condition<- factor(df$condition, levels = c('sensitive' , 'resistant'))
  
  
  # Calculate Tukey-based limits
  q1 <- quantile(df$EMT, 0.25, na.rm = TRUE)
  q3 <- quantile(df$EMT, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr + 0.1
  
  # pval
  wilcox_res <- wilcox.test(EMT ~ condition, data = df)
  
  p<- ggplot(df, aes(x = condition, y = EMT, fill = condition)) +
    geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +  # removes outliers
    stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +   # adds whisker caps
    labs(title = "EMT score", y = "EMT score", x= NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14, face = 'bold', color = 'black'), 
      axis.text.y = element_text(size = 14, face = 'bold', color = 'black'), 
      axis.title.y = element_text(size = 14, face = 'bold', color = 'black'),
      plot.title = element_text(face='bold', size=14, hjust = 0.5),
      legend.position = 'none'
    ) +
    scale_fill_manual(values = rev(col1))+
    annotate(
      "text",
      x = 1.5,
      y = upper_bound,
      label = paste0("p = ", signif(wilcox_res$p.value, 3)),
      fontface = "bold",
      size = 5
    )+
    coord_cartesian(ylim = c(NA, upper_bound))
  
  
  # Output paths
  pdf_path <- file.path(out_dir, paste0(file_name, "_EMT_box.pdf"))
  png_path <- file.path(out_dir, paste0(file_name, "_EMT_box.png"))
  
  # Save PDF
  pdf(pdf_path, width = 4, height = 2.75)
  print(p)
  dev.off()
  
  # Save PNG
  png(png_path, width = 4, height = 2.75, units = "in", res = 300)
  print(p)
  dev.off()
  
  message(paste0(file_name, " - EMT saved."))
}


