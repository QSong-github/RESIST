#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B4_intratumor_heterogeneity.R
# -----------------------------------------------------------------------------
# Purpose  : Quantify transcriptional heterogeneity within the tumour compartment and
#            test whether resistant tumours are more heterogeneous than their
#            treatment-naive counterparts.
# Inputs   : data/*.rds - annotated Seurat objects
# Outputs  : data/results/B/<input_filename>_ITH_box.{pdf,png}
# Usage    : Rscript scripts/B/B4_intratumor_heterogeneity.R
# Origin   : scripts/pipeline4_ITH_replicate.R
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

library(Seurat)
library(tidyverse)
library(patchwork)

# these are the datasets that have pre and post samples
ids_timppint <- resist_datasets("paired_pre_post_seu")
#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

out_dir <- RESULTS_B
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
} 
rds_files <- list.files(PATH_DATA, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)


# Loop through each file
for (file_path in rds_files) { # remove the [13] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_.*", "", file_name)
  
  if (file_name %in% ids_timppint) {
    next
  }
  
  seu <- resist_read_seurat(file_path)
  cell_types <- unique(seu$cell_type)
  if (!"Tumor cells" %in% cell_types) next
  
  malignant_cells <- subset(seu, subset = cell_type == "Tumor cells")
  
  # Extract PCA embeddings (first 20 components)
  pca_mat <- Embeddings(malignant_cells, reduction = "pca")[, 1:20]
  
  # Compute Euclidean distance matrix
  dist_mat <- as.matrix(dist(pca_mat, method = "euclidean"))
  
  # For each cell, calculate the average distance to all other malignant cells
  # Set diagonal to NA to exclude self-distance (0)
  diag(dist_mat) <- NA
  iths <- rowMeans(dist_mat, na.rm = TRUE)
  
  #  Add ITH score to metadata
  malignant_cells$ITH <- iths


  # draw box plot
  df<- malignant_cells@meta.data
  df$condition<- factor(df$condition, levels = c('sensitive' , 'resistant'))
  
  # Calculate Tukey-based limits
  q1 <- quantile(df$ITH, 0.25, na.rm = TRUE)
  q3 <- quantile(df$ITH, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr + 5
  
  # pval
  wilcox_res <- wilcox.test(ITH ~ condition, data = df)
  
  p<- ggplot(df, aes(x = condition, y = ITH, fill = condition)) +
    geom_boxplot(width = 0.4, color = "black", outlier.shape = NA) +  # removes outliers
    stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +   # adds whisker caps
    labs(title = "ITH score", y = "ITH score", x= NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14, face = 'bold', color = 'black'), 
      axis.text.y = element_text(size = 14, face = 'bold', color = 'black'), 
      axis.title.y = element_text(size = 14, face = 'bold', color = 'black'),
      plot.title = element_text(face='bold', size=14, hjust = 0.5),
      legend.position = 'none'
    ) +
    scale_fill_manual(values = rev(col1)) +
    stat_summary(fun = median, geom = "point", shape = 20, size = 2.5, color = "black") +
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
  pdf_path <- file.path(out_dir, paste0(file_name, "_ITH_box.pdf"))
  png_path <- file.path(out_dir, paste0(file_name, "_ITH_box.png"))
  
  # Save PDF
  pdf(pdf_path, width = 4, height = 2.75)
  print(p)
  dev.off()
  
  # Save PNG
  png(png_path, width = 4, height = 2.75, units = "in", res = 300)
  print(p)
  dev.off()
  
  message(paste0(file_name, " - ITH saved."))
}
  


#'-------------------------------------------------------------------------