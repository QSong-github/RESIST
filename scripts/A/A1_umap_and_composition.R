#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module A - Characterization of drug-resistant tumour ecosystems
# A1_umap_and_composition.R
# -----------------------------------------------------------------------------
# Purpose  : Plot existing UMAP embeddings by cluster, condition and cell type,
#            plus grouped cluster-composition bars within each condition.
# Inputs   : paths.data/*.rds or *.RDS; annotated Seurat objects with UMAP
# Outputs  : <results>/A/<input_filename>_UMAP.{pdf,png} (one composite figure)
# Usage    : Rscript scripts/A/A1_umap_and_composition.R
# Origin   : scripts/pipeline1_UMAP_replicate.R
# Notes    : Datasets with paired pre-/on-treatment sampling are skipped; see
#            `datasets.paired_pre_post_seu` in config/config.yaml.
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
outdir <- RESULTS_A
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
rds_files <- list.files(PATH_DATA, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)


#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

# Loop through each file
for (file_path in rds_files) {
  col <- resist_palette("celltype")   # restore the full palette for this dataset
  file_name <- basename(file_path)
  dataset_id <- sub("_.*", "", file_name)
  
  if (file_name %in% ids_timppint) {
    next
  }
  
  seu <- resist_read_seurat(file_path)
  cell_types <- unique(seu$cell_type)
  
  if ("Tumor cells" %in% cell_types) {
    seu$cell_type <- factor(seu$cell_type, levels = c("Tumor cells", setdiff(cell_types, "Tumor cells")))
  } else {
    seu$cell_type <- factor(seu$cell_type)  # keep existing order
    col <- col[-1]  # remove the first color if it's reserved for "Tumor cells"
  }
  
  # UMAP
  p1 <- DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'seurat_clusters', label.size = 4) + ggtitle("UMAP of cell clusters")+
    theme(plot.title = element_text(size = 15, face = "bold"),
          legend.position = 'none')
  p2 <- DimPlot(seu, reduction = "umap", label = FALSE, group.by = 'condition', cols = col1) + ggtitle("UMAP of conditions")+
    theme(plot.title = element_text(size = 15, face = "bold"),
          legend.position = 'none')
  p3 <- DimPlot(seu, reduction = "umap", label = FALSE, group.by = 'cell_type', cols = col) + ggtitle("UMAP of cell types")+
      theme(plot.title = element_text(size = 15, face = "bold"),
            legend.text = element_text(size = 9, color = 'black', face = 'bold'))
  
  # Bar chart data
  seu$condition <- as.factor(seu$condition)
  seu$seurat_clusters <- as.factor(seu$seurat_clusters)
  df_percent <- seu@meta.data %>%
    count(condition, seurat_clusters) %>%
    group_by(condition) %>%
    mutate(percent = n / sum(n) * 100)
  
  # Bar chart
  p_bar <- ggplot(df_percent, aes(x = seurat_clusters, y = percent, fill = condition)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Drug response\ncell composition",
         x = "Cluster",
         y = "Cell percentage") +
    theme_classic() +
    scale_fill_manual(values = col1) +
    theme(axis.text.x = element_text(angle = 45, size = 7, vjust = 1, hjust = 1,face = "bold", color = "black" ),
          axis.title.x = element_text( size = 8,  color = "black" ),
          axis.title.y = element_text(size = 8, face = "bold", color = "black" ),
          plot.title = element_text(size = 15, face = "bold", color = "black",hjust = 0.5),
          plot.margin = margin(0, 0, 0, 0),
          legend.position = 'none')
  
  # Save as PDF
  pdf_path <- file.path(outdir, paste0(file_name, "_UMAP.pdf"))
  pdf(pdf_path, width = 8, height = 5.5)
  print((p1 | p2) / (p3 | p_bar))  # 2x2 grid using patchwork
  dev.off()
  
  # Save as PNG
  png_path <- file.path(outdir, paste0(file_name, "_UMAP.png"))
  png(png_path, width = 8, height = 5.5, units = "in", res = 300)
  print((p1 | p2) / (p3 | p_bar))
  dev.off()
  
  print(paste0(file_name, ' is done!'))
}

#'-------------------------------------------------------------------------