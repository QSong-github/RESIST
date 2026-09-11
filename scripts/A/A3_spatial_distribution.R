#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module A - Characterization of drug-resistant tumour ecosystems
# A3_spatial_distribution.R
# -----------------------------------------------------------------------------
# Purpose  : Map annotated cell types back onto tissue coordinates and summarise how
#            the resistant compartment is spatially organised relative to the
#            surrounding microenvironment.
# Inputs   : data/spatial/*.rds - Seurat objects carrying spatial coordinates
# Outputs  : data/results/A/spatial_plots/*.{pdf,png}
# Usage    : Rscript scripts/A/A3_spatial_distribution.R
# Origin   : scripts/final_spatial_distribution_v2.R
# Notes    : Developed against GSE284989; verify slot and image names before
#            applying it to a new spatial platform.
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
library(ggplot2)
library(patchwork)
library(dplyr)


### custmized and work for GSE284989 data
outdir <- file.path(RESULTS_A, "spatial_plots")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
rds_files <- list.files(PATH_SPATIAL, 
                        pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
# Loop through each spatial dataset
for (file_path in rds_files) {
  col <- resist_palette("celltype")   # restore the full palette for this dataset
  file_name <- basename(file_path)
  
  seu <- readRDS(file_path)
  
  # Standardize metadata column names
  seu@meta.data$final_celltype[seu@meta.data$final_celltype == "Tumor Cells"] <- "Tumor cells"
  colnames(seu@meta.data)[colnames(seu@meta.data) == "final_celltype"] <- "cell_type"
  colnames(seu@meta.data)[colnames(seu@meta.data) == "group"] <- "condition"
  
  cell_types <- unique(seu$cell_type)
  
  # Set factor levels (Tumor cells first if present)
  if ("Tumor cells" %in% cell_types) {
    seu$cell_type <- factor(seu$cell_type, levels = c("Tumor cells", setdiff(cell_types, "Tumor cells")))
  } else {
    seu$cell_type <- factor(seu$cell_type)
  }
  
  # Create NAMED color vector for all cell types
  all_cell_types <- levels(seu$cell_type)
  my_colors <- setNames(rep(unname(resist_palette("celltype")), length.out = length(all_cell_types)), all_cell_types)
  
  # Split by condition
  seu_resistant <- subset(seu, condition == "resistant")
  seu_sensitive <- subset(seu, condition == "sensitive")
  
  # Function to select top 2 samples with most cell type diversity
  # select_top_samples <- function(seu_subset, top_n = 2) {
  #   diversity <- seu_subset@meta.data %>%
  #     group_by(sample) %>%
  #     summarise(n_cell_types = n_distinct(cell_type), .groups = "drop") %>%
  #     arrange(desc(n_cell_types)) %>%
  #     slice(1:min(top_n, n()))
    
  #   return(diversity$sample)
  # }
  select_top_samples <- function(seu_subset, top_n = 2) {
  diversity <- seu_subset@meta.data %>%
    group_by(sample) %>%
    summarise(n_cell_types = n_distinct(cell_type), .groups = "drop") %>%
    arrange(desc(n_cell_types)) %>%
    slice(1:min(top_n, n()))
  
  # Convert sample IDs to image names
  sample_ids <- diversity$sample
  image_names <- paste0("slice1.", sample_ids)
  
  # Get available images in the subset
  available_images <- Images(seu_subset)
  
  # Return only images that exist
  return(image_names[image_names %in% available_images])
  }

  # Get top 2 samples per condition (as character vector)
  top_resistant <- as.character(select_top_samples(seu_resistant, top_n = 2))
  top_sensitive <- as.character(select_top_samples(seu_sensitive, top_n = 2))
  
  # Skip if fewer than 2 samples
  # if (length(top_resistant) < 2 | length(top_sensitive) < 2) {
  #   print(paste0(file_name, ' skipped - insufficient samples'))
  #   next
  # }

  # Get top images per condition
  top_resistant_images <- select_top_samples(seu_resistant, top_n = 2)
  top_sensitive_images <- select_top_samples(seu_sensitive, top_n = 2)

  # Reverse order
  top_resistant_ordered <- rev(top_resistant_images)
  top_sensitive_ordered <- rev(top_sensitive_images)

  # Adjust ncol based on actual number of samples
  ncol_resistant <- length(top_resistant_ordered)
  ncol_sensitive <- length(top_sensitive_ordered)
  
  # Subset using sample names (extract from image names)
  resistant_samples <- gsub("slice1\\.", "", top_resistant_ordered)
  sensitive_samples <- gsub("slice1\\.", "", top_sensitive_ordered)
  
  # Subset to ONLY top 2 samples
  seu_resistant_top <- subset(seu_resistant, sample %in% resistant_samples)
  seu_sensitive_top <- subset(seu_sensitive, sample %in% sensitive_samples)
  ##

  # Create clean titles (remove "slice1." prefix)
  clean_resistant_titles <- gsub("slice1\\.", "", top_resistant_ordered)
  clean_sensitive_titles <- gsub("slice1\\.", "", top_sensitive_ordered)


  # RESISTANT plot: legend on the right, smaller sample titles
  p_resistant <- SpatialDimPlot(seu_resistant_top, 
                                 group.by = "cell_type", 
                                 images = top_resistant_ordered,
                                 cols = my_colors,
                                 pt.size.factor = 1.6,
                                 ncol = 2, ## change from 2 to ncol_resistant
                                 combine = TRUE) &
    theme(plot.title = element_text(size = 7.5),
          legend.position = "right",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          #strip.text.x = element_text(size = 6),  # Smaller sample titles
          legend.key.height = unit(0.3, "lines"),
          legend.spacing.y = unit(0.2, "lines"))
  # Replace titles
  for (i in seq_along(p_resistant)) {
    p_resistant[[i]] <- p_resistant[[i]] + ggtitle(clean_resistant_titles[i])
  }
  # SENSITIVE plot: legend on the right, smaller sample titles
  p_sensitive <- SpatialDimPlot(seu_sensitive_top, 
                                 group.by = "cell_type", 
                                 images = top_sensitive_ordered,
                                 cols = my_colors,
                                 pt.size.factor = 1.6,
                                 ncol = 2, ## change from 2 to ncol_sensitive
                                 combine = TRUE) &
    theme(plot.title = element_text(size = 7.5),
          legend.position = "right",
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 8),
          strip.text.x = element_text(size = 6),  # Smaller sample titles
          legend.key.height = unit(0.3, "lines"),
          legend.spacing.y = unit(0.2, "lines"))
  
  # Replace titles
  for (i in seq_along(p_sensitive)) {
    p_sensitive[[i]] <- p_sensitive[[i]] + ggtitle(clean_sensitive_titles[i])
  }

  # Add main title with LARGER font size
  p_resistant_final <- p_resistant +
    plot_annotation(title = "Spatial distribution of cell types in\nresistant samples",
                    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))
  
  p_sensitive_final <- p_sensitive +
    plot_annotation(title = "Spatial distribution of cell types in\nsensitive samples",
                    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))
  
  # Fixed dimensions for side-by-side layout
  plot_width <- 8
  plot_height <- 2.75
  
  # Save RESISTANT plot
  pdf_path_res <- file.path(outdir, paste0(file_name, "_spatial_cell_distribution_resistant_top2.pdf"))
  pdf(pdf_path_res, width = plot_width, height = plot_height)
  print(p_resistant_final)
  dev.off()
  
  png_path_res <- file.path(outdir, paste0(file_name, "_spatial_cell_distribution_resistant_top2.png"))
  png(png_path_res, width = plot_width, height = plot_height, units = "in", res = 300)
  print(p_resistant_final)
  dev.off()
  
  # Save SENSITIVE plot
  pdf_path_sen <- file.path(outdir, paste0(file_name, "_spatial_cell_distribution_sensitive_top2.pdf"))
  pdf(pdf_path_sen, width = plot_width, height = plot_height)
  print(p_sensitive_final)
  dev.off()
  
  png_path_sen <- file.path(outdir, paste0(file_name, "_spatial_cell_distribution_sensitive_top2.png"))
  png(png_path_sen, width = plot_width, height = plot_height, units = "in", res = 300)
  print(p_sensitive_final)
  dev.off()
  
  print(paste0(file_name, ' spatial plots (resistant & sensitive top 2) done!'))
}
