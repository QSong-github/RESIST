#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C8_rbp_bubble_plot.R
# -----------------------------------------------------------------------------
# Purpose  : Summarise the RBP enrichment results of C7 as circle plots across
#            datasets and cell types.
# Inputs   : data/results/C/*_deg_rbp_enrichment.csv (from C7)
# Outputs  : data/results/C/*_deg_rbp_enrichment_{up,dn}_circle.{pdf,png}
# Usage    : Rscript scripts/C/C8_rbp_bubble_plot.R
# Origin   : scripts/pipeline8_RBP_plot1_bubble_replicate.R
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

library(forcats)
library(dplyr)
library(viridis)
library(packcircles)
library(tibble)
library(readr)
library(ggrepel)


rbp_enrich_files <- list.files(RESULTS_C, pattern = "_deg_rbp_enrichment\\.csv$", full.names = TRUE) # Get list of RBP files
out_dir <- RESULTS_C
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# these are mouse 
mouse_dat<- c('GSE201765_deg_rbp_enrichment.csv',
              'GSE213183_deg_rbp_enrichment.csv',
              'GSE261898_post_deg_rbp_enrichment.csv',
              'GSE261898_pre_deg_rbp_enrichment.csv')




#' #################### human ######################

# Loop through each DEG file
for (file_path in rbp_enrich_files ) {
  
  if (basename(file_path) %in% mouse_dat) next
  
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  # Load rbp_enrich data
  df <- read.csv(file_path)
  df_up<- df %>% filter(FDR_up<0.05)
  df_dn<- df %>% filter(FDR_dn<0.05)
  
  if (nrow(df_up) == 0 & nrow(df_dn) == 0 ) next
  if (nrow(df_up) <5 ) next
  
  df_up$log10FDR<- -log10(df_up$FDR_up)
  df_dn$log10FDR<- -log10(df_dn$FDR_dn)
  
  #'----------
  #' up 
  #' ---------
  packing <- circleProgressiveLayout(df_up$log10FDR, sizetype = "area") # Circle packing layout
  df_up <- cbind(df_up, packing)
  dat_circle <- circleLayoutVertices(packing, npoints = 50) # Generate circle polygon coordinates
  
  
  p_up<- ggplot() +
    geom_polygon(data = dat_circle, aes(x, y, group = id, fill = df_up$log10FDR[id]), color = "black", alpha = 0.8, size = 0.2) +
    geom_text(data = df_up, aes(x, y, label = RBP), size = 1.2, color = "white", fontface = "bold") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.text = element_text(hjust = 0.5, size = 5, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 5, face = "bold")
    )+
    coord_equal() +
    scale_fill_viridis(name = "-log10(Holm-adjusted p)",
                       option = "rocket",
                       direction = 1,
                       begin = 0.15,
                       end = 0.75,
                       guide = guide_colorbar(barheight = unit(1, "cm"),
                                              barwidth = 0.5)) +
    labs(title = "RBP Enrichment for Up-regulated DEGs", subtitle = "")
  
  path<- file.path(out_dir, paste0(dataset_id, '_deg_rbp_enrichment_up_circle.pdf'))
  pdf(path, width = 4, height = 2.75)
  print(p_up )
  graphics.off()
  
  path<- file.path(out_dir, paste0(dataset_id, '_deg_rbp_enrichment_up_circle.png'))
  png(path, width = 4, height = 2.75, units = "in", res = 300)
  print(p_up)
  dev.off()
  
  
  #'----------
  #' dn 
  #' ---------
  
  if (nrow(df_dn) <5 ) next
  
  packing <- circleProgressiveLayout(df_dn$log10FDR, sizetype = "radius") # Circle packing layout
  df_dn <- cbind(df_dn, packing)
  dat_circle <- circleLayoutVertices(packing, npoints = 50) # Generate circle polygon coordinates
  
  
  p_dn<- ggplot() +
    geom_polygon(data = dat_circle, aes(x, y, group = id, fill = df_dn$log10FDR[id]), color = "black", alpha = 0.8, size = 0.2) +
    geom_text(data = df_dn, aes(x, y, label = RBP), size = 2, color = "white", fontface = "bold") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.text = element_text(hjust = 0.5, size = 5, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 5, face = "bold")
    )+
    coord_equal() +
    scale_fill_viridis(name = "-log10(Holm-adjusted p)",
                       option = "viridis",
                       direction = -1,
                       guide = guide_colorbar(barheight = unit(1, "cm"),
                                              barwidth = 0.8)) +
    labs(title = "RBP Enrichment for Down-regulated DEGs", subtitle = "")
  
  
  path<- file.path(out_dir, paste0(dataset_id, '_deg_rbp_enrichment_dn_circle.pdf'))
  pdf(path, width = 4, height = 2.75)
  print(p_dn )
  graphics.off()
  
  path<- file.path(out_dir, paste0(dataset_id, '_deg_rbp_enrichment_dn_circle.png'))
  png(path, width = 4, height = 2.75, units = "in", res = 300)
  print(p_dn)
  dev.off()
  
}
  
  
  






