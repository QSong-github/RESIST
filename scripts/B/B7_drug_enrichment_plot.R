#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module B - Transcriptional programmes of resistance
# B7_drug_enrichment_plot.R
# -----------------------------------------------------------------------------
# Purpose  : Render the connectivity scores of B6 in the layout used by the RESIST
#            browser.
# Inputs   : data/results/B/*_drug_enrichment.csv (from B6)
# Outputs  : data/results/B/*_drug_enrichment_plot.{pdf,png}
# Usage    : Rscript scripts/B/B7_drug_enrichment_plot.R
# Origin   : scripts/pipeline9_drug_enrichment_plot_for_web_replicate.R
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

#'---------------
#'  set up 
#' --------------


library(tidyverse)
color<- c('#BE9F95', '#EA9A9C', '#5FBBD0', '#DDDC94', '#A6D192',
'#56A456', '#BFB1D0', '#AED9E4', '#8C68AA', '#EABACD' )
out_dir<- RESULTS_B



#'---------------
#'  plot 
#' --------------
  # List all CSV files
files <- list.files(RESULTS_B, pattern = "_drug_enrichment\\.csv$",  full.names = TRUE)

for (file in files) {
  file_name <- basename(file)
  df <- read.csv(file)
  
  # Skip if there are no rows
  if (nrow(df) == 0) next
  
  # Subset data
  plot_df <- if (nrow(df) >= 10) df[1:10, ] else df
  
  # Factor to preserve order
  plot_df$pert_id <- factor(plot_df$pert_id, levels = plot_df$pert_id)
  
  # Generate plot
  p <- ggplot(plot_df, aes(x = pert_id, y = es_final, fill = pert_id)) +
  geom_col() +
  scale_fill_manual(values = color) +
  labs(x = "", y = "Enrichment Score", title = "Top 10 Enriched Drugs") +
  coord_cartesian(ylim = c(min(plot_df$es_final), max(plot_df$es_final))) + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black", size = 12),
  axis.title.x = element_text(size = 12, face = "bold", color = "black"),
  axis.text.y = element_text(size = 12, face = "bold", color = "black"),
  axis.title.y = element_text(size = 12, face = "bold", color = "black"),
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
  legend.position = 'none')
  
  ## single figure -------------------
  pdf_path <- file.path(out_dir, paste0(file_name, "_drug_enrichment_bar.pdf"))
  png_path <- file.path(out_dir, paste0(file_name, "_drug_enrichment_bar.png"))

  # Save PDF
  pdf(pdf_path, width = 4, height = 2.75)
  print(p)
  dev.off()

  # Save PNG
  png(png_path, width = 4, height = 2.75, units = "in", res = 300)
  print(p)
  dev.off()

  print(paste0(file_name, 'is done!'))
}