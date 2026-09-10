#'---------------
#'  set up 
#' --------------


library(tidyverse)
color<- c('#BE9F95', '#EA9A9C', '#5FBBD0', '#DDDC94', '#A6D192',
'#56A456', '#BFB1D0', '#AED9E4', '#8C68AA', '#EABACD' )
out_dir<- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'



#'---------------
#'  plot 
#' --------------
  # List all CSV files
files <- list.files("/blue/qsong1/sen.guo/resist_02/resist_test_share/results", pattern = "_drug_enrichment\\.csv$",  full.names = TRUE)

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