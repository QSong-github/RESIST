library(ggplot2)
library(dplyr)
library(stringr)
library(dplyr)
library(forcats)
library(treemap)
library(RColorBrewer)
library(ggsci)
library(treemapify)


colors<- c("#8DD3C7" ,"#FFD700", "#BEBADA", "#FB8072" ,
           "#80B1D3", "#FDB462" ,"#FCCDE5",
           "#D9D9D9" ,"#CCEBC5", "#FFED6F"
)


out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
} 
in_dir <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/results"

up_files <- list.files(in_dir, pattern = "up_miRNA\\.csv$", full.names = TRUE)
dn_files <- list.files(in_dir, pattern = "dn_miRNA\\.csv$", full.names = TRUE)


#'----------------
#' human
#' ---------------
##-----------------------up
for (file_path in up_files) { #remove the [18] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  df_up<- read.csv(file_path)
  df_up <- df_up %>% filter(Score>=90)
  
  if (length(unique(df_up$hgnc_symbol)) <= 5) next
  
  df_up<- df_up %>%
    group_by(hgnc_symbol) %>%
    arrange(desc(Score)) %>%
    slice_head(n = 4)
  
  df_up$hgnc_symbol <- as.character(df_up$hgnc_symbol)
  df_up$miRNA <- as.character(df_up$miRNA)
  df_up$miRNA_label <- str_wrap(df_up$miRNA, width = 10)
  
  
  df_up$Size <- 1 
  
  
  p<- ggplot(df_up, aes(
    area = Size,
    fill = hgnc_symbol,
    label = miRNA,
    subgroup = hgnc_symbol
  )) +
    geom_treemap(color = "black", size = 0.3) +
    geom_treemap_subgroup_border(color = "white", size = 0.8) +  # thicker border for gene groups
    geom_treemap_text(
      fontface = "plain", 
      colour = "black", 
      place = "centre", 
      grow = FALSE,
      reflow = T,
      size = 14
    ) +
    geom_treemap_subgroup_text(
      place = "topleft",
      grow = F,
      alpha = 0.8,
      colour = "black",
      fontface = "bold",
      size = 14,
      min.size = 8
    ) +
    scale_fill_manual(values = colors) +
    labs(title = "miRNA for Up-regulated DEGs") +
    theme(legend.position = "none") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  
  pdf_path <- file.path(out_dir, paste0(dataset_id, "_miRNA_up.pdf"))
  png_path <- file.path(out_dir, paste0(dataset_id, "_miRNA_up.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 5.5)
  print(p)
  graphics.off()
  
  # Save PNG
  png(png_path, width = 8, height = 5.5, units = "in", res = 300)
  print(p)
  dev.off()
  
  message(paste0(file_name, " - miRNA saved."))
}
  
  




##-----------------------dn
for (file_path in dn_files) { #remove the [18] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  df_dn<- read.csv(file_path)
  df_dn <- df_dn %>% filter(Score>=90)
  
  if (length(unique(df_dn$hgnc_symbol)) <= 5) next
  
  df_dn<- df_dn %>%
    group_by(hgnc_symbol) %>%
    arrange(desc(Score)) %>%
    slice_head(n = 4)
  
  df_dn$hgnc_symbol <- as.character(df_dn$hgnc_symbol)
  df_dn$miRNA <- as.character(df_dn$miRNA)
  df_dn$miRNA_label <- str_wrap(df_dn$miRNA, width = 10)
  
  
  df_dn$Size <- 1 
  
  
  p<- ggplot(df_dn, aes(
    area = Size,
    fill = hgnc_symbol,
    label = miRNA,
    subgroup = hgnc_symbol
  )) +
    geom_treemap(color = "black", size = 0.3) +
    geom_treemap_subgroup_border(color = "white", size = 0.8) +  # thicker border for gene grodns
    geom_treemap_text(
      fontface = "plain", 
      colour = "black", 
      place = "centre", 
      grow = FALSE,
      reflow = T,
      size = 14
    ) +
    geom_treemap_subgroup_text(
      place = "topleft",
      grow = F,
      alpha = 0.8,
      colour = "black",
      fontface = "bold",
      size = 14,
      min.size = 8
    ) +
    scale_fill_manual(values = colors) +
    labs(title = "miRNA for Down-regulated DEGs") +
    theme(legend.position = "none") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
  
  
  pdf_path <- file.path(out_dir, paste0(dataset_id, "_miRNA_dn.pdf"))
  png_path <- file.path(out_dir, paste0(dataset_id, "_miRNA_dn.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 5.5)
  print(p)
  graphics.off()
  
  # Save PNG
  png(png_path, width = 8, height = 5.5, units = "in", res = 300)
  print(p)
  dev.off()
  
  message(paste0(file_name, " - miRNA saved."))
}


