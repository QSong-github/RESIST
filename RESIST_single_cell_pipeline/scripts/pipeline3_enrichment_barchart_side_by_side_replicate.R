library(ggplot2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(msigdbr)
source('/blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/color.R')




#'#########
#' GO
#'#########

go_enrich_up_files <- list.files("/blue/qsong1/sen.guo/resist_02/resist_test_share/results", pattern = "_GO_enrichment_up\\.csv$",  full.names = TRUE) # Get list of RBP files
#go_enrich_dn_files <- list.files("./pipeline3_enrichment", pattern  = "_GO_enrichment_dn\\.csv$",  full.names = TRUE) # Get list of RBP files


out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


for (file_path in go_enrich_up_files ) { # remove the [11:12] to run all files
  
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  up<- read.csv(file_path)
  up<- up %>% filter(p.adjust<0.05)
  
  n<- nrow(up)
  if (n<10) {
    df_up <- up[order(up$p.adjust), ][1:n, ] 
  } else {
    df_up <- up[order(up$p.adjust), ][1:10, ] 
  }
  
  df_up$Description<- str_to_sentence(df_up$Description)
  df_up$`-log10padj`<- -log10(df_up$p.adjust)
  df_up$Description<- factor(df_up$Description, levels = rev(df_up$Description))
  
  dn_path <- sub("_up\\.csv$", "_dn.csv", file_path)
  dn<- read.csv(dn_path)
  dn<- dn %>% filter(p.adjust<0.05)

  n<- nrow(dn)
  if (n<10) {
    df_dn <- dn[order(dn$p.adjust), ][1:n, ] 
  } else {
    df_dn <- dn[order(dn$p.adjust), ][1:10, ] 
  }
  
  df_dn$Description<- str_to_sentence(df_dn$Description)
  df_dn$`-log10padj`<- -log10(df_dn$p.adjust)
  df_dn$Description<- factor(df_dn$Description, levels = rev(df_dn$Description))
  
  p_up <- ggplot(df_up, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#87CEEB", high = "#CD5C5C",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Up-regulated GO BP "
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  p_dn <- ggplot(df_dn, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#87CEEB", high = "#CD5C5C",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Down-regulated GO BP "
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  # Output paths
  pdf_path <- file.path(out_dir, paste0(dataset_id, "_GO_barchart_sideByside.pdf"))
  png_path <- file.path(out_dir, paste0(dataset_id, "_GO_barchart_sideByside.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 2.75)
  print(p_up | p_dn)
  dev.off()
  
  # Save PNG
  png(png_path, width = 8, height = 2.75, units = "in", res = 300)
  print(p_up | p_dn)
  dev.off()
  
}
  
  

#'#########
#' KEGG
#'#########

KEGG_enrich_up_files <- list.files("/blue/qsong1/sen.guo/resist_02/resist_test_share/results", pattern = "_KEGG_enrichment_up\\.csv$",  full.names = TRUE) # Get list of RBP files
#KEGG_enrich_dn_files <- list.files("./pipeline3_enrichment", pattern  = "_KEGG_enrichment_dn\\.csv$",  full.names = TRUE) # Get list of RBP files


out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

for (file_path in KEGG_enrich_up_files ) { # remove the [11:12] to run all files
  
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  up<- read.csv(file_path)
  up<- up %>% filter(p.adjust<0.05)
  
  n<- nrow(up)
  if (n<10) {
    df_up <- up[order(up$p.adjust), ][1:n, ] 
  } else {
    df_up <- up[order(up$p.adjust), ][1:10, ] 
  }
  
  df_up$Description<- str_to_sentence(df_up$Description)
  df_up$`-log10padj`<- -log10(df_up$p.adjust)
  df_up$Description<- factor(df_up$Description, levels = rev(df_up$Description))
  
  dn_path <- sub("_up\\.csv$", "_dn.csv", file_path)
  dn<- read.csv(dn_path)

  n<- nrow(dn)
  if (n<10) {
    df_dn <- dn[order(dn$p.adjust), ][1:n, ] 
  } else {
    df_dn <- dn[order(dn$p.adjust), ][1:10, ] 
  }

  df_dn$Description<- str_to_sentence(df_dn$Description)
  df_dn$`-log10padj`<- -log10(df_dn$p.adjust)
  df_dn$Description<- factor(df_dn$Description, levels = rev(df_dn$Description))
  
  p_up <- ggplot(df_up, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#175D8B", high = "#F9FBBB",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Up-regulated KEGG"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  p_dn <- ggplot(df_dn, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#175D8B", high = "#F9FBBB",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Down-regulated KEGG"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  # Output paths
  pdf_path <- file.path(out_dir, paste0(dataset_id, "_KEGG_barchart_sideByside.pdf"))
  png_path <- file.path(out_dir, paste0(dataset_id, "_KEGG_barchart_sideByside.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 2.75)
  print(p_up | p_dn)
  dev.off()
  
  # Save PNG
  png(png_path, width = 8, height = 2.75, units = "in", res = 300)
  print(p_up | p_dn)
  dev.off()
  
}










#'#########
#' HALLMARK
#'#########

HALLMARK_enrich_up_files <- list.files("/blue/qsong1/sen.guo/resist_02/resist_test_share/results", pattern = "_Hallmark_enrichment_up\\.csv$",  full.names = TRUE) # Get list of RBP files
#HALLMARK_enrich_dn_files <- list.files("./pipeline3_enrichment", pattern  = "_HALLMARK_enrichment_dn\\.csv$",  full.names = TRUE) # Get list of RBP files


out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

for (file_path in HALLMARK_enrich_up_files   ) { # remove the [11:12] to run all files
  
  file_name <- basename(file_path)
  dataset_id <- sub("_deg.*", "", file_name)
  
  up<- read.csv(file_path)
  up<- up %>% filter(p.adjust<0.05)
  
  n<- nrow(up)
  if (n<10) {
    df_up <- up[order(up$p.adjust), ][1:n, ] 
  } else {
    df_up <- up[order(up$p.adjust), ][1:10, ] 
  }
  
  df_up$Description<- str_to_sentence(df_up$Description)
  df_up$`-log10padj`<- -log10(df_up$p.adjust)
  df_up$Description<- factor(df_up$Description, levels = rev(df_up$Description))
  
  dn_path <- sub("_up\\.csv$", "_dn.csv", file_path)
  dn<- read.csv(dn_path)
  dn<- dn %>% filter(p.adjust<0.05)
  
  n<- nrow(dn)
  if (n<10) {
    df_dn <- dn[order(dn$p.adjust), ][1:n, ] 
  } else {
    df_dn <- dn[order(dn$p.adjust), ][1:10, ] 
  }
  
  df_dn$Description<- str_to_sentence(df_dn$Description)
  df_dn$`-log10padj`<- -log10(df_dn$p.adjust)
  df_dn$Description<- factor(df_dn$Description, levels = rev(df_dn$Description))
  
  p_up <- ggplot(df_up, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#7B2990", high = "#EC9160",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Up-regulated HALLMARK"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  p_dn <- ggplot(df_dn, aes(x = `-log10padj`, y = Description, fill = Count)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.3) +  # Add dashed line
    scale_fill_gradient(low = "#7B2990", high = "#EC9160",
                        guide = guide_colorbar(
                          direction = "horizontal", 
                          barheight = 0.3,  # height of the color bar (in "npc" units for vertical bar)
                          barwidth = 5,
                          title.position = "left",
                          title.hjust = 0.5)) +  # optional gradient
    labs(
      x = expression(-log[10]("adjusted p-value")),
      y = NULL,
      title = "Down-regulated HALLMARK"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5, face = "bold", color = "black"),
      axis.title = element_text(size = 6, face = "bold",color = "black"),
      plot.title = element_text(hjust = 0.5, size = 6, face = "bold"),
      legend.text = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.title = element_text(hjust = 0.5, size = 3, face = "bold"),
      legend.position = "bottom",
      legend.justification = c(0, 0),  # Bottom-left corner
      legend.box.just = "left"  
    )
  
  # Output paths
  pdf_path <- file.path(out_dir, paste0(dataset_id, "_HALLMARK_barchart_sideByside.pdf"))
  png_path <- file.path(out_dir, paste0(dataset_id, "_HALLMARK_barchart_sideByside.png"))
  
  # Save PDF
  pdf(pdf_path, width = 8, height = 2.75)
  print(p_up | p_dn)
  dev.off()
  
  # Save PNG
  png(png_path, width = 8, height = 2.75, units = "in", res = 300)
  print(p_up | p_dn)
  dev.off()
  
}











#'-------------------------------------------------------------------------