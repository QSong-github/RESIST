library(Seurat)
library(tidyverse)
library(patchwork)
source('/blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/color.R')

# these are the datasets that have pre and post samples
ids_timppint<- c('GSE120575_anti_PD1_seu.rds',
                 'GSE120575_anti_CTLA4_PD1_seu.rds',
                 'GSE123813_bcc_pembrolizumab_seu.rds',
                 'GSE123813_scc_pembrolizumab_seu.rds',
                 'GSE123813_scc_cemiplimab_seu.rds',
                 'GSE261898_seu.rds')



#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
} 
rds_files <- list.files('/blue/qsong1/sen.guo/resist_02/resist_test_share/data', pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)


# Loop through each file
for (file_path in rds_files) { # remove the [13] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_.*", "", file_name)
  
  if (file_name %in% ids_timppint) {
    next
  }
  
  seu<- readRDS(file_path)
  seu@meta.data$celltype[seu@meta.data$celltype == "Malignant cells"] <- "Tumor cells" ### the input rds data use "Malignant cells"
  colnames(seu@meta.data)[colnames(seu@meta.data) == "celltype"] <- "cell_type" ## change the column name to "cell_type" for consistency
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