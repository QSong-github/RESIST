library(GSVA)
library(GSEABase)
library(tidyverse)
library(Seurat)
source('/blue/qsong1/sen.guo/resist_02/resist_test_share/scripts/color.R')

gmt_human <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/ref_data/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt"
gene_sets_human <- getGmt(gmt_human)
gene_sets_human<- geneIds(gene_sets_human)



# these are the datasets that have pre and post samples
ids_timppint <- c(
  "GSE111014_seurat_afterAnno.RDS",
  "GSE161195_seurat_afterAnno.RDS",
  "GSE161801_IMiD_seurat_afterAnno.RDS",
  "GSE161801_PI_seurat_afterAnno.RDS",
  "GSE162117_seurat_afterAnno.RDS",
  "GSE199333_seurat_afterAnno.RDS"
)


#'-------------------------------------------------------------------------
#' part 1: only pre sensitive and post resistant
#'-------------------------------------------------------------------------

out_dir <- '/blue/qsong1/sen.guo/resist_02/resist_test_share/results'
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
#rds_files <- list.files('../DRM_data', full.names = TRUE)
rds_files <- list.files('./data', pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)

# Loop through each file
for (file_path in rds_files) { ## remove the [10:33] to run all files
  file_name <- basename(file_path)
  dataset_id <- sub("_.*", "", file_name)
  
  if (file_name %in% ids_timppint) {
    next
  }
  
  seu<- readRDS(file_path)
  celltypes <- unique(seu$celltype)
  if (!"Malignant cells" %in% celltypes) next
  
  malignant_cells <- subset(seu, subset = celltype == "Malignant cells")
  
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


