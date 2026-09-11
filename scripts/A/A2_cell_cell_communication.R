#!/usr/bin/env Rscript
# =============================================================================
# RESIST | Module A - Characterization of drug-resistant tumour ecosystems
# A2_cell_cell_communication.R
# -----------------------------------------------------------------------------
# Purpose  : Infer ligand-receptor mediated signalling between the cell types of the
#            tumour ecosystem separately in the sensitive and resistant states, and
#            contrast the two networks to expose interactions that are gained or lost
#            with resistance.
# Inputs   : data/spatial/*.rds - annotated Seurat objects
# Outputs  : data/results/A/<dataset>_cellchat*.{pdf,png}
# Usage    : Rscript scripts/A/A2_cell_cell_communication.R
# Origin   : scripts/pipeline10_cellchat_batch.R
# Notes    : Uses the vendored `lib/netVisual_bubble.R`, a patched CellChat
#            bubble plot that drops NA interactions.
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
library(CellChat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(stringr)
library(grid)
# Set stringsAsFactors to FALSE globally
options(stringsAsFactors = FALSE)
###
source(file.path(RESIST_HOME, "scripts", "lib", "netVisual_bubble.R")) ## modified function to remove NAs

# these are the datasets that have pre and post samples
ids_timppint <- resist_datasets("paired_pre_post_seu")
outdir <- RESULTS_A
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
rds_files <- list.files(PATH_DATA, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)




#'-------------------------------------------------------------------------
#' 
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
  #seu@meta.data$final_celltype[seu@meta.data$final_celltype == "Tumor Cells"] <- "Tumor cells" ### the input rds data use " Tumor Cells" to represent tumor cells, change it to "Tumor cells" for consistency
  #colnames(seu@meta.data)[colnames(seu@meta.data) == "cell_type"] <- "cell_type" ## change the column name to "cell_type" for consistency
  #colnames(seu@meta.data)[colnames(seu@meta.data) == "sctype_classification"] <- "cell_type" ## change the column name to "condition" for consistency
  #seu$cell_type <- seu$sctype_classification
  cell_types <- unique(seu$cell_type)
  
  #seu@meta.data$cell_type[seu@meta.data$cell_type == "Cancer cells"] <- "Tumor cells"

  if ("Tumor cells" %in% cell_types) {
    seu$cell_type <- factor(seu$cell_type, levels = c("Tumor cells", setdiff(cell_types, "Tumor cells")))
  } else {
    seu$cell_type <- factor(seu$cell_type)  # keep existing order
    col <- col[-1]  # remove the first color if it's reserved for "Tumor cells"
  }
  

  #


#---------------------------------------------------
### Cellchat analysis
#---------------------------------------------------

# Use non-integrated assay
# Split once by condition
obj.list <- SplitObject(seu, split.by = "condition")  # obj.list$sensitive, obj.list$resistant

# Build CellChat for a Seurat subset (no DefaultAssay/Idents changes needed)
mkCellChat_group <- function(seu_group,
                             assay = "SCT",
                             celltype_col = "cell_type",
                             sample_col = "sample",
                             min.cells = 10,
                             db = CellChatDB.human) {

  stopifnot(all(c(celltype_col, sample_col) %in% colnames(seu_group@meta.data)))
  cells <- colnames(seu_group)

  # Expression from specified assay
  expr <- GetAssayData(seu_group, assay = assay, slot = "data")[, cells, drop = FALSE]

  # Metadata with proper factors and a samples column
  ct  <- seu_group@meta.data[, celltype_col]
  smp <- seu_group@meta.data[, sample_col]
  meta <- data.frame(
    cell_type = droplevels(factor(ct)),
    samples   = droplevels(factor(smp)),
    row.names = cells
  )

  # Drop very rare cell types to avoid downstream errors
  keep_types <- names(which(table(meta$cell_type) >= min.cells))
  keep_cells <- rownames(meta)[meta$cell_type %in% keep_types]
  meta <- meta[keep_cells, , drop = FALSE]
  expr <- expr[, keep_cells, drop = FALSE]

  # CRITICAL: Re-drop factor levels after filtering !!!!!!!!!!!
  meta$cell_type <- droplevels(meta$cell_type)
  meta$samples   <- droplevels(meta$samples)
  
  if (nlevels(meta$cell_type) < 2) {
    stop(sprintf("Fewer than 2 cell types with >= %d cells after filtering.", min.cells))
  }

  # Build CellChat
  cc <- createCellChat(object = expr, meta = meta, group.by = "cell_type")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, population.size = TRUE)
  cc <- filterCommunication(cc, min.cells = min.cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

# Create CellChat objects
cellchat.sens <- mkCellChat_group(obj.list$sensitive, assay = "SCT", db = CellChatDB.human)
cellchat.res  <- mkCellChat_group(obj.list$resistant,  assay = "SCT", db = CellChatDB.human)

#
#save.image(file.path(outdir, paste0(file_name, "_cellchat_objects.RData"))) ## save RData for future use


#-------------------------------------------------
# Visualize the  interactions for each condition
plot_circle <- function(cc, title) {
  group.size <- as.numeric(table(cc@idents)[levels(cc@idents)])
  netVisual_circle(cc@net$count,
                   vertex.weight = group.size,
                   weight.scale = TRUE,
                   label.edge = FALSE,
                   title.name = title)
}




###
# Save Sensitive (right) and Resistant (left) circle plots in one PDF/PNG
save_circle_pair <- function(cc_left, cc_right,
                             titles = c("Resistant", "Sensitive"),
                             pdf_file, png_file,
                             width = 8, height = 2.75) {
  get_w <- function(cc) as.numeric(table(cc@idents)[rownames(cc@net$count)])

  wL <- get_w(cc_left)
  wR <- get_w(cc_right)

  # PDF
  pdf(pdf_file, width = width, height = height)
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5))  # Increased top margin for titles
  
  netVisual_circle(cc_left@net$count, 
                   vertex.weight = wL, 
                   weight.scale = TRUE,vertex.label.cex = 0.8,
                   label.edge = FALSE)
  title(main = titles[1], line = 1, cex.main = 1.2)  # Add title explicitly
  
  netVisual_circle(cc_right@net$count, 
                   vertex.weight = wR, 
                   weight.scale = TRUE,vertex.label.cex = 0.8,
                   label.edge = FALSE)
  title(main = titles[2], line = 1, cex.main = 1.2)  # Add title explicitly
  
  dev.off()

  # PNG
  png(png_file, width = width, height = height, units = "in", res = 300)
  par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5))  # Increased top margin for titles
  
  netVisual_circle(cc_left@net$count, 
                   vertex.weight = wL, 
                   weight.scale = TRUE,vertex.label.cex = 0.8,
                   label.edge = FALSE)
  title(main = titles[1], line = 1, cex.main = 1.2)  # Add title explicitly
  
  netVisual_circle(cc_right@net$count, 
                   vertex.weight = wR, 
                   weight.scale = TRUE,vertex.label.cex = 0.8 ,
                   label.edge = FALSE)
  title(main = titles[2], line = 1, cex.main = 1.2)  # Add title explicitly
  
  dev.off()
}

# Use it to save circle plots
save_circle_pair(cellchat.res, cellchat.sens,
                 pdf_file = file.path(outdir, paste0(file_name, "circle_sensitive_resistant.pdf")),
                 png_file = file.path(outdir, paste0(file_name, "circle_sensitive_resistant.png")))

#-------------------------------------------------
# Merge and compare two CellChat objects
#-------------------------------------------------
# To ensure fair comparison, only keep cell types present in both conditions
###
#load(file.path(outdir, paste0(file_name, "_cellchat_objects.RData")) ) ## load RData if needed
cellchat.sens2 <- cellchat.sens
cellchat.res2 <- cellchat.res

# Find cell types present in BOTH conditions
shared_types <- intersect(levels(cellchat.sens2@idents), 
                          levels(cellchat.res2@idents))
# Keep only those with cells in both
shared_types <- shared_types[table(cellchat.sens2@idents)[shared_types] > 0 & 
                              table(cellchat.res2@idents)[shared_types] > 0]

# Subset both objects to shared types
cellchat.sens2 <- subsetCellChat(cellchat.sens2, idents.use = shared_types)
cellchat.res2 <- subsetCellChat(cellchat.res2, idents.use = shared_types)

# Then merge
cellchat.merged <- mergeCellChat(
  #list(Sensitive = cellchat.sens2, Resistant = cellchat.res2),
  #add.names = c("Sensitive","Resistant")
    list(Resistant = cellchat.res2, Sensitive = cellchat.sens2),
    add.names = c("Resistant","Sensitive")

)
#####


  #------------------------------------------------------
  # Identify differentially weighted signaling pathways
  #-------------------------------------------------------

  # Visualize the overall differential interactions
  # Save combined diffInteraction plots (Overall vs Tumor cells as source)

# 
save_diffInteraction_combined_advanced <- function(cellchat.obj,
                                                   pdf_file, png_file,
                                                   width = 8, height = 2.75,
                                                   title_left = "Overall Interactions",
                                                   title_right = "Tumor Cells Signaling",
                                                   title.size = 1,
                                                   mar_custom = c(1, 1, 2.5, 1)) {
  
  # PDF
  pdf(pdf_file, width = width, height = height)
  par(mfrow = c(1, 2), 
      mar = mar_custom,
      xpd = TRUE)
  
  netVisual_diffInteraction(cellchat.obj, weight.scale = TRUE, title.name = "")
  title(main = title_left, line = 1.2, cex.main = title.size, font.main = 1.5)
  
  netVisual_diffInteraction(cellchat.obj, weight.scale = TRUE, sources.use = "Tumor cells", title.name = "")
  title(main = title_right, line = 1.2, cex.main = title.size, font.main = 1.5)
  
  dev.off()
  
  # PNG
  png(png_file, width = width, height = height, units = "in", res = 300)
  par(mfrow = c(1, 2), 
      mar = mar_custom,
      xpd = TRUE)
  
  netVisual_diffInteraction(cellchat.obj, weight.scale = TRUE, title.name = "")
  title(main = title_left, line = 1.2, cex.main = title.size, font.main = 1.5)
  
  netVisual_diffInteraction(cellchat.obj, weight.scale = TRUE, sources.use = "Tumor cells", title.name = "")
  title(main = title_right, line = 1.2, cex.main = title.size, font.main = 1.5)
  
  dev.off()
  
  cat("Saved:", pdf_file, "\n")
  cat("Saved:", png_file, "\n")
}

# Example usage with custom titles and sizing
save_diffInteraction_combined_advanced(
  cellchat.obj = cellchat.merged,
  pdf_file = file.path(outdir, paste0(file_name, "_cellchat_compare_interactions_combined.pdf")),
  png_file = file.path(outdir, paste0(file_name, "_cellchat_compare_interactions_combined.png")),
  width = 8,
  height = 2.75,
  title_left = "All Cell-Cell Interactions (Resistant vs Sensitive)",
  title_right = "Tumor Cell Signaling (Resistant vs Sensitive)",
  title.size = 1,
  mar_custom = c(1, 1, 2.5, 1)
)


 #----------------------------------------------------------------------------------
    ## Identify and visualize the top differentially weighted ligand–receptor pairs
  #----------------------------------------------------------------------------------
  # 1) Build Top 10 ligand–receptor pairs by absolute difference (Resistant - Sensitive)
  comm_res  <- subsetCommunication(cellchat.res2)
  comm_sens <- subsetCommunication(cellchat.sens2)
  key_cols <- c("source","target","interaction_name")
  m <- merge(comm_res[, c(key_cols, "prob")],
        comm_sens[, c(key_cols, "prob")],
        by = key_cols, all = TRUE, suffixes = c(".res", ".sens"))
  m$prob.res[is.na(m$prob.res)]   <- 0
  m$prob.sens[is.na(m$prob.sens)] <- 0
  m$diff <- m$prob.res - m$prob.sens  # Resistant minus Sensitive

  # OPTIONAL: focus on a sender/recipient subset (e.g., Tumor as source)
  sources <- "Tumor cells"
  targets <- NULL
  if (!is.null(targets)) {
      m <- m[m$source %in% sources & m$target %in% targets, ]
  } else {
      m <- m[m$source %in% sources, ]
  }


# Top 10 LR by absolute difference across the selected edges
top10_lr <- head(m[order(-abs(m$diff)), "interaction_name"], 10)
pairLR.use <- data.frame(interaction_name = unique(top10_lr), stringsAsFactors = FALSE)



# 2) Visualize the differences for these top 10 LR pairs
# Save as PDF
pdf_path <- file.path(outdir, paste0(file_name, "_cellchat_diff_LR_pair_top.pdf"))
pdf(pdf_path, width = 4, height = 2.75)

p <- netVisual_bubble(
    object = cellchat.merged,        
    pairLR.use = pairLR.use,         
    comparison = c(1, 2),            
    sources.use = sources,           
    targets.use = targets,           
    remove.isolate = TRUE,
    font.size = 8, 
    angle.x = 45, 
    color.text = c('#E4C66F', '#5B9BD5'),
    title.name = "Top L-R: Resistant vs. Sensitive",show.legend = T,
    return.data = FALSE
)
# Extract plot and data


# Modify the plot to remove condition labels from x-axis only
if(inherits(p, "ggplot")) {
    # Extract and modify x-axis labels
    p_build <- ggplot_build(p)
    x_labels <- p_build$layout$panel_params[[1]]$x$get_labels()
    x_labels_clean <- gsub(" \\(Resistant\\)| \\(Sensitive\\)", "", x_labels)
    
    p <- p + 
        scale_x_discrete(labels = x_labels_clean) +theme(
        legend.key.size = unit(0.4, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0.2, "cm"),        # Increased spacing to show all items
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.15, "cm")      # Space between legend boxes
    )
    
    
    print(p)
    # Helper: draw a custom condition legend at bottom-right of the current plot device
    add_condition_legend_right <- function(
    x_right_mm = 19.5,      # distance from the right edge
    y_bottom_mm = 16,     # distance up from the bottom edge
    c_res = "#E4C66F", c_sens = "#5B9BD5",
    size_mm = 3, text_size = 8, spacing_mm = 5
    ) {
    x <- grid::unit(1, "npc") - grid::unit(x_right_mm, "mm")
    y <- grid::unit(0, "npc") + grid::unit(y_bottom_mm, "mm")

    # Resistant
    grid::grid.points(x = x, y = y,
                        pch = 18, size = grid::unit(size_mm, "mm"),
                        gp = grid::gpar(col = c_res))
    grid::grid.text("Resistant",
                    x = x + grid::unit(size_mm + 2, "mm"),
                    y = y,
                    just = "left",
                    gp = grid::gpar(col = c_res, fontsize = text_size))

    # Sensitive (stacked below)
    grid::grid.points(x = x, y = y - grid::unit(spacing_mm, "mm"),
                        pch = 16, size = grid::unit(size_mm, "mm"),
                        gp = grid::gpar(col = c_sens))
    grid::grid.text("Sensitive",
                    x = x + grid::unit(size_mm + 2, "mm"),
                    y = y - grid::unit(spacing_mm, "mm"),
                    just = "left",
                    gp = grid::gpar(col = c_sens, fontsize = text_size))
    }
    add_condition_legend_right(x_right_mm = 19.5, y_bottom_mm = 16,
                            c_res = '#E4C66F', c_sens = '#5B9BD5')
} else {
    print(p)
}

dev.off()

# Save as PNG
png_path <- file.path(outdir, paste0(file_name, "_cellchat_diff_LR_pair_top.png"))
png(png_path, width = 4, height = 2.75, units = "in", res = 300)

p <- netVisual_bubble(
    object = cellchat.merged,        
    pairLR.use = pairLR.use,         
    comparison = c(1, 2),            
    sources.use = sources,           
    targets.use = targets,           
    remove.isolate = TRUE,
    font.size = 10, 
    angle.x = 45, 
    color.text = c('#E4C66F', '#5B9BD5'),
    title.name = "Top L-R: Resistant vs. Sensitive",show.legend = T,
    return.data = FALSE
)

# Modify the plot
if(inherits(p, "ggplot")) {
    # Extract and modify x-axis labels
    p_build <- ggplot_build(p)
    x_labels <- p_build$layout$panel_params[[1]]$x$get_labels()
    x_labels_clean <- gsub(" \\(Resistant\\)| \\(Sensitive\\)", "", x_labels)
    
    p <- p + 
        scale_x_discrete(labels = x_labels_clean) + theme(
        legend.key.size = unit(0.4, "cm"),
        legend.key.height = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.spacing.y = unit(0.2, "cm"),        # Increased spacing to show all items
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.15, "cm")      # Space between legend boxes
    )
    print(p)
    add_condition_legend_right(x_right_mm = 19.5, y_bottom_mm = 16,
                            c_res = '#E4C66F', c_sens = '#5B9BD5')
    # Add custom legend with colored dots
   
} else {
    print(p)
}

dev.off()

print(paste0(file_name, ' is done!'))
}
