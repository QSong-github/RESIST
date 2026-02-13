###############################################################
# Single-cell Mutation Enrichment Analysis + AF Heatmap
#
# Upstream input (from cellSNP):
#   - cellSNP.cells.vcf.gz   (VCF with per-cell AD/DP in FORMAT fields)
#   - (other cellSNP outputs such as *.mtx, samples.tsv are not required here)
#
# Step 0 (VCF -> long-format AF table):
#   trans.py converts cellSNP.cells.vcf(.gz) into a per-variant-per-cell table:
#
#     Input : cellSNP.cells.vcf.gz
#     Output: cell_AF_table.tsv with columns:
#       - ID      : chrom_pos_ref_alt (e.g., chr22_43162920_A_G)
#       - Cell    : cell barcode (sample column name in VCF)
#       - s_reads : ALT-supporting reads (from AD)
#       - t_reads : total reads (from DP)
#       - AF      : allele frequency = s_reads / t_reads
#
#   Example:
#     python trans.py cellSNP.cells.vcf.gz cell_AF_table.tsv
#
# Step 1 (AF table -> Fisher test):
#   This R script merges mutation tables (one or multiple samples),
#   assigns cells into groups (sensitive vs resistant),
#   binarizes mutation status per cell using AF/read thresholds,
#   and performs Fisher's exact test per mutation ID.
#
#   Outputs:
#     - fisher_results.tsv  (ID, odds_ratio, p_value)
#
# Step 2 (Visualization):
#   Filters significant mutations and generates an AF heatmap across cells:
#     - heatmap.pdf
#     - heatmap.png
#
# Optional (Variant annotation):
#   format.py + mutation_annotation.py can be used to annotate significant loci
#   via dbSNP/ClinVar using Chr/Pos extracted from fisher_results.tsv.
###############################################################
# Load all required packages
###############################################################
#library(dplyr)
#library(readr)
#library(stringr)
#library(tidyverse)
#library(broom)

###############################################################
# Step 1: Load mutation files and standardize cell names
###############################################################

# Input directory containing mutation tables
data_dir <- "/home/liangjialu/orange_qsong1/liangjialu/GSE230538-result/gp1"

# Collect all .txt files
data_files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)

all_data <- data.frame()

for (file in data_files) {
  
  # Use filename as sample prefix
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_prefix <- gsub("[^a-zA-Z0-9]", "_", sample_name)
  
  df <- read_tsv(file, col_types = cols(.default = col_guess()))
  
  # Ensure required column exists
  if (!"Cell" %in% colnames(df)) {
    stop(paste("Missing 'Cell' column in:", file))
  }
  
  # Add sample prefix to cell barcode to ensure global uniqueness
  df$Cell <- paste0(sample_prefix, "_", df$Cell)
  
  all_data <- bind_rows(all_data, df)
}

###############################################################
# Step 2: Assign resistant/sensitive groups
###############################################################

all_data <- all_data %>%
  mutate(
    sample_prefix = str_extract(Cell, "^[^_]+"),
    Group = case_when(
      str_detect(sample_prefix, "resistant") ~ "resistant",
      str_detect(sample_prefix, "sensitive") ~ "sensitive",
      TRUE ~ NA_character_
    )
  )

infile <- all_data

###############################################################
# Step 3: Binarize mutation status per cell
###############################################################

# Convert AF to numeric
infile <- infile %>% mutate(AF = as.numeric(AF))

# Define mutation presence
# Criteria:
#   - Allele frequency >= 0.1
#   - Total supporting reads >= 10
#   - Variant-supporting reads >= 3
infile <- infile %>%
  mutate(
    mut = if_else(AF >= 0.1 & s_reads + t_reads >= 10 & s_reads >= 3, 1, 0)
  )

###############################################################
# Step 4: Perform Fisher's Exact Test per mutation
###############################################################

fisher_res <- infile %>%
  group_by(ID) %>%
  nest() %>%
  mutate(
    test = map(data, ~ {
      
      tbl <- table(.x$Group, .x$mut)
      
      # Skip invalid contingency tables
      if(ncol(tbl) < 2 | nrow(tbl) < 2){
        tibble(odds_ratio = NA, p_value = NA)
      } else {
        ft <- fisher.test(tbl, alternative = "two.sided")
        
        tidy(ft) %>%
          select(estimate, p.value) %>%
          rename(
            odds_ratio = estimate,
            p_value = p.value
          )
      }
    })
  ) %>%
  select(ID, test) %>%
  unnest(test)

# Save Fisher test results
fisher_res %>%
  arrange(p_value) %>%
  write.table(
    file = str_c(data_dir,"/fisher_results.tsv"),
    sep = "\t",
    row.names = FALSE
  )

###############################################################
# Step 5: Select significant and biologically relevant mutations
###############################################################

selected_mutations <- fisher_res %>%
  filter(
    p_value < 0.05,
    !is.na(odds_ratio),
    odds_ratio != 0,
    is.finite(odds_ratio)
  ) %>%
  pull(ID)

qualified_ids <- infile %>%
  filter(ID %in% selected_mutations) %>%
  group_by(ID) %>%
  summarise(
    prop_mut = mean(mut == 1, na.rm = TRUE),
    prop_AF = mean(AF[mut == 1] > 0.2, na.rm = TRUE)
  ) %>%
  filter(prop_mut > 0.2, prop_AF > 0.2) %>%
  pull(ID)

df_filtered <- infile %>%
  filter(ID %in% qualified_ids)

###############################################################
# Step 6: Construct mutation × cell AF matrix
###############################################################

library(reshape2)
library(grid)

df_filtered$AF <- as.numeric(df_filtered$AF)

# Remove duplicated mutation-cell pairs
df_filtered_dedup <- df_filtered %>%
  group_by(ID, Cell) %>%
  slice(1) %>%
  ungroup()

# Create wide matrix
heatmap_data <- reshape2::dcast(
  df_filtered_dedup,
  ID ~ Cell,
  value.var = "AF",
  fill = 0
)

rownames(heatmap_data) <- heatmap_data$ID
heatmap_data <- heatmap_data[, -1]

###############################################################
# Step 7: Prepare column annotations
###############################################################

cell_group <- df_filtered %>%
  select(Cell, Group) %>%
  distinct() %>%
  arrange(Group)

heatmap_data <- heatmap_data[, cell_group$Cell]

col_annotation <- data.frame(Group = cell_group$Group)
rownames(col_annotation) <- cell_group$Cell

###############################################################
# Step 8: Generate clustered heatmap
###############################################################

ph <- pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  silent = TRUE
)

row_order <- ph$tree_row$order
heatmap_data_clustered <- heatmap_data[row_order, ]

group_colors <- c(
  resistant = "#E3C07B",
  sensitive = "#3B80B5"
)

pheatmap(
  as.matrix(heatmap_data_clustered),
  annotation_col = col_annotation,
  annotation_colors = list(Group = group_colors),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "pink", "red"))(100),
  fontsize_row = 3,
  fontsize_col = 3,
  width = 4,
  height = 2.75,
  filename = file.path(data_dir, "heatmap.pdf")
)

pheatmap(
  as.matrix(heatmap_data_clustered),
  annotation_col = col_annotation,
  annotation_colors = list(Group = group_colors),
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "pink", "red"))(100),
  fontsize_row = 3,
  fontsize_col = 3,
  width = 4,
  height = 2.75,
  filename = file.path(data_dir, "heatmap.png")
)
