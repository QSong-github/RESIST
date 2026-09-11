#!/usr/bin/env Rscript
# Read-only prerequisite checks. This script does not execute analyses.
home <- Sys.getenv("RESIST_HOME")
if (!nzchar(home)) stop("Run through run_A.sh, run_B.sh, run_C.sh, or run_D.sh.")
steps <- strsplit(Sys.getenv("RESIST_STEPS"), ",", fixed = TRUE)[[1]]
required <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
required(requireNamespace("yaml", quietly = TRUE), "Install R package yaml; see TUTORIAL.md, Step 1.")
required(requireNamespace("jsonlite", quietly = TRUE), "Install R package jsonlite; see TUTORIAL.md, Step 1.")
source(file.path(home, "scripts/lib/config.R"))
source(file.path(home, "scripts/lib/input.R"))
CFG <- resist_load_config(home)
paths <- lapply(CFG$paths, resist_resolve, home = home)
PATH_REF <- paths$ref_data
PATH_RESULTS <- paths$results
spec <- jsonlite::read_json(file.path(home, "config/steps.json"), simplifyVector = FALSE)
packages <- unique(c(unlist(spec$shared_packages), unlist(lapply(steps, function(s) spec$steps[[s]]$packages))))
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
required(!length(missing), paste("Missing R packages:", paste(missing, collapse = ", "),
                                "\nInstall on the HPC using config/setup/install_r_packages.R."))
suppressPackageStartupMessages(library(Seurat))
cat("Configuration:", normalizePath(Sys.getenv("RESIST_CONFIG")), "\n")
cat("Input:", paths$data, "\nReference:", PATH_REF, "\nOutput:", PATH_RESULTS, "\n")
for (s in steps) {
  refs <- unlist(spec$steps[[s]]$references)
  for (key in refs) {
    value <- CFG$reference_files[[key]]
    required(!is.null(value), paste("Missing reference_files key:", key))
    path <- file.path(PATH_REF, value)
    required(file.exists(path), paste("Missing reference:", path,
      "\nOn the HPC, unpack data/reference_bundle.tar.gz or obtain the external file; see docs/references.md."))
  }
}

primary_steps <- intersect(steps, c("A1", "A2", "B1", "B4", "B5"))
if (length(primary_steps)) {
  files <- list.files(paths$data, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  required(length(files) > 0L, paste("No Seurat .rds/.RDS files in", paths$data,
    "\nDownload the example on the HPC; see TUTORIAL.md, Step 2."))
  for (s in primary_steps) {
    exclusions <- unlist(CFG$datasets[[if (s == "B5") "paired_pre_post_anno" else "paired_pre_post_seu"]])
    if (s == "B1") exclusions <- c(exclusions, unlist(CFG$datasets$harmony_integrated))
    required(any(!basename(files) %in% exclusions), paste(s, "has no eligible input after configured exclusions."))
  }
  for (f in files) {
    active <- Filter(function(s) {
      exc <- unlist(CFG$datasets[[if (s == "B5") "paired_pre_post_anno" else "paired_pre_post_seu"]])
      if (s == "B1") exc <- c(exc, unlist(CFG$datasets$harmony_integrated))
      !basename(f) %in% exc
    }, primary_steps)
    if (!length(active)) next
    x <- resist_read_seurat(f)
    cat("\nObject:", basename(f), "| cells:", ncol(x), "| active assay:", DefaultAssay(x), "\n")
    print(table(x$cell_type, x$condition))
    cat("Input MD5:", unname(tools::md5sum(f)), "\n")
    if ("A1" %in% active) {
      required("umap" %in% Reductions(x), paste("A1 requires an existing umap reduction:", f))
      required("seurat_clusters" %in% names(x[[]]), paste("A1 requires seurat_clusters:", f))
    }
    if ("B1" %in% active) {
      ct <- table(x$cell_type, factor(x$condition, levels = c("sensitive", "resistant")))
      required(any(rowSums(ct >= 3) == 2), paste("B1 needs >=3 cells in each condition for at least one cell type:", f))
    }
    if (any(c("B4", "B5") %in% active)) {
      tumor <- x$cell_type == "Tumor cells"
      required(any(tumor), paste("B4/B5 need Tumor cells or Malignant cells:", f))
      required(length(unique(x$condition[tumor])) == 2L, paste("B4/B5 need both conditions among tumor cells:", f))
    }
    if ("B4" %in% active) {
      required("pca" %in% Reductions(x), paste("B4 requires a pca reduction:", f))
      required(ncol(Embeddings(x, "pca")) >= 20L, paste("B4 requires at least 20 existing PCs:", f))
    }
    if ("A2" %in% active) {
      required("sample" %in% names(x[[]]), paste("A2 requires study-verified sample labels; do not invent replicate IDs:", f))
      required("SCT" %in% Assays(x), paste("The inherited A2 method expects normalized SCT assay data:", f))
      counts <- table(x$cell_type, factor(x$condition, levels = c("sensitive", "resistant")))
      required(all(colSums(counts >= 10) >= 2), paste("A2 needs >=2 cell types with >=10 cells each in each condition:", f))
      required("Tumor cells" %in% rownames(counts) && all(counts["Tumor cells", ] >= 10),
               paste("The inherited tumor-sender A2 figures need tumor cells in both conditions:", f))
    }
    rm(x)
  }
}

if ("A3" %in% steps) {
  files <- list.files(paths$spatial_data, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  required(length(files) > 0L, "A3 requires spatial Seurat objects in paths.spatial_data; the GSE104987 example is not spatial.")
  for (f in files) {
    x <- readRDS(f)
    required(inherits(x, "Seurat"), paste("A3 expected Seurat:", f))
    required(all(c("final_celltype", "group", "sample") %in% names(x[[]])),
             "A3 is a study-specific template requiring final_celltype, group, and sample metadata; see docs/module-a.md.")
    required(length(Images(x)) > 0L, "A3 requires spatial image/coordinate data.")
  }
}

needs_table <- function(consumer, producer, module, pattern) {
  if (consumer %in% steps && !producer %in% steps)
    required(length(list.files(file.path(PATH_RESULTS, module), pattern = pattern)) > 0L,
             paste(consumer, "requires outputs from", producer, "in", file.path(PATH_RESULTS, module)))
}
needs_table("B2", "B1", "B", "_deg\\.sig\\.csv$")
needs_table("B3", "B2", "B", "_enrichment_up\\.csv$")
needs_table("B6", "B1", "B", "_deg\\.sig\\.csv$")
needs_table("B7", "B6", "B", "_drug_enrichment\\.csv$")
needs_table("C7", "B1", "B", "_deg\\.csv$")
needs_table("C8", "C7", "C", "_deg_rbp_enrichment\\.csv$")
needs_table("C9", "B1", "B", "_deg\\.sig\\.csv$")
needs_table("C10", "C9", "C", "_miRNA\\.csv$")

if ("D3" %in% steps) {
  apa_path <- Sys.getenv("RESIST_APA_CONFIG")
  required(file.exists(apa_path), paste("APA configuration missing:", apa_path))
  apa <- yaml::read_yaml(apa_path)
  for (k in c("gff_path", "annotation_seurat_path", "sample_sheet")) {
    required(!is.null(apa[[k]]) && !grepl("/path/to/", apa[[k]], fixed = TRUE),
             paste("Replace", k, "in", apa_path, "with an actual HPC path; see TUTORIAL.md, Step 7."))
    apa[[k]] <- resist_resolve(apa[[k]], home)
    required(file.exists(apa[[k]]), paste("Missing APA input:", apa[[k]]))
  }
  samples <- read.csv(apa$sample_sheet, stringsAsFactors = FALSE)
  required(nrow(samples) > 0L && all(c("sample_id", "group", "tenx_dir", "txs_rds") %in% names(samples)),
           "APA sample sheet needs rows and sample_id, group, tenx_dir, txs_rds columns.")
  required(!anyNA(samples) && !any(samples$sample_id == ""), "APA sample sheet contains missing values.")
  for (k in c("tenx_dir", "txs_rds")) for (v in samples[[k]])
    required(file.exists(resist_resolve(v, home)), paste("Missing APA sample input:", v))
  cat("APA input paths found. Barcode mapping and TXS contents still require HPC validation.\n")
}
cat("\nPrerequisite checks passed. This is not an end-to-end scientific validation.\n")
cat("R version:", R.version.string, "\n")
for (p in packages) cat(p, as.character(packageVersion(p)), "\n")
cat("RESIST_OUTPUT_ROOT=", PATH_RESULTS, "\n", sep = "")
