# Input compatibility and validation; no biological relabeling beyond synonyms.
resist_read_seurat <- function(path) {
  seu <- readRDS(path)
  if (!inherits(seu, "Seurat")) stop("Expected a Seurat object: ", path)
  md <- seu[[]]
  candidates <- intersect(c("cell_type", "celltype"), names(md))
  if (!length(candidates)) stop("Missing celltype/cell_type metadata in ", path)
  canonical <- function(v) {
    v <- as.character(v)
    v[v == "Malignant cells"] <- "Tumor cells"
    v
  }
  ct <- canonical(md[[candidates[1]]])
  if (length(candidates) > 1L && !identical(ct, canonical(md[[candidates[2]]])))
    stop("Conflicting celltype and cell_type annotations in ", path)
  if (anyNA(ct) || any(!nzchar(ct))) stop("Missing cell-type labels in ", path)
  if (!"condition" %in% names(md)) stop("Missing condition metadata in ", path)
  cond <- as.character(md$condition)
  if (anyNA(cond) || !all(cond %in% c("sensitive", "resistant")))
    stop("condition must contain only sensitive/resistant; inspect study metadata: ", path)
  # Drop the alias in memory to avoid duplicate column names in downstream code.
  if ("celltype" %in% names(md)) seu$celltype <- NULL
  seu$cell_type <- ct
  seu
}
