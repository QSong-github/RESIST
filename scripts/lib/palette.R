# =============================================================================
# RESIST | scripts/lib/palette.R
# -----------------------------------------------------------------------------
# Shared colour definitions. These are the palettes used across every RESIST
# figure and across the web interface, so changing them here changes the whole
# resource consistently.
#
# Migrated verbatim from the original `scripts/color.R`.
# =============================================================================

#' Cell-type palette
#'
#' Sixteen qualitative colours. Position 1 is reserved for "Tumor cells", which
#' the characterisation entry points place first in the factor level order; a
#' dataset without a tumour compartment drops the first colour so that the rest
#' of the assignment is unchanged.
RESIST_COL_CELLTYPE <- c(
  "#E377C2FF",
  "#FF7F0EFF",
  "#2CA02CFF",
  "#9467BDFF",
  "#8C564BFF",
  "#FF9896FF",
  "#BCBD22FF",
  "#17BECFFF",
  "#AEC7E8FF",
  "#FFBB78FF",
  "#98DF8AFF",
  "#C5B0D5FF",
  "#C49C94FF",
  "#F7B6D2FF",
  "#DBDB8DFF",
  "#9EDAE5FF"
)

#' Condition palette: sensitive vs resistant
RESIST_COL_CONDITION <- c("#E4C66F", "#5B9BD5")

#' Return a named palette
#'
#' Entry points that mutate `col` inside a per-dataset loop call
#' `resist_palette("celltype")` at the top of each iteration to restore the
#' full vector. This replaces the original pattern of re-sourcing `color.R`
#' from inside the loop.
#'
#' @param which "celltype" or "condition".
resist_palette <- function(which = c("celltype", "condition")) {
  which <- match.arg(which)
  switch(which,
         celltype  = RESIST_COL_CELLTYPE,
         condition = RESIST_COL_CONDITION)
}
