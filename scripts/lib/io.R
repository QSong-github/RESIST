# =============================================================================
# RESIST | scripts/lib/io.R
# -----------------------------------------------------------------------------
# Output helpers shared by all modules.
# =============================================================================

#' Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Write a figure as PDF and PNG side by side
#'
#' RESIST publishes a vector copy for figure assembly and a raster copy for the
#' web interface. `expr` is re-evaluated once per device, which keeps
#' base-graphics and grid-based figures (ComplexHeatmap, CellChat) behaving the
#' same way as ggplot objects.
#'
#' @param expr        Expression that draws the figure, or a ggplot/grid object.
#' @param path_prefix Output path without extension.
#' @param width,height Device size in inches.
#' @param dpi         Raster resolution; defaults to `figures.dpi` in the config.
#' @return The prefix, invisibly.
resist_save_figure <- function(expr, path_prefix, width = 8, height = 5.5,
                               dpi = NULL) {
  if (is.null(dpi)) dpi <- CFG$figures$dpi %||% 300
  q   <- substitute(expr)
  env <- parent.frame()

  draw <- function() {
    obj <- eval(q, env)
    if (inherits(obj, c("ggplot", "patchwork", "grob", "gtable"))) print(obj)
    invisible(NULL)
  }

  grDevices::pdf(paste0(path_prefix, ".pdf"), width = width, height = height)
  draw()
  grDevices::dev.off()

  grDevices::png(paste0(path_prefix, ".png"), width = width, height = height,
                 units = "in", res = dpi)
  draw()
  grDevices::dev.off()

  invisible(path_prefix)
}

#' Create a directory if needed and return it
resist_dir <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Timestamped progress message
resist_log <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
}
