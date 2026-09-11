# =============================================================================
# RESIST | scripts/lib/config.R
# -----------------------------------------------------------------------------
# Configuration loading and path resolution.
# =============================================================================

#' Read `config/config.yaml`
#'
#' @param home Repository root.
#' @param file Optional explicit path; overrides the default location and the
#'   `RESIST_CONFIG` environment variable.
#' @return Parsed configuration as a nested list.
resist_load_config <- function(home, file = NULL) {
  if (is.null(file)) {
    file <- Sys.getenv("RESIST_CONFIG", unset = "")
    if (!nzchar(file)) file <- file.path(home, "config", "config.yaml")
  }
  if (!file.exists(file)) {
    stop("RESIST configuration not found: ", file, call. = FALSE)
  }
  yaml::read_yaml(file)
}

#' Resolve a configured path against the repository root
#'
#' Absolute paths are returned unchanged, which lets a deployment point
#' `paths.data` or `paths.ref_data` at a shared filesystem location without
#' editing any script.
#'
#' @param p    Path from the configuration.
#' @param home Repository root.
resist_resolve <- function(p, home = RESIST_HOME) {
  if (is.null(p) || !nzchar(p)) stop("empty path in configuration", call. = FALSE)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) normalizePath(p, mustWork = FALSE)
  else normalizePath(file.path(home, p), mustWork = FALSE)
}

#' Absolute path to a reference file declared under `reference_files`
#'
#' @param key Key in `config.yaml -> reference_files`.
#' @param must_exist Stop with an actionable message when the file is absent.
resist_ref <- function(key, must_exist = TRUE) {
  val <- CFG$reference_files[[key]]
  if (is.null(val)) {
    stop("unknown reference_files key in config.yaml: '", key, "'", call. = FALSE)
  }
  p <- file.path(PATH_REF, val)
  if (must_exist && !file.exists(p)) {
    stop("required reference file is missing:\n  ", p,
         "\nSee docs/references.md for how to obtain it.", call. = FALSE)
  }
  p
}

#' Per-module output directory, created on first use
#'
#' @param module Module directory name, e.g. "B".
#' @param ...    Optional sub-directory components.
resist_results_dir <- function(module, ...) {
  d <- file.path(PATH_RESULTS, module, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' A configured dataset hold-out list
#'
#' @param key Key in `config.yaml -> datasets`.
resist_datasets <- function(key) {
  val <- CFG$datasets[[key]]
  if (is.null(val)) {
    stop("unknown datasets key in config.yaml: '", key, "'", call. = FALSE)
  }
  unlist(val, use.names = FALSE)
}

#' Input Seurat objects for a module
#'
#' @param dir     Directory to scan; defaults to the single-cell input directory.
#' @param exclude Character vector of file names to skip.
#' @return Full paths, sorted, of `*.rds` / `*.RDS` objects.
resist_input_objects <- function(dir = PATH_DATA, exclude = character()) {
  if (!dir.exists(dir)) {
    stop("input directory does not exist: ", dir,
         "\nSee TUTORIAL.md and docs/input-data.md.", call. = FALSE)
  }
  f <- list.files(dir, pattern = "\\.rds$", full.names = TRUE, ignore.case = TRUE)
  f <- sort(f)
  if (length(exclude)) f <- f[!basename(f) %in% exclude]
  if (!length(f)) {
    warning("no Seurat objects found in ", dir, call. = FALSE)
  }
  f
}
