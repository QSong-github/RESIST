# =============================================================================
# RESIST | scripts/lib/init.R
# -----------------------------------------------------------------------------
# Single entry point for the shared R layer. Every module script sources this
# file through the standard bootstrap block and, in return, receives:
#
#   RESIST_HOME  repository root (character)
#   CFG          parsed `config/config.yaml` (list)
#   PATH_DATA    single-cell input directory
#   PATH_SPATIAL spatial input directory
#   PATH_REF     reference-data directory
#   RESULTS_A..D per-module output paths (created by analysis scripts)
#   col, col1    cell-type and condition colour palettes
#   plus the helpers defined in io.R, config.R and palette.R.
#
# Output directories are created by the analysis that uses them.
# =============================================================================

if (!exists(".RESIST_INITIALISED", envir = globalenv())) {

  if (!exists(".resist_home", envir = globalenv())) {
    stop("scripts/lib/init.R must be sourced through the RESIST bootstrap block; ",
         ".resist_home is not defined.", call. = FALSE)
  }

  RESIST_HOME <- normalizePath(get(".resist_home", envir = globalenv()),
                               mustWork = TRUE)

  source(file.path(RESIST_HOME, "scripts", "lib", "packages.R"))
  source(file.path(RESIST_HOME, "scripts", "lib", "config.R"))
  source(file.path(RESIST_HOME, "scripts", "lib", "palette.R"))
  source(file.path(RESIST_HOME, "scripts", "lib", "io.R"))

  source(file.path(RESIST_HOME, "scripts", "lib", "input.R"))

  CFG <- resist_load_config(RESIST_HOME)

  PATH_DATA    <- resist_resolve(CFG$paths$data,         RESIST_HOME)
  PATH_SPATIAL <- resist_resolve(CFG$paths$spatial_data, RESIST_HOME)
  PATH_REF     <- resist_resolve(CFG$paths$ref_data,     RESIST_HOME)
  PATH_RESULTS <- resist_resolve(CFG$paths$results,      RESIST_HOME)

  RESULTS_A <- file.path(PATH_RESULTS, "A")
  RESULTS_B <- file.path(PATH_RESULTS, "B")
  RESULTS_C <- file.path(PATH_RESULTS, "C")
  RESULTS_D <- file.path(PATH_RESULTS, "D")

  if (!is.null(CFG$project$seed)) set.seed(CFG$project$seed)

  col  <- resist_palette("celltype")
  col1 <- resist_palette("condition")

  .RESIST_INITIALISED <- TRUE

  message(sprintf("[RESIST] repository root : %s", RESIST_HOME))
  message(sprintf("[RESIST] input data      : %s", PATH_DATA))
  message(sprintf("[RESIST] reference data  : %s", PATH_REF))
  message(sprintf("[RESIST] results         : %s", PATH_RESULTS))
}
