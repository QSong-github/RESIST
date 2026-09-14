#!/usr/bin/env Rscript
# Install required R software, then check it using small in-memory objects.
# No study data, reference datasets, or RESIST analysis workflows are run.
# Usage: Rscript --vanilla config/setup/install_r_packages.R [--check] [--help]
args <- commandArgs(trailingOnly = TRUE)
if (any(args %in% c('--help', '-h'))) {
  cat('Usage: Rscript --vanilla config/setup/install_r_packages.R [--check]\n',
      'Default: install missing dependencies and the pinned CellChat revision, then check software.\n',
      '--check: check installed software only; no downloads or installations.\n',
      'Reports: data/results/setup/<UTC timestamp>-<process id>/\n', sep = '')
  quit(status = 0)
}
if (length(setdiff(args, '--check'))) stop('Unknown argument. Use --help.', call. = FALSE)
check_only <- '--check' %in% args
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1])
# Some Rscript builds encode spaces in --file as ~+~ in commandArgs().
if (!file.exists(script)) script <- gsub('~+~', ' ', script, fixed = TRUE)
home <- normalizePath(file.path(dirname(script), '../..'), mustWork = TRUE)

# CellChat's exact source is fixed; transitive packages are recorded after resolution.
cellchat_sha <- '75253cd0c9e68410e6e721a6d3a0419a1d7e358f'
cellchat_repo <- paste0('jinworks/CellChat@', cellchat_sha)
cran <- c('BiocManager', 'remotes', 'tidyverse', 'Seurat', 'yaml', 'jsonlite',
          'RColorBrewer', 'ggsci', 'ggrepel', 'cowplot', 'patchwork', 'pheatmap',
          'viridis', 'reshape2', 'scales', 'broom', 'msigdbr', 'treemap',
          'treemapify', 'packcircles', 'R.utils', 'devtools')
# Depends / Imports / LinkingTo from the pinned CellChat DESCRIPTION.
cellchat_cran <- c('dplyr', 'igraph', 'ggplot2', 'future', 'future.apply', 'pbapply',
                  'irlba', 'NMF', 'ggalluvial', 'stringr', 'svglite', 'Matrix',
                  'ggrepel', 'circlize', 'RColorBrewer', 'cowplot', 'RSpectra',
                  'Rcpp', 'RcppEigen', 'reticulate', 'scales', 'sna', 'reshape2',
                  'FNN', 'shape', 'magrittr', 'patchwork', 'colorspace', 'plyr',
                  'ggpubr', 'ggnetwork', 'plotly', 'shiny', 'bslib', 'collapse')
cellchat_bioc <- c('BiocNeighbors', 'BiocGenerics', 'ComplexHeatmap')
bioc <- unique(c(cellchat_bioc, 'SingleCellExperiment', 'SummarizedExperiment',
                 'GenomicRanges', 'rtracklayer', 'clusterProfiler', 'GSVA',
                 'GSEABase', 'fgsea', 'biomaRt', 'org.Hs.eg.db', 'org.Mm.eg.db',
                 'cmapR', 'cytolib', 'flowCore'))
cran <- unique(c(cran, cellchat_cran))
all_packages <- unique(c(cran, bioc, 'CellChat'))
minimum <- c(NMF = '0.23.0', circlize = '0.4.12')

# Use the active conda R library explicitly. Do not accidentally repair base R
# or satisfy requirements with packages from a personal library outside the env.
prefix <- Sys.getenv('CONDA_PREFIX')
if (nzchar(prefix)) {
  expected <- normalizePath(file.path(prefix, 'lib/R'), mustWork = FALSE)
  if (!identical(normalizePath(R.home()), expected))
    stop('Rscript is outside the active conda environment. Activate resist and check which Rscript.', call. = FALSE)
  lib <- normalizePath(.Library)
  .libPaths(lib, include.site = FALSE)
} else {
  lib <- .libPaths()[1]
}
cat('R: ', R.version.string, '\nR home: ', R.home(),
    '\nInstallation library: ', lib, '\n', sep = '')
if (!check_only && file.access(lib, 2L) != 0L)
  stop('The installation library is not writable: ', lib, call. = FALSE)

# Record all attempted checks, including a failing installation.
report_dir <- file.path(home, 'data/results/setup',
                       paste0(format(Sys.time(), '%Y%m%dT%H%M%SZ', tz = 'UTC'), '-', Sys.getpid()))
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
status <- 'failed'
checks <- character()
bioc_release <- NA_character_
record <- function() {
  cols <- c('Package', 'Version', 'LibPath', 'Built')
  pkg <- as.data.frame(installed.packages()[, cols, drop = FALSE], stringsAsFactors = FALSE)
  write.csv(pkg, file.path(report_dir, 'installed-packages.csv'), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(report_dir, 'sessionInfo.txt'))
  writeLines(c(paste('Status:', status), paste('Mode:', if (check_only) 'check' else 'install'),
               paste('R library:', lib), paste('Bioconductor:', bioc_release),
               paste('CellChat requested:', cellchat_repo), checks),
             file.path(report_dir, 'software-checks.txt'))
  cat('Software report: ', report_dir, '\n', sep = '')
}

# Test loading in a child R session: do not retain DLLs/namespaces that may need
# reinstalling in this process. Failed loads report the underlying error.
load_failures <- function(packages) {
  code <- paste0('p <- ', paste(deparse(packages, width.cutoff = 500L), collapse = ''),
    '; .libPaths(', paste(deparse(.libPaths()), collapse = ''), ', include.site=FALSE); ',
    'for (x in p) tryCatch(loadNamespace(x), error=function(e) ',
    'cat("RESIST_MISSING:",x,"\\n",conditionMessage(e),"\\n",sep=""))')
  out <- system2(file.path(R.home('bin'), 'Rscript'),
                 c('--vanilla', '-e', shQuote(code)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, 'status')) && attr(out, 'status') != 0L)
    stop('Package-load check could not start:\n', paste(out, collapse = '\n'), call. = FALSE)
  missing <- sub('^RESIST_MISSING:', '', grep('^RESIST_MISSING:', out, value = TRUE))
  if (length(missing)) cat(paste(out, collapse = '\n'), '\n')
  missing
}
needs_minimum <- function() {
  names(minimum)[vapply(names(minimum), function(p) {
    loc <- find.package(p, quiet = TRUE)
    !length(loc) || utils::packageVersion(p) < package_version(minimum[[p]])
  }, logical(1))]
}
cellchat_matches <- function() {
  loc <- find.package('CellChat', quiet = TRUE)
  if (!length(loc)) return(FALSE)
  desc <- read.dcf(file.path(loc, 'DESCRIPTION'))
  'RemoteSha' %in% colnames(desc) && identical(unname(desc[1, 'RemoteSha']), cellchat_sha)
}

main <- function() {
  options(timeout = max(600, getOption('timeout', 60)),
          repos = c(CRAN = 'https://cloud.r-project.org'))
  if (!check_only) {
    if (!length(find.package('BiocManager', quiet = TRUE)))
      install.packages('BiocManager', lib = lib)
    if (!requireNamespace('BiocManager', quietly = TRUE))
      stop('BiocManager installation failed; inspect the install output.', call. = FALSE)
    # Preserve the compatible release selected for this R installation. Do not
    # force the current Bioconductor release onto an older R minor version.
    bioc_release <<- as.character(BiocManager::version())
    options(repos = BiocManager::repositories(version = bioc_release))
    cat('Bioconductor repositories: ', bioc_release, '\n', sep = '')
    # Bioconductor prerequisites first, including BiocNeighbors. CRAN-only
    # install.packages("BiocNeighbors") cannot resolve this dependency.
    miss <- load_failures(bioc)
    if (length(miss)) BiocManager::install(miss, lib = lib, version = bioc_release,
                                         ask = FALSE, update = FALSE, force = TRUE)
    miss <- unique(c(load_failures(cran), needs_minimum()))
    if (length(miss)) install.packages(miss, lib = lib, dependencies = NA)
    miss <- load_failures(unique(c(cran, bioc)))
    if (length(miss)) stop('Dependencies still missing or unloadable: ',
                          paste(miss, collapse = ', '), call. = FALSE)
    if (length(needs_minimum())) stop('CellChat minimum dependency versions were not met.', call. = FALSE)
    if (!cellchat_matches() || length(load_failures('CellChat'))) {
      remotes::install_github(cellchat_repo, lib = lib, dependencies = NA,
                              upgrade = 'never', force = TRUE,
                              build_vignettes = FALSE, build_manual = FALSE)
    }
  }
  miss <- load_failures(all_packages)
  if (length(miss)) stop('Packages still missing or unloadable: ', paste(miss, collapse = ', '),
                        '. Run this installer without --check to install them.', call. = FALSE)
  if (length(needs_minimum())) stop('CellChat minimum dependency versions were not met.', call. = FALSE)
  if (!cellchat_matches()) stop('CellChat revision differs from the specified source. Run the installer without --check.', call. = FALSE)
  bioc_release <<- as.character(BiocManager::version())
  checks <<- c(checks, 'PASS: all required R package namespaces load; CellChat revision matches.')

  # Small API/compiled-code checks. These are not biological acceptance tests.
  set.seed(1)
  counts <- matrix(rpois(30L * 8L, lambda = 4), nrow = 30,
                   dimnames = list(paste0('gene', seq_len(30)), paste0('cell', seq_len(8))))
  seu <- Seurat::CreateSeuratObject(counts = counts)
  stopifnot(inherits(seu, 'Seurat'), ncol(seu) == 8L)
  checks <<- c(checks, 'PASS: Seurat creates a count-matrix object.')
  nn <- BiocNeighbors::findKNN(matrix(seq_len(24), nrow = 8), k = 2)
  stopifnot(identical(dim(nn$index), c(8L, 2L)), all(is.finite(nn$distance)))
  meta <- data.frame(cell_type = rep(c('typeA', 'typeB'), each = 4),
                     samples = factor(rep('sample1', 8)), row.names = colnames(counts))
  cc <- CellChat::createCellChat(log1p(counts), meta = meta, group.by = 'cell_type')
  stopifnot(inherits(cc, 'CellChat'), length(cc@idents) == 8L)
  checks <<- c(checks, 'PASS: BiocNeighbors compiled kNN and CellChat object construction.')
  required_api <- list(CellChat = c('createCellChat', 'subsetData', 'identifyOverExpressedGenes',
    'identifyOverExpressedInteractions', 'computeCommunProb', 'filterCommunication',
    'computeCommunProbPathway', 'aggregateNet', 'mergeCellChat', 'netVisual_circle'),
    GSVA = c('gsvaParam', 'gsva'), clusterProfiler = c('enrichGO', 'enrichKEGG', 'enricher'),
    cmapR = 'parse_gctx', biomaRt = 'getBM', rtracklayer = 'import')
  for (p in names(required_api)) for (f in required_api[[p]])
    stopifnot(is.function(getExportedValue(p, f)))
  gs <- GSVA::gsvaParam(log1p(counts), list(example_set = rownames(counts)[1:10]),
                        minSize = 1, maxSize = 30)
  stopifnot(methods::is(gs, 'gsvaParam'))
  gene_ranges <- GenomicRanges::GRanges('chr1', IRanges::IRanges(1:2, width = 10))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  stopifnot(length(gene_ranges) == 2L, ncol(sce) == 8L)
  checks <<- c(checks, 'PASS: module APIs, GSVA parameter object, genomic ranges and single-cell container.')
  cat(paste(checks, collapse = '\n'), '\n')
  status <<- 'passed'
  cat('R software checks passed. Study-data workflows still require HPC validation.\n')
}
tryCatch(main(), finally = record())
