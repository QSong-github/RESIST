# Release notes — v3_02

## Response to the organizational review

| Concern | Change |
|---|---|
| Too many top-level folders | Reduced 12 to 4: `config`, `scripts`, `data`, `docs` |
| Mixed quick-start commands | Four root entry points, A/B/C/D, each with an explicit default scope |
| No worked example | Root `TUTORIAL.md`, using the observed GSE104987 file and metadata |
| Unclear results | Quick-start product summaries, exact example filenames, and `docs/outputs.md` |
| Large references on a laptop | Supplied 75.3 MB reference bundle compressed to 16.8 MB; large optional downloads deferred to HPC |
| Difficult to trace changes | Original step IDs preserved; file-level migration/checksum map and static validation record |

## Organization

The four analysis modules now live under `scripts/A`–`D`, with analysis scripts
flattened within each module. Shared R code and manual utilities are grouped under
`scripts/lib` and `scripts/utils`. Installation and SLURM material live under
`config`. Small example metadata and sequences are in `data/templates`. Inputs,
unpacked references, and outputs are created on the HPC only when needed.

The source directory supplied by the user was not modified. No original analysis
entry point was deliberately removed. Missing C2 and TF code remain clearly
identified; empty placeholder folders were not retained to imply functionality.

## Execution and documentation repairs

- Updated root discovery, shared-library paths, documentation links, and result
  paths for the new layout. The launcher resolves a config path from the caller's
  directory and runs scripts with a stable package root.
- Added selected-step prerequisite checks and logs that distinguish failures from
  successful commands that produce no new files. These checks are intended for
  later HPC use; they do not claim comprehensive data validation.
- Added metadata alias handling for A1/A2/B1/B4/B5 without changing saved inputs.
- Moved A2 computation inside its dataset loop, corrected its input source from
  the spatial directory to the configured single-cell directory, and relocated
  the custom CellChat helper. These changes require HPC runtime review.
- Replaced A3's undefined plotting-palette variable with the shared palette;
  the study-specific spatial interface remains documented.
- Corrected C8 plot paths to join the output directory and filename, and corrected
  its mouse basename comparison. Made C7's inherited Holm adjustment explicit
  and labeled that method in the C8 legends; no BH substitution was made.
- Connected D3's annotation loading to its existing metadata helper and resolved
  GTF/annotation/10x/TXS relative paths against the package root. Ambiguous RData
  containers now fail instead of silently selecting among several Seurat objects.
- Added Python/jsonlite to the launcher's installation requirements and corrected
  dependency declarations for optional R tools. Package-resolution/runtime testing
  remains for the HPC.
- Replaced overstated or inaccurate documentation: A1 plots cluster composition;
  B4/B5 currently export figures rather than score CSVs; B1 is cell-level;
  C7's legacy FDR-named columns hold Holm-adjusted values; D5 is not a generic
  local script; D2 and LINCS parsing inconsistencies are disclosed.

## Deliberately preserved scientific choices

The analysis definitions, cell-count/significance thresholds, dataset exclusion
lists, default active-assay behavior, human-reference branches, and C8 plot gates
remain inherited. This package does not perform a new methodological validation,
recalculate published statistics, harmonize cohort designs, infer missing C2
filtering, or implement TF motifs.

## Validation boundary

Per the user's request, the analysis workflows and large external reference
downloads were not run for this release. An earlier metadata-only inspection of
the shared object established the tutorial's input facts; the temporary large
example file was then removed to conserve local storage. No real or synthetic
analysis results are presented as completed output. See [validation.md](validation.md)
for the exact static checks performed and the remaining HPC acceptance work.
