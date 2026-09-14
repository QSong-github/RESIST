# Package update notes

The package now presents four module launchers under `scripts/`, uses the supplied
PDF as the workflow figure, and installs CellChat’s required dependencies through
the correct R repositories. The README describes the products of A, B, C and D
without exposing numbered script identifiers in the quick start.

## Workflow and organization

- Moved the four launchers to `scripts/run_A.sh` through `scripts/run_D.sh` and
  updated the tutorial, module guides, SLURM template and command examples.
- Removed the separate input-preparation check script from active execution.
  The shared dispatcher reads configuration, records software, and runs each
  selected analysis in its own R process. Dataset checks remain within the
  analysis code and tutorial input checkpoints.
- Retired the manual result-renaming utility, preserving its source in the
  historical archive. Moved the useful annotation-table helper beside Module D
  as `scripts/D/prepare_cell_annotations.sh`; removed the extra utilities folder.
- Included `docs/figures/Fig1.pdf` unchanged and rendered it to a PNG for GitHub.
  The README displays the PNG and links to the original PDF.
- Retained four top-level directories and the concise key-tools table. Public
  guides use a generic package name; release metadata remains in `VERSION`,
  `CITATION.cff`, and `config/steps.json`.

## Installation and reproducibility

The R installer explicitly includes CellChat’s CRAN and Bioconductor dependencies,
including `BiocNeighbors`, `BiocGenerics`, and `ComplexHeatmap`. It installs the
Bioconductor prerequisites first using a release compatible with the active R,
then resolves CRAN dependencies before installing a fixed CellChat commit through
`remotes`. Compilers and common compiled dependencies are included in the conda
specification. No large reference database is downloaded by the installer.

In an active conda environment, both installation and module execution verify
that R belongs to that environment. The installer prints its library destination.
It stops when required packages cannot load, the CellChat revision differs, or
one of its small software checks fails. `--check` belongs to the installer and
repeats only these software checks without installation.

Setup reports contain package versions, library paths, the CellChat source
revision, the R session, and check status. Module runs keep software metadata,
configuration snapshots, script/shared-helper hashes, step logs, and hashes of
new or changed outputs. Input/reference checksums and separate external-tool
versions must be retained with the HPC run.

## Scientific scope

Analysis expressions, default selections, output filenames, statistical
contrasts, thresholds, and reference data remain unchanged. The shared example
supports the core A–C walkthrough. CellChat needs a separate dataset with multiple
cell types and sample labels; D requires independently prepared APA inputs.

The missing VCF conversion and TF motif implementation, the cohort-specific
spatial/APA assumptions, and the LINCS input-format requirements remain documented.
The figure depicts the wider resource; it does not establish that every panel is
an automated branch of the quick start.

The [migration map](migration.tsv) accounts for all original supplied files,
including the unchanged compressed references and archived support code.
Historical records are preserved in [source-history.tar.gz](source-history.tar.gz).
See [validation.md](validation.md) for the local checks and the work deferred to
HPC acceptance. Installation specifications and small API checks do not establish
numerical reproduction of the original biological analyses.
