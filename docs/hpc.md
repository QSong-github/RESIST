# HPC preparation and submission

Transfer the compact package, install dependencies in an HPC environment, and
place large data/reference files in project or shared storage. Nothing in the
quick start requires storing those large files on a laptop.

## Prepare once

1. Extract the package in persistent project space.
2. Create and activate the R/Python environment using tutorial Step 1.
3. Download the example and unpack references on the HPC.
4. Use `config/example.yaml` for the worked example; copy `config/config.yaml`
   for another cohort and edit only that cohort's paths/settings.
5. Confirm the R installer passes its software checks; review the tutorial’s input
   checkpoint on an appropriate compute node.

The supplied software specifications are not a tested lockfile. Capture the
resolved installation for the run:

```bash
conda env export > /your/HPC/project/resist-environment-resolved.yml
conda list --explicit > /your/HPC/project/resist-conda-explicit.txt
Rscript -e 'writeLines(capture.output(sessionInfo()), "/your/HPC/project/resist-sessionInfo.txt")'
```

## Submit one module at a time

From the package root, with the environment activated:

```bash
sbatch --export=ALL,RESIST_HOME="$PWD" config/slurm/submit_module.sbatch \
  A --config config/example.yaml
sbatch --export=ALL,RESIST_HOME="$PWD" config/slurm/submit_module.sbatch \
  B --config config/example.yaml
```

Submit C **after B succeeds**. To enforce that relationship in SLURM:

```bash
job_b=$(sbatch --parsable --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch B --config config/example.yaml)
sbatch --dependency="afterok:${job_b%%;*}" --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch C --config config/example.yaml
```

The example submission template requests one node, four CPUs, 16 GB RAM, and four
hours as initial placeholders, not measured requirements. Edit account/partition
settings for your cluster. It inherits the activated environment and does not
hard-code an institution-specific R module. Check that both `Rscript` and `python3`
are available inside the job.

B4 memory grows quadratically with tumor-cell count. B6 needs a much larger
allocation; the source comments describe tens of GB and include a 90 GB
submission request. This is not a measured requirement for your dataset. An example override is:

```bash
sbatch --mem=90G --time=12:00:00 --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch B --config config/example.yaml --steps B6,B7
```

Measure actual peak memory and elapsed time with the scheduler and adjust the
request. The `--cpus-per-task` allocation does not itself make every R step
parallel. Submit D3 only once its independent cohort inputs are prepared; GPU
structural prediction should use the institution's separate GPU workflow.

## Paths, logs, and repeated runs

Every module launcher can be called from another working directory by absolute path.
Use an absolute `--config` path in that case. The wrapper passes the package root
to the scripts; configured relative data paths still resolve against that root.
Direct R entry points also retain their root-discovery bootstrap, but direct
invocation bypasses the launcher's software records and provenance logs.

Each runner logs to `<results>/logs/` and writes analysis products to
`<results>/A`, `/B`, `/C`, or `/D`. Configure a new results root for a new analysis
comparison. Same-name biological output files can be overwritten by repeated
steps; only the run log directory is always unique. Downstream steps scan their
configured results directories, so stale tables should not be mixed between
cohorts.

## First-run acceptance

For the GSE104987 example, check the observed input metadata, inspect A's UMAP,
verify B's CSV schemas and effect-size direction, and inspect C's adjusted-p-value
method and plotting gates. Then enable optional steps individually. The installer’s
software checks verify package loading and selected APIs; numerical validation
requires execution on the appropriate study inputs. Record any package
compatibility or scientific-method changes made during HPC testing before using
results in a manuscript.

## CellChat installation

`BiocNeighbors` is a **Bioconductor** package. Installing it only from CRAN can
report “not available for this version of R,” even when a compatible Bioconductor
package exists. CellChat also requires `BiocGenerics` and `ComplexHeatmap`,
alongside CRAN dependencies such as `NMF`, `circlize`, `igraph`, `RcppEigen`,
`RSpectra`, `ggalluvial`, and `svglite`. The installer includes CellChat’s direct
required dependencies and resolves their dependencies through the appropriate
repositories before installing CellChat. See the [pinned CellChat dependency
manifest](https://github.com/jinworks/CellChat/blob/75253cd0c9e68410e6e721a6d3a0419a1d7e358f/DESCRIPTION)
and [Bioconductor installation guidance](https://bioconductor.org/install/).

For the environment in which the earlier installation failed, enter this
package’s directory and rerun its updated installer:

```bash
conda activate resist
which Rscript
Rscript --vanilla config/setup/install_r_packages.R
```

Do not rerun `conda env create` merely to repair a missing R package in an
existing environment. The installer checks that R belongs to the active conda
environment and uses its `lib/R/library` directory. A path under
`.../envs/resist/lib/R/library/CellChat` is the expected installation destination.
The separate conda `pkgs/` directory is a download/extraction cache.

The supplied environment requests R 4.5. Bioconductor 3.21 and 3.22 support
R 4.5; newer Bioconductor releases may require newer R. The installer retains
the compatible release selected by `BiocManager` and records it, rather than
forcing the latest release. See the official [R/Bioconductor compatibility
table](https://bioconductor.org/about/release-announcements/).

CellChat is installed from commit
`75253cd0c9e68410e6e721a6d3a0419a1d7e358f` using
[`remotes::install_github`](https://remotes.r-lib.org/reference/install_github.html).
This replaces the deprecated devtools wrapper and avoids following a moving
GitHub branch. If an installed CellChat has a different or unrecorded source
revision, the installation command replaces it with this revision; `--check`
only reports the mismatch. Dependencies already meeting requirements are not
bulk-updated. The environment remains an installation specification until the
HPC run produces a resolved software record.

Repeat only the software checks after installation:

```bash
Rscript --vanilla config/setup/install_r_packages.R --check
```

This checks required namespaces, CellChat source identity, selected module APIs,
BiocNeighbors nearest-neighbor computation, and small Seurat/CellChat/GSVA and
single-cell/genomic-range objects. It does not run communication inference,
pathway queries, an APA cohort, or any large analysis. Review
`data/results/setup/` for success/failure and exact installed versions.

The conda environment includes compilers and common compiled dependencies for
building remaining R packages on the HPC. If compilation fails, retain the first
compiler/linker error and check the cluster’s supported toolchain. The
`invalid uid/gid … nobody` tar warnings do not explain a missing BiocNeighbors
package. A nonzero installation exit status or “Packages still missing” must be
resolved before enabling that analysis.

The earlier X11 `ClobberError` and package-cache `SafetyError` are separate
conda issues. Adding explicit CellChat dependencies does not establish that
those conflicts are resolved. Do not suppress file conflicts or remove a shared
HPC package cache as an installation shortcut; retain the solver/transaction
log and resolve the conflicting packages or cache with the cluster administrator.

## Preserve the resolved environment

Keep the conda export and explicit package list together with the R installer’s
package inventory and source revision. Conda exports alone may not describe R
packages installed from CRAN, Bioconductor, or GitHub. External tools such as
cellSNP-lite, scUTRquant, OptiType, ColabFold and PyMOL use their own environments;
record their versions and commands when those branches are used. The shared R
installer does not claim to install or validate those separate tools.
