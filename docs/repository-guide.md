# Repository guide

Read the root [README](../README.md) and [tutorial](../TUTORIAL.md) first.
The four root runners share one Python dispatcher and one R preflight script.

| Location | Responsibility |
|---|---|
| `config/config.yaml` | Production paths, reference keys, dataset exclusions |
| `config/example.yaml` | Same schema, isolated GSE104987 input/result locations |
| `config/apa_config.yaml`, `config/apa_samples.csv` | Separate APA cohort inputs |
| `config/steps.json` | Selectable step IDs, explicit defaults, dependency/output descriptions |
| `config/setup/` | Environment and R/Python installation specifications |
| `config/slurm/submit_module.sbatch` | One-module SLURM submission template |
| `scripts/A/` … `scripts/D/` | Analysis entry points, flattened within each module |
| `scripts/lib/` | R configuration, metadata compatibility, palettes, common I/O, CellChat plotting helper |
| `scripts/utils/` | Legacy manual utilities; not invoked automatically |
| `data/reference_bundle.tar.gz` | Compressed references from v2-2 |
| `data/templates/` | Small sequencing/APA/structure examples, not production inputs |
| `docs/` | Human-readable method, output, HPC, and release records |

Default steps: A1; B1/B4; C7/C8; D3. These choices make the A–C example route
explicit while keeping additional analyses accessible with `--steps`. There is no
`run_all.sh` with ambiguous coverage.

The shared R initialization resolves paths and exposes `PATH_DATA`, `PATH_SPATIAL`,
`PATH_REF`, `PATH_RESULTS`, `RESULTS_A` through `RESULTS_D`, `resist_ref()`, and
output/palette helpers. Result paths are created when used. `RESIST_CONFIG` and
`RESIST_APA_CONFIG` support direct-script deployments; the root launchers set them
from the command-line options.

To add a step, retain the bootstrap used by the neighboring script, write outputs
under the correct module result directory, add its contract to `config/steps.json`,
and document prerequisites and products. Do not introduce cohort-specific absolute
paths into an analysis script. Record scientific method changes separately from
path/organization changes.

The historical v2 mapping is retained as [migration-v2.md](migration-v2.md).
It describes an earlier release and is not a current execution guide. The v3 map
is [migration.tsv](migration.tsv); [release notes](release-notes.md) explain the
behavioral changes.
