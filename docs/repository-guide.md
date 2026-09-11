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
| `data/reference_bundle.tar.gz` | Compressed supplied references |
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

## Release identity and traceability

`VERSION` is the source of the package label displayed by all four launchers
and saved in `run.json`. Read it without loading data or R packages:

```bash
bash run_A.sh --version
```

`CITATION.cff` and `config/steps.json` carry matching release metadata.
Reference-database releases, software versions, and sequencing chemistry
identifiers are independent of the RESIST package label and must retain
their exact values. Public usage examples use a generic `RESIST` directory;
substitute the actual location of your extracted package.

The [migration map](migration.tsv) traces each original supplied file to its
current location using `source_path`, `current_path`, `source_sha256`, and
`current_sha256`. A target written as `archive.tar.gz::member` refers to a
file inside an archive. Retired wrappers and empty placeholders have no
current target and carry an explicit reason. Original filenames in the
source column remain unchanged for provenance.

The [source-history archive](source-history.tar.gz) preserves earlier layout,
migration, release, and validation records. Its contents describe historical
states and are not current execution instructions. Use the active guides in
this directory and the [release notes](release-notes.md) for this package.

`SHA256SUMS.txt` covers every other distributed file, including the archives
and migration map. Verify an untouched extracted package on the HPC with
`sha256sum -c SHA256SUMS.txt` before editing its configuration.
