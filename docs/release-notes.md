# Release notes — v4

This is a complete RESIST package containing the four analysis modules, shared
code, configuration, compressed references, workflow figure, and documentation.
The public README uses the resource name without a release-number banner.

## Changes in this release

| Area | Current behavior |
|---|---|
| Release identity | `VERSION`, `CITATION.cff`, and `config/steps.json` consistently identify the package as `v4` |
| Launcher reporting | All four launchers read `VERSION` for their execution banner and `run.json`; `--version` displays it without loading analysis inputs or R packages |
| Usage documentation | Tutorial, module guides, data instructions, and SLURM messages use the current layout and generic package paths |
| Workflow overview | The supplied main workflow figure remains prominently embedded in the README |
| Quick start | Four A/B/C/D commands state their exact default steps and expected products |
| Worked example | The GSE104987 tutorial retains download, configuration, input inspection, A–C commands, output checkpoints, and a separate APA route |
| APA descriptions | Script comments and the Module D guide accurately describe the implemented annotation loading and remaining cohort checks |
| Structure image | D6 writes `MHC_Peptide_Interaction.png`; script, module guide, and output catalog agree |
| Historical records | Earlier layout, migration, release, and validation records are preserved in `docs/source-history.tar.gz` |
| File traceability | `docs/migration.tsv` uses stable source/current column names; current hashes and package checksums are regenerated |

The package retains four top-level folders: `config`, `scripts`, `data`, and
`docs`. The default route is A1 → B1/B4 → C7/C8. D defaults to D3 and requires its
own prepared APA inputs. Optional analyses retain their original step IDs and
remain documented in the module guides.

## Compatibility and provenance

Analysis definitions, thresholds, statistical adjustments, reference contents,
and default step selections are unchanged from the preceding package. The D6
image filename is the only analysis-output naming change; update downstream
commands that expect its previous filename. The original source directory and
previous local package remain untouched.

Reference-database releases, sequencing chemistry, and external tool/model
versions are independent of the RESIST release. Identifiers such as `10xv3`,
`miRDB_v6.0`, Hallmark `v2024.1`, and ColabFold/AlphaFold model versions retain
their correct values. Historical source filenames remain in provenance records
and source-origin comments. They must not be relabeled as RESIST versions.

The migration map accounts for all 188 original supplied files, including 119
reference files stored unchanged inside the compressed reference bundle. Retired
wrappers and empty placeholders are identified explicitly. Historical documents
inside the source-history archive describe earlier states; use the active guides
for execution. The root `SHA256SUMS.txt` covers all other distributed files.

## Validation and HPC testing

This release was checked for syntax, launcher interfaces, version consistency,
documentation links, source traceability, and packaging integrity. The analysis
workflows, environment installation, SLURM submissions, and large external
reference downloads were not run. See the [validation record](validation.md)
and [HPC acceptance guidance](hpc.md#first-run-acceptance).

The package retains the documented implementation boundaries: C2 and TF motif
code are absent; A2/A3/D3 need suitable independent inputs; D2 read roles require
review; LINCS parsing and dependency compatibility require HPC validation.
No biological outputs or claims of numerical reproduction are supplied.
