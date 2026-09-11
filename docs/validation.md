# Validation record — v3_02

**Scope:** organization, documentation, static syntax, path/manifest consistency,
and packaging integrity. No full analysis workflows, SLURM jobs, or large external
reference downloads were executed for this release.

## Checks completed

| Check | Result |
|---|---|
| Top-level layout | 12 original folders consolidated into 4 |
| Python parsing | 6 files parsed with the Python AST parser |
| R parsing | 25 files parsed; analysis expressions not evaluated |
| Shell parsing | 10 shell/SLURM files checked with `bash -n` |
| Root runner interfaces | A/B/C/D `--list` and `--help` exercised from `/tmp`, using the package path containing spaces |
| Invalid requests | Cross-module step and reversed B2/B1 order rejected before R execution |
| Step manifest | 15 selectable R entries resolve to existing scripts |
| Local documentation | Relative file and heading links checked; the labeled historical v2 map is excluded |
| Source traceability | 188 original files accounted for in the migration map, including explicitly retired wrappers/placeholders |
| Bundled references | 119 files preserved byte-for-byte in the compressed reference archive |
| Package contents | No populated example-data folder, expanded reference tree, or generated biological results |

The exact static-check output is saved in [static-checks.txt](static-checks.txt).
Environment/YAML and archive checks are also recorded there when performed.

## Observed example metadata

Before the user clarified that execution would occur later on the HPC, the shared
`GSE104987_seurat_afterAnno.RDS` was downloaded temporarily and inspected as a
Seurat object. It contained 2,669 cells, with 1,597 resistant and 1,072 sensitive,
all annotated as malignant cells; RNA/SCT assays; UMAP/PCA reductions; and 50 PCA
components. This inspection did not run any RESIST biological analysis step.
The temporary approximately 291 MB file was removed to conserve local storage.

The observed metadata informs the tutorial's input checkpoint. There is no
reference set of DEG counts, enrichment results, figures, or runtimes from a
completed v3_02 run.

## Deferred to HPC acceptance testing

- Resolve the specified package environment and record its exact versions.
- Download and inspect the current shared example; preserve its checksum.
- Execute A1, then B1/B4, then C7/C8, inspecting each table and figure.
- Check empty/significance-filtered results and active-assay behavior.
- Enable B2/B3/B5/C9/C10 only after the relevant dependencies are prepared.
- Validate A2/A3 and D3 on appropriate multi-cell-type/spatial/APA inputs.
- Review D2 read assignments, recover or validate C2, and define the missing TF
  workflow before claiming those analyses are reproducible.
- Obtain LINCS only when needed, validate its parser/schema, and measure HPC memory.
- Compare numerical outputs with an authoritative reference analysis before
  describing this release as reproducing published results.

Static success demonstrates a coherent package structure and syntactically
readable code. It does not establish package-version compatibility, correct
barcode mapping, numerical equivalence, or scientific validity of every method.
