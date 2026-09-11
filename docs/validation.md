# Validation record

**Scope:** documentation, release identity, static syntax,
launcher interfaces, path/manifest consistency, and packaging integrity.
No analysis workflows, environment installation, SLURM jobs, or large external
reference downloads were executed during this package update.

## Checks completed

| Check | Result |
|---|---|
| Release consistency | `VERSION`, citation metadata, and the step manifest agree; all A–D `--version` responses identify the same package |
| Version cleanup | Active guides use no RESIST release-number labels; package metadata agrees; independent tool/reference/chemistry identifiers retain their exact values |
| README clarity | No numbered script IDs; A–D descriptions state default outputs, file formats, destinations, and input dependencies |
| Analysis preservation | All analysis and launcher scripts are byte-identical to the preceding package |
| Python parsing | 6 files parsed with the Python AST parser; analysis scripts not imported or executed |
| R parsing | 25 files parsed without evaluating analysis expressions |
| Shell parsing | 10 shell/SLURM files checked with `bash -n` |
| Root runner interfaces | A/B/C/D `--list`, `--help`, and `--version` exercised from `/tmp`, including a package path containing spaces |
| Invalid requests | Cross-module steps and reversed B2/B1 order rejected before R execution |
| Step manifest | 15 selectable R entries resolve to existing scripts; default selections are unchanged |
| Structured files | Environment/configuration/template YAML, citation metadata, and step JSON parse successfully |
| Local documentation | Active Markdown file/heading links and HTML image/navigation links resolve |
| Main workflow figure | Identical to the supplied image |
| Source traceability | 188 original files accounted for in the migration map, with current target hashes verified |
| Bundled references | 119 files preserved byte-for-byte inside the compressed reference archive |
| Historical records | Earlier documentation retained in a separate provenance archive |
| Package contents | Exactly 4 top-level folders; no populated example, input, expanded-reference, or biological-result directories |

The detailed local check output is saved in [static-checks.txt](static-checks.txt).
The complete ZIP is checked against the package files, including CRC integrity
and the root checksum manifest. Its SHA-256 is provided in the adjacent
`.zip.sha256` file. To verify the extracted package before changing configuration:

```bash
sha256sum -c SHA256SUMS.txt
```

## Example metadata used by the tutorial

During an earlier documentation-preparation step, the shared
`GSE104987_seurat_afterAnno.RDS` was temporarily inspected as a Seurat object.
It contained 2,669 cells: 1,597 resistant and 1,072 sensitive, all annotated as
malignant cells; RNA/SCT assays; UMAP/PCA reductions; and 50 PCA components.
This metadata inspection did not execute a RESIST biological analysis step.
The temporary approximately 291 MB input was removed to conserve local storage.

The tutorial carries forward these observed input facts. This documentation update
neither downloaded nor re-inspected that object. No reference DEG counts,
enrichment results, analysis figures, or measured workflow runtimes are supplied.

## Deferred to HPC acceptance testing

- Resolve the specified package environment and record its exact versions.
- Download and inspect the current shared example; preserve its checksum.
- Execute A1, then B1/B4, then C7/C8, inspecting each table and figure.
- Check empty/significance-filtered results and active-assay behavior.
- Enable B2/B3/B5/C9/C10 after their relevant dependencies are prepared.
- Validate A2/A3 and D3 on appropriate multi-cell-type/spatial/APA inputs.
- Review D2 read assignments, recover or validate C2, and define the missing TF
  workflow before claiming those analyses are reproducible.
- Obtain LINCS only when needed, validate its parser/schema, and measure HPC memory.
- Compare numerical outputs with an authoritative reference analysis before
  describing a run as reproducing published results.

Static success establishes a consistent package layout and readable syntax.
It does not establish dependency compatibility, correct barcode mapping,
numerical equivalence, or scientific validity for a particular cohort.
