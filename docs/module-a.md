# Module A — Characterization

Module A visualizes existing annotations and reductions and, when suitable inputs
are available, infers communication or displays spatial distributions. It does
not create a fully annotated dataset from raw reads.

| Step | Input | Output under `<results>/A/` | Availability |
|---|---|---|---|
| A1 | Seurat object; `celltype`/`cell_type`, `condition`, `seurat_clusters`, existing `umap` | `<input_filename>_UMAP.pdf` and `.png` | Default; applicable to the shared GSE104987 example |
| A2 | SCT data; `sample`, `condition`, multiple cell types in each group | Condition-specific circle figures and CellChat comparison/L–R figures | Additional inputs and CellChat required |
| A3 | Spatial Seurat objects with `final_celltype`, `group`, `sample`, images/coordinates | `spatial_plots/*_spatial_cell_distribution_{sensitive,resistant}_top2.pdf/png` | Study-specific spatial template |

```bash
bash run_A.sh --config config/example.yaml          # A1
bash run_A.sh --config config/config.yaml --steps A2 --check
bash run_A.sh --config config/config.yaml --steps A2
bash run_A.sh --config config/config.yaml --steps A3 --check
```

A1's composition panel compares cluster percentages within condition using grouped
bars. It does not write a separate cell-composition table or carry out a
sample-level abundance test. It uses the input UMAP, not a newly fitted embedding.
The reader harmonizes `Malignant cells` to `Tumor cells` in memory.

A2 reads `paths.data`, uses the human CellChat database and SCT assay,
and requires at least two retained cell types per condition (at least ten cells
per type). Its tumor-sender figures also require tumor cells. The shared example
has only malignant cells and no `sample` column, so A2 is not applicable to it.
Do not invent other cell types or replicate IDs to make the checks pass.
The bundled `scripts/lib/netVisual_bubble.R` retains the customized plotting
function supplied with the source code. A2 loads this helper and processes each
eligible dataset within its loop. CellChat output and comparison direction
still require HPC review.

A3 was written for a particular spatial layout. It expects image names derived as
`slice1.<sample>` and selects up to two samples by cell-type diversity per
condition. Adapt and validate that mapping for a new cohort rather than presenting
it as a general spatial preprocessing pipeline. Plots use the shared palette
defined under `scripts/lib/`.

A1/A2 preserve the configured `paired_pre_post_seu` exclusions. Their meaning
must be reviewed for new studies. See [the input contract](input-data.md) and
[scientific limitations](methods-and-limitations.md).
