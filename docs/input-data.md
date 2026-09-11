# Input data contract

RESIST begins with processed objects and branch-specific prepared inputs. It does
not replace upstream sequencing QC, normalization, sample integration, or annotation.

## Single-cell objects

Set `paths.data` in the selected YAML. The directory is scanned non-recursively
for `.rds`/`.RDS` files. Use one annotated Seurat object per dataset and preserve
its accession-based filename. Keep reference RDS files and intermediate objects
elsewhere. Dataset IDs are inferred from filename conventions in the inherited
scripts; do not rename files casually, and inspect output basenames for collisions.

| Field | Contract |
|---|---|
| Object class | `Seurat` |
| Cell type | `celltype` or `cell_type`, nonmissing; contradictory duplicate columns are rejected |
| Condition | Exact `sensitive` / `resistant` labels; derive them from study metadata |
| A1 clusters | `seurat_clusters` |
| A1 embedding | Existing `umap` reduction |
| B4 embedding | Existing `pca`, at least 20 components |
| Expression | Correctly prepared assay appropriate to the inherited method; B1 uses the active assay |
| A2 additional metadata | Verified `sample` labels, multiple cell types, normalized SCT data |

`Malignant cells` is standardized to `Tumor cells` in memory for A1/A2/B1/B4/B5.
No cell-type predictions or experimental labels are created. A1 can display one
condition, but differential comparisons require both. B1 requires at least three
cells per group within an eligible cell type.

The exact shared object is documented in [the tutorial](../TUTORIAL.md). Its
observed dimensions and metadata are a useful input checkpoint, not a substitute
for checking study design or the original preprocessing provenance.

## Other branches

| Input root/configuration | Needed for |
|---|---|
| `paths.spatial_data` | A3 spatial objects with the study-specific metadata/image mapping |
| B1 tables in the same `<results>/B/` | B2/B6 and C7/C9 |
| `paths.ref_data` plus `reference_files` | Shared and separately obtained references |
| C1 command-line arguments | BAM, barcodes, reference SNP VCF |
| C3 command-line arguments | Validated per-cell AF tables from the missing/external C2 conversion |
| `config/apa_config.yaml` and its sample sheet | D3 annotation Seurat, GTF, filtered 10x matrices, TXS RDS |
| D4 arguments | OptiType and appropriate FASTQ files |
| D5/D6 separate environment | Candidate sequences, ColabFold outputs, chosen model |

Absolute paths allow large files to remain in shared storage. A relative
`--config` argument is resolved from the shell's current directory; relative paths
inside the configuration are resolved from the RESIST package root.
