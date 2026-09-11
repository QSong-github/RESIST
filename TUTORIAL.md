# Tutorial: from the GSE104987 example to RESIST results

This tutorial explains the commands to run **later on the HPC**. It begins with
the shared annotated Seurat object described below. It does not require a LINCS matrix,
miRDB download, BAM, FASTQ, or GPU for the core A–C walkthrough.

**Execution status:** the example object's structure was inspected during preparation;
the analysis commands below await HPC execution and validation. Expected
output names come from the source code, not from an invented successful run.
HPC execution and numerical validation remain to be performed.

## What this example can teach

| Item | Observed in the shared example |
|---|---|
| Filename | `GSE104987_seurat_afterAnno.RDS` |
| Download size | Approximately 291.3 MB, as shown by the shared folder |
| Object | Annotated Seurat object; RNA and SCT assays; SCT active |
| Cells | 2,669 total: 1,597 `resistant`, 1,072 `sensitive` |
| Cell annotation | All 2,669 labeled `Malignant cells` |
| Existing reductions | UMAP and PCA; 50 PCA components |
| Source labels | `orig.ident`: `C70R` and `MCF7` |
| Not present | A dedicated `sample` column; multiple cell types; spatial images/APA inputs |

These are observations of the supplied processed object, not a new reconstruction
of the original experiment. Preserve the filename and verify the source study
before assigning treatment, replicate, or patient identities. In particular,
`orig.ident` is not automatically an independent biological-replicate identifier.

**Core walkthrough:** A — UMAP and cluster composition; B — differential
expression and heterogeneity; C — RBP target enrichment and plots. C uses
B's full DEG table. **Extensions:** pathway enrichment, EMT scoring, and miRNA
target enrichment with their additional dependencies. **Separate data:**
communication, spatial and variant analyses, and Module D. Technical script
identifiers are provided alongside the detailed methods below.

## Step 1 — Prepare the HPC environment

Transfer the RESIST package ZIP to your HPC project space, extract it, and enter
the extracted directory. The path below is a placeholder for that directory.
Follow your institution's rules for environment installation and
compute allocation; perform analysis on a compute node.

```bash
cd /your/HPC/project/RESIST
conda env create -f config/setup/environment.yml
conda activate resist
Rscript config/setup/install_r_packages.R
python3 --version
Rscript --version
```

Replace `/your/HPC/project/RESIST` once with your actual location. The supplied
environment is an installation specification, not a fully resolved lockfile.
The R installer adds missing CRAN/Bioconductor packages and CellChat. It downloads
software dependencies, not the large reference datasets. Resolve any installation
errors with the HPC's supported compiler/R environment before running analyses.

Python 3 is used by the module launcher. The dependencies in
`config/setup/requirements.txt` are needed only for the separate Python analysis
tools; install them later if using C4–C6 or D1. D5 has its own ColabFold environment.

## Step 2 — Download the example directly onto the HPC

The example is in the [RESIST shared example folder](https://drive.google.com/drive/folders/1ve9xCxtnPBiF7JgU9kDDzOgzA5dCQl2F).
The shared file is
[GSE104987_seurat_afterAnno.RDS](https://drive.google.com/file/d/1orpBafB6Ii-xW8dMOxtdCI9P2c-OtT-I/view).

On an HPC node that permits downloads, use the public-file downloader
[gdown](https://github.com/wkentaro/gdown):

```bash
python -m pip install gdown
mkdir -p data/example
gdown 'https://drive.google.com/uc?id=1orpBafB6Ii-xW8dMOxtdCI9P2c-OtT-I' \
  -O data/example/GSE104987_seurat_afterAnno.RDS
ls -lh data/example/GSE104987_seurat_afterAnno.RDS
sha256sum data/example/GSE104987_seurat_afterAnno.RDS
```

Record the SHA-256 with your run. No publisher-supplied immutable checksum was
provided with the example, so this records the file you downloaded rather than
proving identity with every earlier copy. If Drive reports a quota or permission
error, use the shared folder or ask the data owner to restore access. Do not save
a login/HTML page with an `.RDS` extension and treat it as an object.

Only place this Seurat input in `data/example/`. Reference RDS files such as
`mart_h.rds` belong in `data/reference/`; the analysis discovers all `.rds`/`.RDS`
files immediately inside the configured input directory.

## Step 3 — Unpack the bundled references and inspect the input

The archive contains the supplied RESIST reference files. It is about
16.8 MB compressed and 75.3 MB expanded. Extract it **on the HPC**:

```bash
tar -xzf data/reference_bundle.tar.gz -C data
```

For a shared reference location, see [reference storage](docs/references.md).
Do not download LINCS or miRDB for the core walkthrough.

`config/example.yaml` already selects these paths:

```yaml
paths:
  data: data/example
  spatial_data: data/spatial
  ref_data: data/reference
  results: data/results/example
```

The exact keys are `spatial_data` and `ref_data`. Keep the other settings in the
file, including reference keys and configured dataset exclusions. To inspect the
object yourself:

```bash
Rscript - <<'RS'
suppressPackageStartupMessages(library(Seurat))
x <- readRDS("data/example/GSE104987_seurat_afterAnno.RDS")
print(x)
print(table(x$condition, useNA = "ifany"))
print(table(x$celltype, useNA = "ifany"))
print(Reductions(x))
print(dim(Embeddings(x, "pca")))
print(DefaultAssay(x))
RS
```

**Checkpoint:** compare the metadata with the table above. Stop and investigate
if conditions are missing, the file is not a Seurat object, or the expected
embeddings are absent. This package starts from processed, annotated objects; it
does not perform raw-read QC, normalization, clustering, or cell-type annotation
as part of A1. Do not regenerate those steps merely to match a screenshot.

## Step 4 — A: inspect the embedding and cluster composition

```bash
bash run_A.sh --config config/example.yaml --check
bash run_A.sh --config config/example.yaml
```

`--check` reads prerequisites and reports the object and package versions without
running A1. The second command executes A1 and writes:

```text
data/results/example/A/
  GSE104987_seurat_afterAnno.RDS_UMAP.pdf
  GSE104987_seurat_afterAnno.RDS_UMAP.png
```

The retained `.RDS` in the figure basename is intentional compatibility with the
source script. This is one composite figure with four panels: existing UMAP colored
by cluster, condition, and cell type, plus grouped bars showing cluster percentages
within each condition. The bars describe **cluster composition**, not independent
patient-level estimates or a cell-type abundance test.

The reader standardizes `Malignant cells` to `Tumor cells` in memory and accepts
either `celltype` or `cell_type`; it does not change the saved input object. Because
this example has one annotated cell type, a single tumor label is expected.

**Inspect:** the PDF opens, both conditions appear, and cluster labels agree with
the input. Do not infer cell–cell communication or spatial localization from this
plot. A2 and A3 are not part of the default run.

## Step 5 — B: obtain differential expression and heterogeneity

```bash
bash run_B.sh --config config/example.yaml --check
bash run_B.sh --config config/example.yaml
```

This runs B1 followed by B4. Expected products are:

| Filename under `data/results/example/B/` | Meaning |
|---|---|
| `GSE104987_deg.csv` | Full table returned by the per-cell-type comparison |
| `GSE104987_deg.sig.csv` | Rows with `p_val_adj < 0.05` and `abs(avg_log2FC) >= 1` |
| `GSE104987_deg_Tumor_cells_volcano.pdf` / `.png` | Tumor-cell effect sizes and adjusted p-values |
| `GSE104987_seurat_afterAnno.RDS_ITH_box.pdf` / `.png` | Distribution of the PCA-distance heterogeneity score |

B1 uses `FindMarkers(ident.1 = "resistant", ident.2 = "sensitive")` with the
object's active assay. Positive `avg_log2FC` indicates higher expression in
resistant cells; negative values indicate higher expression in sensitive cells.
`pct.1` and `pct.2` refer to those groups in that order. A cell type with fewer
than three cells in either group is skipped. The full table is the tested output,
not a guarantee that every gene in the assay was tested. Seurat's adjustment and
filtering are described in its [FindMarkers reference](https://satijalab.org/seurat/reference/findmarkers).

Inspect the CSVs after the HPC run:

```bash
Rscript - <<'RS'
d <- read.csv("data/results/example/B/GSE104987_deg.csv")
s <- read.csv("data/results/example/B/GSE104987_deg.sig.csv")
stopifnot(all(c("gene", "cell_type", "avg_log2FC", "p_val_adj") %in% names(d)))
cat("Tested rows:", nrow(d), " Significant rows:", nrow(s), "\n")
print(head(d[order(d$p_val_adj), c("gene", "avg_log2FC", "p_val_adj")], 10))
RS
```

No DEG count or list of top genes is prescribed here: those are results of your
future run. An empty significant table can be a valid outcome.

B4 takes the first 20 existing PCs, computes pairwise Euclidean distances among
all tumor cells, and assigns each cell its mean distance to the others. Its plot
compares the score distributions by condition. It does **not** currently export
per-cell ITH CSVs. This calculation can require substantial RAM for large cohorts
because the pairwise distance matrix grows with the square of the cell count.

**Interpretation:** cell-level contrasts do not by themselves establish a
replicated treatment effect. For new cohorts, use a design appropriate to their
biological replicates; see [methods and limitations](docs/methods-and-limitations.md).

## Step 6 — C: examine RBP target-set enrichment

C uses the B1 table in the **same output profile**. Run B first; changing the
results path between B and C breaks that dependency.

```bash
bash run_C.sh --config config/example.yaml --check
bash run_C.sh --config config/example.yaml
```

Expected main table:

```text
data/results/example/C/GSE104987_deg_rbp_enrichment.csv
```

C7 compares up- and down-regulated tumor DEGs with each bundled human RBP target
set. Inspect `RBP`, `Targets`, `Up_DEG`, `Up_overlap`, `Up_pval`, `Down_DEG`,
`Down_overlap`, `Down_pval`, `FDR_up`, and `FDR_dn`.

**Naming caveat:** the columns `FDR_up` and `FDR_dn` contain
**Holm-adjusted p-values**. C7 explicitly uses `p.adjust(method="holm")`,
preserving the source method and column names. These values are not BH FDR.

C8 can additionally write:

```text
GSE104987_deg_rbp_enrichment_up_circle.pdf / .png
GSE104987_deg_rbp_enrichment_dn_circle.pdf / .png
```

The C8 plotting gate requires at least five significant up-regulated RBP
sets before producing the up plot; it then requires at least five significant
down-regulated sets for the down plot. Thus a table may exist without either
plot, and a down plot can be suppressed by the earlier up-set gate. The launcher
reports `no_new_outputs` for a plotting step that writes nothing. Read the table
and gate before treating an absent figure as an execution failure.

**Interpretation:** enrichment prioritizes regulator hypotheses. It does not
measure RBP activity, prove direct regulation, or establish a causal mechanism.

## Step 7 — Module D with a separately prepared APA cohort

The shared GSE104987 `.RDS` is insufficient for D. Keep this part separate from
the A–C example rather than fabricating FASTQ, APA, or HLA results.

1. Prepare the cohort's Cell Ranger filtered 10x matrices and obtain scUTRquant
   `.txs.Rds` outputs using a compatible target annotation.
2. Identify the exact GTF used for the poly(A)-site mapping and a Seurat annotation
   object whose cell barcodes can be matched to that cohort.
3. Edit `config/apa_config.yaml`: set `cohort_id`, `gff_path`,
   `annotation_seurat_path`, and `sample_sheet` to real HPC paths.
4. Replace every placeholder row in `config/apa_samples.csv`. Required columns are
   `sample_id,group,tenx_dir,txs_rds`; use one row per actual sample.
5. Check before executing, then run D3:

```bash
bash run_D.sh --config config/config.yaml \
  --apa-config config/apa_config.yaml --check
bash run_D.sh --config config/config.yaml \
  --apa-config config/apa_config.yaml
```

Relative APA paths resolve against the package root. The sample sheet's `group`
column is not the source of the final sensitive/resistant mapping in
D3; the annotation object's mapped conditions determine that assignment.
Inspect mapping failures and `Unknown` assignments before interpreting contrasts.

Expected products under `data/results/D/APA/<cohort_id>/` include
`RE_outputs/*_RE_gene_by_cell.rds`, `RE_outputs/*_RE_mean.csv`,
`sce_barcode_to_seurat_mapping.csv`, `mapping_diagnostics.csv`, and
`RE_Differential_Analysis_Results.csv`, with eligible ECDF/violin figures.
The path check does not validate all TXS internals or barcode mapping. That is
part of your HPC acceptance testing. See [Module D](docs/module-d.md) for the
separate HLA and structural workflows.

## Step 8 — Add only the analyses needed for your question

| Extension | Command after preparing its inputs | Additional requirement |
|---|---|---|
| Pathway tables and plots | `bash run_B.sh --config config/example.yaml --steps B2,B3` | B1 outputs; human annotations; database/network access |
| EMT plot | `bash run_B.sh --config config/example.yaml --steps B5` | GSVA/GSEABase and the bundled human EMT gene set |
| miRNA tables and treemaps | `bash run_C.sh --config config/example.yaml --steps C9,C10` | B1 outputs; separately downloaded miRDB reference |
| LINCS connectivity | `bash run_B.sh --config config/example.yaml --steps B6,B7` | Large LINCS RDS and a high-memory HPC allocation |

Append `--check` first for any extension. See [references](docs/references.md)
for what is bundled and what must be obtained separately. Variant, communication,
spatial, and immunogenomic workflows need additional data; merely selecting their
step IDs does not create those inputs.

## Step 9 — Review and preserve the run

Each executed module creates a unique directory at
`<results>/logs/<UTC_timestamp>-<module>-<id>/`. It contains `preflight.txt`, a
configuration snapshot, per-step `.log` files, and `run.json`. The JSON records
completed/failed/no-output steps and the paths, sizes, and SHA-256 hashes of files
written in that run. Previous result files stay in the output directory; use a
new results path for a clean comparison between configurations.

For the first HPC run, confirm that A produces its composite figure; B produces
readable DEG tables and interpretable volcano/ITH figures; and C produces the RBP
table, with plotting gates recorded when figures are absent. Confirm dependencies
and sample metadata before enabling extensions. Save environment versions and
input/reference checksums alongside the logs. Use [the HPC guide](docs/hpc.md)
to submit each module as a job.

## Troubleshooting

| Symptom | Check and next action |
|---|---|
| No `.rds` inputs found | Check `paths.data`, the directory level, and filename extension. |
| R reports an invalid input format | Inspect the file size/type; a Drive HTML response is not an RDS. |
| Conflicting cell-type columns | Reconcile `celltype` and `cell_type` with the source annotation; do not keep contradictory labels. |
| Missing UMAP/PCA | Confirm the processed example was downloaded; raw count files are different inputs. |
| SCT multiple-model or joined-layer error | Review the Seurat version and assay preparation with the analysis owner; preserve the scientific contrast. |
| Missing reference | Extract the supplied archive on the HPC or obtain the specific external reference. |
| C cannot find B output | Use the same configuration/results root for B and C. |
| A2 rejects GSE104987 | Expected: it contains only one cell type and lacks a `sample` column. |
| D3 rejects `/path/to/…` | Replace the APA placeholders with actual cohort inputs. |
| A plot is absent | Inspect the run status and method-specific eligibility/significance gate. |
| Job runs out of memory | Inspect scheduler usage; B4 is quadratic in cell count and B6 loads a large LINCS object. |

The future numerical results should be judged against the source methods and
cohort design. A completed shell command alone is not scientific validation.
