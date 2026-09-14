# Output catalog

Paths below are relative to the configured results root. For the worked example, that root is `data/results/example/`. Outputs are expected products from the source code; no numerical results are asserted. Figures may be conditional on eligible data.

## Module A

| Step | Expected products | Required inputs/limitations |
|---|---|---|
| A1 | `A/` → `<input_filename>_UMAP.pdf/png` | Seurat object with condition, cell type, clusters and UMAP |
| A2 | `A/` → `*circle_sensitive_resistant.pdf/png`; `*_cellchat_*.pdf/png` | SCT assay, sample labels, multiple cell types in both conditions |
| A3 | `A/` → `spatial_plots/*.pdf/png` | Study-specific spatial objects with image/coordinate data |

## Module B

| Step | Expected products | Required inputs/limitations |
|---|---|---|
| B1 | `B/` → `<dataset>_deg.csv`; `<dataset>_deg.sig.csv`; `<dataset>_deg_Tumor_cells_volcano.pdf/png` | Annotated Seurat object; >=3 cells per condition within a cell type |
| B2 | `B/` → `*_deg.sig_{GO,KEGG,Hallmark}_enrichment_{up,dn}.csv` | B1 significant tables; human annotations; internet/database access |
| B3 | `B/` → `*_{GO,KEGG,HALLMARK}_barchart_sideByside.pdf/png` | B2 tables; nonempty eligible pathways |
| B4 | `B/` → `<input_filename>_ITH_box.pdf/png` | Tumor cells in both conditions; >=20 existing PCA components |
| B5 | `B/` → `<input_filename>_EMT_box.pdf/png` | Tumor cells; normalized expression; human EMT GMT |
| B6 | `B/` → `*_drug_enrichment.csv` | B1 significant tables; large LINCS RDS; high-memory node |
| B7 | `B/` → `*_drug_enrichment_bar.pdf/png` | B6 tables |

## Module C

| Step | Expected products | Required inputs/limitations |
|---|---|---|
| C7 | `C/` → `*_deg_rbp_enrichment.csv` | B1 full DEG tables; bundled human RBP targets |
| C8 | `C/` → `*_deg_rbp_enrichment_{up,dn}_circle.pdf/png` | C7 tables; significance/count filters may produce no plot |
| C9 | `C/` → `*_deg.sig_{up,dn}_miRNA.csv` | B1 significant tables; miRDB and mapping references |
| C10 | `C/` → `*_miRNA_{up,dn}.pdf/png` | C9 tables; nonempty significant miRNA results |

## Module D

| Step | Expected products | Required inputs/limitations |
|---|---|---|
| D3 | `D/` → `APA/<cohort>/RE_outputs/*`; `mapping_diagnostics.csv`; `RE_Differential_Analysis_Results.csv`; `ECDF/violin plots` | Actual scUTRquant TXS, 10x matrices, annotation Seurat, GTF and APA sample sheet |

## Separately invoked tools

| Step | Output location/product | Prerequisite |
|---|---|---|
| C1 | Requested output directory; cellSNP VCF and supporting files | BAM, barcodes, reference VCF, cellSNP-lite |
| C2 | No script/output supplied | Missing conversion must be recovered or independently validated |
| C3 | Requested directory or `<results>/C/variants/`; `fisher_results.tsv`, `heatmap.pdf/png` | Validated condition-labeled AF tables |
| C4/C5 | Requested significant-locus and annotation CSV paths | Fisher table; annotation services and correct genome context |
| C6 | Existing `result/group*` directories; `input.csv`, `out.csv` | Legacy directory convention; scans its current directory |
| D1 | Requested CSV, normally `kit.csv` | Cell Ranger web-summary HTML paths |
| D2 | Requested FASTQ output directory | Verified SRA read roles; link/merge mapping differs in inherited code |
| D4 | OptiType output directory; result TSV | Appropriate FASTQ inputs and OptiType |
| D5 | Colab/ColabFold job output; candidate structural models | External environment, sequences, model settings/resources |
| D6 | Current directory; `MHC_Peptide_Interaction.png` | Chosen model with verified chains A/B/C |
| TF motifs | None | No implementation supplied |

## Run record

`<results>/logs/<run_id>/` contains `config.yaml`, `software.json`, one `.log` per executed step, and `run.json`; D also copies its APA configuration. The JSON records script/config/output hashes and exact files written during each step. The software record contains the R session, library paths, and installed-package versions. Input and reference checksums must be recorded separately as described in the tutorial; the runner does not read large inputs merely to hash them.

## Filenames and absent outputs

`<input_filename>` includes the original `.RDS` extension for A1/B4/B5. B1 removes the `_seu…` suffix when deriving its dataset ID. These inherited naming conventions are retained, so the example has both `GSE104987_…` table names and `GSE104987_seurat_afterAnno.RDS_…` figure names. B4/B5 do not write per-cell score CSVs.

C8 can write no figures because of its significance/count gates; consult the tutorial. Similar filtering and missing features can affect other plots. A launcher record of `no_new_outputs` identifies that no files were newly written; it does not by itself prove a biologically null result.

See the [tutorial](../TUTORIAL.md) for exact example filenames and the [module guides](../README.md#analysis-modules) for branch-specific detail.

The R installer writes `installed-packages.csv`, `sessionInfo.txt`, and
`software-checks.txt` under `data/results/setup/<UTC timestamp>-<process id>/`.
A failed check produces a report with failed status, not a success certificate.
