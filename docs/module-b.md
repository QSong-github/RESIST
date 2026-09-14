# Module B — Transcriptional programs

B1 supplies the DEG tables consumed by pathway, drug-connectivity, RBP, and miRNA
steps. B4 and B5 independently read the annotated Seurat object; they do not
require an A1 output file.

| Step | Products under `<results>/B/` | Prerequisite |
|---|---|---|
| B1 | `*_deg.csv`, `*_deg.sig.csv`, `*_deg_Tumor_cells_volcano.pdf/png` | Processed Seurat objects; both conditions |
| B2 | `*_deg.sig_{GO,KEGG,Hallmark}_enrichment_{up,dn}.csv` | B1 significant tables; human annotations and database access |
| B3 | `*_{GO,KEGG,HALLMARK}_barchart_sideByside.pdf/png` | B2 result pairs and eligible terms |
| B4 | `<input_filename>_ITH_box.pdf/png` | Tumor cells, both conditions, at least 20 existing PCs |
| B5 | `<input_filename>_EMT_box.pdf/png` | Tumor cells, normalized expression, GSVA and human EMT GMT |
| B6 | `*_drug_enrichment.csv` | B1 significant tables; external LINCS object; high-memory node |
| B7 | `*_drug_enrichment_bar.pdf/png` | B6 results |

```bash
bash scripts/run_B.sh --config config/example.yaml             # B1,B4
bash scripts/run_B.sh --config config/example.yaml --steps B2,B3
bash scripts/run_B.sh --config config/example.yaml --steps B5
# Only after obtaining LINCS and a suitable HPC allocation:
bash scripts/run_B.sh --config config/example.yaml --steps B6,B7
```

B1 tests resistant versus sensitive within each cell type using Seurat's active
assay. It skips groups with fewer than three cells. Significant rows require
`p_val_adj < 0.05` and `abs(avg_log2FC) >= 1`; positive log fold change means higher
expression in resistant cells. The tumor volcano is drawn only when tumor results
are available. The script does not implement a sample-aware pseudobulk model.

B2 focuses on tumor-cell DEGs and currently selects human GO/KEGG/Hallmark
annotations. It does not supply a measured-gene background to the enrichment
functions. B3's plotted results are subject to the inherited filters and may be
empty. Human settings must be reviewed before analyzing mouse cohorts. The
`msigdbr` schema and Seurat/GSVA interfaces can differ across versions; preserve
package versions and confirm compatibility on the HPC.

B4 uses mean pairwise Euclidean distances across all tumor cells in the first 20
PCs. B5 scores the human Hallmark EMT set. Both currently write box plots, not
per-cell score CSVs, and compare conditions using cell-level Wilcoxon tests.
They retain their original scientific definitions. B5's separate
`paired_pre_post_anno` exclusion list remains unchanged.

B6 loads the LINCS object in memory. The source notes describe memory use in the
tens of GB; treat any allocation as a starting point and measure actual usage.
Candidate negative-connectivity signatures are hypotheses for follow-up, not
validated treatment recommendations.

See [references](references.md) before B5/B6 and [methods and limitations](methods-and-limitations.md)
for interpretive boundaries and implementation caveats.
