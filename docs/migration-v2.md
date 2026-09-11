> Historical v2-2 mapping, preserved for provenance. Paths below describe that earlier release, not v3_02. See `migration.tsv` for current paths.

# Migration Map

Every file of the two original bundles and where it went, with the
changes applied. Analysis logic was migrated verbatim; the changes
listed here are confined to configuration, paths and provenance
headers. Nothing that affects a computed result was altered.

## File-by-file

### `Module_A_Characterization/A1_umap_and_composition.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline1_UMAP_replicate.R`

- `ids_timppint` literal vector -> `resist_datasets("paired_pre_post_seu")`
- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- in-loop `source(color.R)` -> `col <- resist_palette("celltype")`
- 1 x `<abs>/results` -> `RESULTS_A`
- 1 x `<abs>/data` -> `PATH_DATA`

### `Module_A_Characterization/A2_cell_cell_communication.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline10_cellchat_batch.R`

- `ids_timppint` literal vector -> `resist_datasets("paired_pre_post_seu")`
- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- in-loop `source(color.R)` -> `col <- resist_palette("celltype")`
- `source(<abs>/netVisual_bubble.R)` -> repository-relative path
- 1 x `<abs>/results` -> `RESULTS_A`
- 1 x `<abs>/data` -> `PATH_SPATIAL`

### `Module_A_Characterization/A3_spatial_distribution.R`

*From* `RESIST_single_cell_pipeline/scripts/final_spatial_distribution_v2.R`

- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- in-loop `source(color.R)` -> `col <- resist_palette("celltype")`
- 1 x `<abs>/results/spatial_plots` -> `file.path(RESULTS_A, "spatial_plots")`
- 1 x `<abs>/data` -> `PATH_SPATIAL`

### `Module_A_Characterization/lib/netVisual_bubble.R`

*From* `RESIST_single_cell_pipeline/scripts/netVisual_bubble.R`

- provenance header added; code unchanged

### `Module_B_Transcriptional/B1_differential_expression.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline2_DEG_replicate.R`

- `ids_timppint` literal vector -> `resist_datasets("paired_pre_post_seu")`
- `ids_harmony` literal vector -> `resist_datasets("harmony_integrated")`
- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- 3 x `<abs>/results` -> `RESULTS_B`
- 1 x `<abs>/data` -> `PATH_DATA`

### `Module_B_Transcriptional/B2_pathway_enrichment.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline3_enrichment_replicate.R`

- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- 2 x `<abs>/results` -> `RESULTS_B`

### `Module_B_Transcriptional/B3_pathway_enrichment_barplot.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline3_enrichment_barchart_side_by_side_replicate.R`

- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- 6 x `<abs>/results` -> `RESULTS_B`

### `Module_B_Transcriptional/B4_intratumor_heterogeneity.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline4_ITH_replicate.R`

- `ids_timppint` literal vector -> `resist_datasets("paired_pre_post_seu")`
- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- 1 x `<abs>/results` -> `RESULTS_B`
- 1 x `<abs>/data` -> `PATH_DATA`

### `Module_B_Transcriptional/B5_emt_score.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline5_EMT_replicate.R`

- `ids_timppint` literal vector -> `resist_datasets("paired_pre_post_anno")`
- top-level `source(color.R)` removed (palettes come from common/R/init.R)
- 1 x `<abs>/ref_data/...` -> `file.path(PATH_REF, ...)`
- 1 x `<abs>/results` -> `RESULTS_B`
- `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt` -> `resist_ref("emt_gmt_human")`
- `list.files('./data',` -> `list.files(PATH_DATA,`

### `Module_B_Transcriptional/B6_drug_enrichment.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline9_drug_enrichment_replicate.R`

- 1 x `<abs>/ref_data/...` -> `file.path(PATH_REF, ...)`
- 2 x `<abs>/results` -> `RESULTS_B`
- `GSE70138_LINCS_Level4.rds` -> `resist_ref("lincs_level4_rds")`

### `Module_B_Transcriptional/B7_drug_enrichment_plot.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline9_drug_enrichment_plot_for_web_replicate.R`

- 2 x `<abs>/results` -> `RESULTS_B`

### `Module_C_Regulatory/mirna/C10_mirna_treemap_plot.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline7_miRNA_plot_replicate.R`

- 2 x `<abs>/results` -> `RESULTS_C`

### `Module_C_Regulatory/mirna/C9_mirna_target_enrichment.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline7_miRNA_replicate.R`

- `deg_dir <- "/blue/qsong1/sen.guo/resist_02/resist_test_share/results"` -> `deg_dir <- RESULTS_B`
- 5 x `<abs>/ref_data/...` -> `file.path(PATH_REF, ...)`
- 1 x `<abs>/results` -> `RESULTS_C`
- `miRDB_v6.0_prediction_result.txt` -> `resist_ref("mirdb")`
- `mart_h.rds` -> `resist_ref("biomart_human")`
- `mart_m.rds` -> `resist_ref("biomart_mouse")`
- `all_mappings_h.csv` -> `resist_ref("gene_map_human")`
- `all_mappings_m.csv` -> `resist_ref("gene_map_mouse")`

### `Module_C_Regulatory/rbp/C7_rbp_target_enrichment.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline8_RBP1_enrichment_replicate.R`

- `list.files('/blue/qsong1/sen.guo/resist_02/resist_test_share/results', pattern = "deg\\.csv$"` -> `list.files(RESULTS_B, pattern = "deg\\.csv$"`
- `file = paste0("/blue/qsong1/sen.guo/resist_02/resist_test_share/results/", deg_file_short, "_rbp_enrichment.csv")` -> `file = file.path(RESULTS_C, paste0(deg_file_short, "_rbp_enrichment.csv"))`
- 1 x `<abs>/ref_data/...` -> `file.path(PATH_REF, ...)`
- 1 x `<abs>/results` -> `RESULTS_C`
- `All_RBP_TargetGene` -> `resist_ref("rbp_target_dir")`

### `Module_C_Regulatory/rbp/C8_rbp_bubble_plot.R`

*From* `RESIST_single_cell_pipeline/scripts/pipeline8_RBP_plot1_bubble_replicate.R`

- 2 x `<abs>/results` -> `RESULTS_C`

### `Module_C_Regulatory/variants/C1_call_variants_cellsnp.sh`

*From* `README.md (PART I Step 1, documented command)`

- command from the original README captured as an executable, argument-checked script

### `Module_C_Regulatory/variants/C3_fisher_enrichment_heatmap.R`

*From* `RESIST-main/Variant_detection/scMutation_fisher_heatmap.R`

- commented-out `library()` calls restored and completed (tidyr/purrr/pheatmap were used but never declared)
- hard-coded `data_dir` -> command-line argument; results now go to `results/Module_C_Regulatory/variants/` instead of the input tree
- three output paths redirected from `data_dir` to `out_dir`

### `Module_C_Regulatory/variants/C4_select_significant_variants.py`

*From* `RESIST-main/Variant_detection/format.py`

- provenance header added; code unchanged

### `Module_C_Regulatory/variants/C5_annotate_variants.py`

*From* `RESIST-main/Variant_detection/mutation_annotation.py`

- provenance header added; code unchanged

### `Module_C_Regulatory/variants/C6_batch_annotate.py`

*From* `RESIST-main/Variant_detection/annotation.py`

- provenance header added; code unchanged

### `Module_D_Immunogenomic/apa/D1_detect_kit_version.py`

*From* `RESIST-main/APA_analysis/find_kit_version.py`

- provenance header added; code unchanged

### `Module_D_Immunogenomic/apa/D2_prepare_cellranger_fastqs.sh`

*From* `RESIST-main/APA_analysis/Untitled-1.sh`

- scratch file rewritten as a documented, argument-driven script
- hard-coded SRR accessions and absolute FASTQ paths replaced by arguments
- REMOVED: an unrelated RNAmigos inference command (out of scope for RESIST)
- REMOVED: four near-identical `join` commands that merged Cell Ranger clustering and t-SNE tables into annots.csv, each overwriting the last; see workflows/utils/merge_cellranger_annotations.sh for a corrected form
- REMOVED: a free-text note to a code reviewer left at the end of the file

### `Module_D_Immunogenomic/apa/D3_apa_quantification_and_plots.R`

*From* `RESIST-main/APA_analysis/APA_processing_graph.R`

- hard-coded GTF path, annotation object, output directory and the inline two-sample `tribble()` replaced by `config/apa_config.yaml` plus a sample-sheet CSV
- output now written under `results/Module_D_Immunogenomic/APA/<cohort_id>/` instead of a relative directory in the working directory

### `Module_D_Immunogenomic/apa/metadata/kit.csv`

*From* `RESIST-main/APA_analysis/kit.csv`

- moved (example output of D1; absolute cluster paths replaced by placeholders)

### `Module_D_Immunogenomic/apa/metadata/list.txt`

*From* `RESIST-main/APA_analysis/list.txt`

- moved (example input to D1)

### `Module_D_Immunogenomic/apa/metadata/raw_data_info.txt`

*From* `RESIST-main/APA_analysis/raw_data_info.txt`

- moved (example sample-to-group mapping)

### `Module_D_Immunogenomic/apa/scutrquant/config.yaml`

*From* `RESIST-main/APA_analysis/config.yaml`

- moved; input to the external scUTRquant workflow, not read by RESIST
- `sample_file` and `cell_annots` absolute cluster paths replaced by placeholders; all analysis settings unchanged
- provenance header added

### `Module_D_Immunogenomic/apa/scutrquant/sample_sheet.csv`

*From* `RESIST-main/APA_analysis/sample_sheet.csv`

- moved; sample ids and file types unchanged, absolute BAM paths replaced by placeholders

### `Module_D_Immunogenomic/hla/D4_run_optitype.sh`

*From* `README.md (PART II Step 1, documented command)`

- command from the original README captured as an executable, argument-checked script

### `Module_D_Immunogenomic/neoantigen/D5_fold_hla_peptide_complex.py`

*From* `RESIST-main/Neoantigen_visualization/alphafold2.py`

- provenance header added; notebook body unchanged
- `query_sequence` and `jobname` form-field defaults set to the worked HLA:B2M:peptide example (neoantigen/example/); no executable statement altered
- header corrected: the query is three chains (heavy:B2M:peptide), not two

### `Module_D_Immunogenomic/neoantigen/D6_visualize_complex.pml`

*From* `RESIST-main/Neoantigen_visualization/pymol_command`

- renamed to the conventional `.pml` extension; provenance header added; commands unchanged
- header chain assignment corrected to A = heavy chain, B = B2M, C = peptide, which is what the selections in the body already assume

### `common/R/palette.R`

*From* `RESIST_single_cell_pipeline/scripts/color.R`

- palette vectors migrated verbatim; library() calls moved to common/R/packages.R; accessor `resist_palette()` added

## Files added in the reorganization

| Path | Why |
|------|-----|
| `config/config.yaml` | Single source of truth for data, reference and results paths, dataset hold-out lists, and figure defaults. |
| `config/apa_config.yaml`, `config/apa_samples.csv` | Externalise the cohort settings that Module D's APA script previously carried inline. |
| `common/R/init.R` | Loaded by the bootstrap header of every R entry point; sources the four files below. |
| `common/R/config.R` | Reads the YAML, resolves paths against the repository root, exposes `resist_ref()`, `resist_datasets()`, `resist_results_dir()`. |
| `common/R/packages.R` | The library calls that used to sit at the top of `color.R` and each pipeline script. |
| `common/R/palette.R` | The palettes from `color.R`, plus `resist_palette()`. |
| `common/R/io.R` | `resist_save_figure()` (paired PDF/PNG), `resist_input_objects()`, `resist_log()`. |
| `workflows/run_module_{a,b,c,d}.sh`, `workflows/run_all.sh` | Ordered, path-independent runners; replace the absolute `Rscript` list in `run.sh`. |
| `workflows/README.md` | Coverage table: which steps each runner covers, and why the per-sample and external-tool steps have none. |
| `workflows/slurm/submit_all.sbatch` | The SLURM header from `run.sh`, without the hard-coded project directory. |
| `workflows/utils/merge_cellranger_annotations.sh` | The `join`/`sed` incantation from the scratch shell file, parameterised and no longer self-overwriting. |
| `install/environment.yml`, `install/install_r_packages.R`, `install/requirements.txt` | Declared dependencies; previously only a comment in `run.sh` saying that some R packages may need installing. |
| `.gitignore` | Keeps data, results and the three large reference downloads out of version control. |
| `Module_D_Immunogenomic/neoantigen/example/` | Worked pMHC class I query (heavy chain : B2M : `ITDVGSGMY`) in FASTA and `:`-joined form, so `D5` and `D6` can be run end to end without supplying data. |
| Module and directory `README.md` files | Per-module documentation, including the two gaps recorded below. |

## Deliberate removals

| Removed | Reason |
|---------|--------|
| RNAmigos inference command in `APA_analysis/Untitled-1.sh` | Belongs to an unrelated project; nothing in RESIST calls it. |
| Four repeated `join` commands in the same file | All four wrote to the same `annots.csv`, so only the last survived. Replaced by `workflows/utils/merge_cellranger_annotations.sh`, which names the output after the accession. |
| A free-text note to a code reviewer at the end of the same file | Not code. The defect it reports is recorded in `docs/repository-guide.md` under Known issues. |
| `ref_data/readme_refdata.txt`, `data/example_data_download_link.txt` | Superseded by `ref_data/README.md` and `data/README.md`. |

## Not carried over, because it was never there

| Expected | Status |
|----------|--------|
| `trans.py` (VCF to per-cell AF table, step C2) | Referenced by the original README and by the header of the Fisher script, but absent from the handover. Not reconstructed; the contract it must satisfy is documented in `Module_C_Regulatory/README.md`. |
| TF motif enrichment (Figure 1, Module C panel 3) | No implementation was supplied. `Module_C_Regulatory/tf_motif/README.md` records the intended interface. |
