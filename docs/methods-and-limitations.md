# Scientific scope and limitations

RESIST combines analyses with distinct data requirements and statistical
assumptions. Choose the branches appropriate to the study design and validate
their numerical behavior on the HPC. The static package checks do not
establish reproduction of published results.

## Interpretation by analysis

| Result | What it supports | What it does not establish |
|---|---|---|
| UMAP/cluster percentages | Description of processed cells and annotations | A de novo annotation or replicated differential abundance |
| Cell-level DEGs | Within-cell-type expression association between labels | A sample-level causal treatment effect |
| PCA-distance ITH / EMT scores | Exploratory differences under inherited score definitions | Independent-patient inference or universal biological thresholds |
| Pathway/RBP/miRNA enrichment | Candidate programs or target-set associations | Direct regulator activity or causality |
| Variant enrichment/annotation | Candidate locus associations requiring validation | Somatic status, pathogenicity, or causal resistance by itself |
| LINCS connectivity | Candidate perturbation signatures for follow-up | Clinical drug effectiveness |
| APA/HLA/pMHC outputs | Prepared-input molecular summaries or model hypotheses | Validated antigen presentation, immunogenicity, or therapy response |

## Statistical and implementation details to retain

- **Contrast:** B1 compares resistant to sensitive. Seurat `p_val_adj` is an
  adjusted p-value, not a column that should automatically be renamed FDR.
  Its [reference](https://satijalab.org/seurat/reference/findmarkers) describes
  the default Bonferroni adjustment and gene filtering. Record your actual
  Seurat version and active assay.
- **Experimental unit:** the original B1/B4/B5 comparisons operate on cells.
  Cells from the same sample are not independent biological replicates.
  Where the design contains replicates, choose an appropriate sample-aware
  analysis; Seurat's [differential-expression vignette](https://satijalab.org/seurat/articles/de_vignette)
  provides a pseudobulk example. No pseudobulk model is added in this release.
- **Exclusions:** the distinct `paired_pre_post_seu`, `paired_pre_post_anno`, and
  `harmony_integrated` filename lists are preserved. Their original comments
  are not proof that every listed dataset must be excluded in a new study.
  Reconcile cohort design and naming before reuse.
- **Species:** B2 and B5 retain active human branches; C7 excludes named mouse
  datasets. A2 selects the human CellChat database. A supplied mouse reference
  does not mean the relevant script automatically switches species.
- **Gene universe:** B2 uses enrichment-function defaults; C7 uses the genes
  returned in the tumor DEG table. These choices affect enrichment and should
  be reviewed rather than silently described as all measured genes.
- **RBP adjustment:** C7's legacy `FDR_up`/`FDR_dn` columns contain Holm-adjusted
  p-values. C7 explicitly uses `method="holm"`, preserving the source method.
  C8 legends identify this adjustment as Holm.
- **Conditional figures:** C8 and other plotting scripts filter their inputs.
  A lack of a figure can be a valid no-output step; it is not evidence of a
  negative biological result without inspecting the underlying table.

## Remaining blockers and assumptions

1. **C2 is absent.** The cellSNP-to-AF conversion cannot be reproduced from the
   supplied code. No replacement filtering is invented.
2. **TF motif code is absent.** There is no runnable TF motif result.
3. **A3 is a cohort-specific template.** Metadata and image-name assumptions
   remain; they require adaptation for other spatial objects.
4. **D2 read mappings differ between branches.** Verify chemistry and FASTQ
   roles before use; the release preserves the original mapping choices.
5. **D3 needs cohort validation.** Annotation loading and package-relative paths
   are implemented, but TXS schema, barcode matching, group assignment, and
   statistical interpretation require HPC testing. RData inputs must
   contain an unambiguous Seurat annotation object.
6. **D5 is Colab-specific.** Use the external ColabFold CLI for an HPC workflow
   and record its options/environment; no equivalence of model outputs is claimed.
7. **Dependencies are not pinned completely.** Older Seurat/GSVA/msigdbr/CellChat
   interfaces and cached objects may need version-specific adjustments.
8. **All branches have not been executed.** Static parsing and path/packaging
   checks are documented in the release validation record; they are not an
   end-to-end analysis benchmark.

The resource counts in `CITATION.cff` match the supplied resource description
and workflow figure. No fresh census of the RESIST database was performed here. Complete publication metadata,
authorship, dataset provenance, and redistribution/license information with the
project authors before public release; do not invent a DOI or software license.
