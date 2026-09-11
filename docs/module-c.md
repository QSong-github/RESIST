# Module C — Regulatory associations

The default C runner covers the RBP branch. miRNA enrichment is selectable; the
variant branch uses its own per-sample files and commands. No TF motif
implementation was supplied.

## RBP and miRNA branches

| Step | Reads | Writes under `<results>/C/` |
|---|---|---|
| C7 | B1 `*_deg.csv`; bundled human RBP target files | `*_deg_rbp_enrichment.csv` |
| C8 | C7 tables | `*_deg_rbp_enrichment_{up,dn}_circle.pdf/png` when gates pass |
| C9 | B1 `*_deg.sig.csv`; miRDB, cached mappings and biomaRt objects | `*_deg.sig_{up,dn}_miRNA.csv` |
| C10 | C9 tables | `*_miRNA_{up,dn}.pdf/png` when eligible |

```bash
bash run_C.sh --config config/example.yaml             # C7,C8; run B first
bash run_C.sh --config config/example.yaml --steps C9,C10 --check
bash run_C.sh --config config/example.yaml --steps C9,C10
```

C7 uses genes returned in the tumor-cell DEG table as its universe and tests target
overlap by direction with a hypergeometric test. The columns `FDR_up` and `FDR_dn`
retain **Holm-adjusted** p-values, as in the source's default `p.adjust()` behavior.
The current circle labels name that method correctly. C8 requires at least five
significant up sets before proceeding, then five down sets for the down plot.
A missing circle figure does not establish that the table failed.

C7 preserves its mouse-dataset exclusions. C9's mappings and species handling
should be checked for each cohort. Enriched target sets suggest candidate
regulators; they do not demonstrate RBP/miRNA expression, activity, or causality.

## Variant branch: supplied code and the missing bridge

| Step | Input → output | Status |
|---|---|---|
| C1 | BAM, barcode list, reference SNP VCF → cellSNP outputs | External cellSNP-lite required |
| C2 | cellSNP VCF → long per-cell AF table | **Not supplied; cannot reproduce this bridge** |
| C3 | Condition-labeled AF tables → `fisher_results.tsv`, `heatmap.pdf/png` | Can run only with appropriately prepared AF tables |
| C4 | Fisher table → significant-locus CSV | Supplied Python script |
| C5 | Locus CSV → annotation CSV | External annotation services; genome-build/locus review required |
| C6 | Tree containing `result/group*/fisher_results.tsv` → `input.csv`, `out.csv` | Legacy batch utility; works from its current directory |

The C3 input table requires `ID`, `Cell`, `s_reads`, `t_reads`, and `AF`.
`ID` has the form `chrom_pos_ref_alt`; `AF` is ALT-supporting reads divided by
total reads. Input filenames must include `sensitive` or `resistant`, because C3
infers the condition from their names. That naming convention must agree with
sample metadata. The missing C2 filtering cannot be inferred from the downstream
column names.

Example interfaces, to use only with your own validated sample paths:

```bash
bash scripts/C/C1_call_variants_cellsnp.sh \
  /path/to/sample.bam /path/to/barcodes.txt /path/to/reference_snps.vcf \
  /path/to/cellsnp_output 8

# C2 must be obtained from the authors or separately implemented and validated.
# Once validated AF tables are available:
Rscript scripts/C/C3_fisher_enrichment_heatmap.R \
  /path/to/condition_labeled_AF_tables /path/to/variant_results
python scripts/C/C4_select_significant_variants.py \
  --input /path/to/variant_results/fisher_results.tsv \
  --output /path/to/variant_results/input.csv
python scripts/C/C5_annotate_variants.py \
  /path/to/variant_results/input.csv /path/to/variant_results/annotated_output.csv
```

C3–C5 are not invoked by `run_C.sh`. Confirm the reference genome, chromosome
convention, allele context, statistical thresholds, and annotation service
versions before reporting variant associations.

## TF motifs

No TF motif-enrichment script or runnable TF workflow was included in v2-2.
This release preserves that limitation instead of adding an empty folder or
claiming a TF result. Establish the intended input regions/gene mapping, motif
database, background, and method with the authors before implementing it.
