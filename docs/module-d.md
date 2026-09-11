# Module D — Immunogenomic analyses

Module D groups three distinct workflows: APA summaries, HLA typing, and candidate
peptide–MHC structural inspection. They require their own data; the shared
GSE104987 Seurat object cannot substitute for raw reads or scUTRquant outputs.

## APA: D1–D3

| Step | Supplied material | Role |
|---|---|---|
| D1 | `scripts/D/D1_detect_kit_version.py` | Read Cell Ranger web-summary metadata to report chemistry |
| D2 | `scripts/D/D2_prepare_cellranger_fastqs.sh` | Legacy FASTQ linking/merging utility |
| External | Templates in `data/templates/` | Example metadata and scUTRquant configuration |
| D3 | `scripts/D/D3_apa_quantification_and_plots.R` | Relative-expression summaries and comparisons from prepared inputs |

[scUTRquant](https://github.com/Mayrlab/scUTRquant) is an external isoform
quantification workflow. Its target annotation and TXS output must match the
assumptions of D3. The small supplied templates are examples, not completed
cohort data.

**D2 requires review before use:** the inherited `link` branch maps SRA read 2 to
R1 and read 3 to R2, while the `merge` branch maps indices 1/2/3 to R1/R2/I1.
These mappings differ. Preserve raw FASTQs and confirm chemistry/read roles from
the actual submission before choosing or adapting a branch. This reorganization
does not guess the correct read assignment.

Configure `config/apa_config.yaml` and `config/apa_samples.csv` as described in
[tutorial Step 7](../TUTORIAL.md#step-7--module-d-with-a-separately-prepared-apa-cohort):

```bash
bash run_D.sh --config config/config.yaml \
  --apa-config config/apa_config.yaml --check
bash run_D.sh --config config/config.yaml \
  --apa-config config/apa_config.yaml
```

The root D runner invokes **D3 only**. Relative paths in both APA configuration
and sample registry resolve against the package root. The annotation input may
be an RDS Seurat object or an RData containing one unambiguous Seurat object.
The final condition mapping comes from the annotation metadata, not directly
from the sample sheet's `group` column.

Outputs under `<results>/D/APA/<cohort_id>/`:

| Product | Interpretation/check |
|---|---|
| `RE_outputs/*_RE_gene_by_cell.rds` | Gene-by-cell relative-expression matrix |
| `RE_outputs/*_RE_mean.csv` | Per-cell mean relative-expression summaries |
| `sce_barcode_to_seurat_mapping.csv` | Trace how barcodes map to annotations |
| `mapping_diagnostics.csv` | Review missing/ambiguous overlap before interpretation |
| `RE_Differential_Analysis_Results.csv` | Eligible differential comparisons |
| `ECDF_by_celltype/`, `violin_by_celltype/` | Distribution figures, when eligible |

D3 loads the configured annotation object through the shared annotation helper
and assigns it to `annot_info` for barcode mapping and diagnostics. Validate
TXS structure, barcode overlap, and mapped conditions on the HPC before
interpreting APA results.

## HLA typing: D4

With a separately prepared OptiType installation and suitable reads:

```bash
bash scripts/D/D4_run_optitype.sh \
  --pipeline /path/to/OptiTypePipeline.py --type rna \
  --out /path/to/hla_results /path/to/hla_reads_1.fastq.gz /path/to/hla_reads_2.fastq.gz
```

The wrapper expects the predicted HLA-A/B/C result table under the requested
output directory. It does not perform upstream read fishing or supply reference
indices. See [OptiType](https://github.com/FRED-2/OptiType).

## Candidate peptide–MHC structures: D5/D6

Candidate generation with [DIPAN](https://github.com/YY-TMU/DIPAN) is external and
has no automated wrapper here. The supplied D5 is a ColabFold notebook export
that imports `google.colab`; do not run it as an ordinary HPC Python program.
The original example sequence is retained in the small templates.

For a separate local ColabFold installation, prepare a **single FASTA record**
whose sequence contains the three chains separated by `:`:

```text
>hla_peptide_ITDVGSGMY
<HLA-heavy-chain>:<beta-2-microglobulin>:<peptide>
```

A ready record derived from the supplied sequence is at
`data/templates/hla_peptide_complex_colabfold.fasta`. After reviewing the sequence
and installed [ColabFold](https://github.com/sokrypton/ColabFold) options, a local
command has this shape:

```bash
colabfold_batch data/templates/hla_peptide_complex_colabfold.fasta \
  /path/to/neoantigen_models
```

Do not assume a three-record FASTA will be interpreted as one three-chain complex.
Model type, MSA, templates, relaxation, software version, and resources should be
chosen and recorded for the actual study. No structures or runtime guarantees
are supplied with this release.

D6 expects the loaded model to contain heavy chain A, B2M chain B, and peptide
chain C. Run it from the directory where you want the image saved:

```bash
pymol -cq /path/to/chosen_model.pdb /path/to/RESIST/scripts/D/D6_visualize_complex.pml
```

It writes `MHC_Peptide_Interaction.png` in the current working directory.
Confirm the chain identities and examine model confidence and interface geometry.
A predicted complex does not establish binding affinity, antigen presentation,
immunogenicity, or therapeutic efficacy.
