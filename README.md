# RESIST – Neoantigen and Mutation Analysis Workflow

This repository contains a complete workflow for:

- Single-cell mutation enrichment analysis  
- HLA typing  
- Neoantigen prediction  
- HLA–peptide structural modeling  
- PyMOL visualization  

---

# PART I – Single-Cell Mutation Enrichment Pipeline

## Workflow Overview

```
BAM + barcode list + reference SNPs
→ cellSNP-lite
→ cellSNP.cells.vcf.gz
→ trans.py
→ cell_AF_table.tsv
→ R workflow (Fisher test)
→ fisher_results.tsv
→ format.py
→ input.csv (Chr, Pos)
→ mutation_annotation.py
→ annotated_output.csv
```

---

## Step 1 – Variant Calling with cellSNP-lite

### Example command

```bash
cellsnp-lite \
    -s sensitive_merged.sorted.bam \
    -b sensitive_barcodes.txt \
    -O vcf_sensitive \
    -R reference_snps.vcf \
    --cellTAG CB \
    --UMItag UB \
    --gzip \
    --genotype \
    --minMAF 0.1 \
    --minCOUNT 20 \
    -p 8
```

### Parameter description

- `-s` : Sorted BAM file  
- `-b` : Cell barcode list  
- `-O` : Output directory  
- `-R` : Reference SNP VCF  
- `--cellTAG CB` : Cell barcode tag in BAM  
- `--UMItag UB` : UMI tag in BAM  
- `--genotype` : Perform per-cell genotyping  
- `--minMAF 0.1` : Minor allele frequency threshold  
- `--minCOUNT 20` : Minimum total read depth  
- `-p 8` : Number of threads  

### Output used in downstream analysis

```
cellSNP.cells.vcf.gz
```

---

## Step 2 – Convert VCF to AF Table

Run:

```bash
python trans.py cellSNP.cells.vcf.gz cell_AF_table.tsv
```

### Output format

| ID | Cell | s_reads | t_reads | AF |
|----|------|---------|---------|----|

Where:

- `ID = chrom_pos_ref_alt`  
- `s_reads = ALT reads`  
- `t_reads = total reads`  
- `AF = s_reads / t_reads`  

---

## Step 3 – Fisher Enrichment Analysis in R

Follow the provided R workflow script.

The R analysis performs:

- Merge mutation tables from multiple samples  
- Standardize cell naming by adding sample prefix  
- Assign group labels (resistant or sensitive)  
- Binarize mutation per cell using thresholds:  
  - AF ≥ 0.1  
  - total reads ≥ 10  
  - ALT reads ≥ 3  
- Perform Fisher’s exact test per mutation ID  

### Output files

- `fisher_results.tsv`  
- `heatmap.pdf`  
- `heatmap.png`  

---

## Step 4 – Variant Annotation

### Extract genomic coordinates

```bash
python format.py --input fisher_results.tsv --output input.csv
```

### Annotate variants

```bash
python mutation_annotation.py input.csv annotated_output.csv
```

This step retrieves:

- rsID  
- Gene name  
- Functional consequence  
- Clinical significance  
- ClinVar annotations  

---

# PART II – HLA Typing and Neoantigen Prediction

## Step 1 – HLA Typing using OptiType

Follow the OptiType instructions:

https://github.com/FRED-2/OptiType/issues/141

### Example (official CLI format)

```bash
python /path/to/OptiTypePipeline.py \
  -i sample_fished_1.fastq sample_fished_2.fastq \
  (--rna | --dna) \
  --outdir /path/to/out_dir/
```

Notes:

- `-i` accepts one FASTQ (single-end) or two FASTQs (paired-end).
- Use `--rna` for RNA-seq reads and `--dna` for DNA-seq reads.
- The main output is typically a TSV file containing predicted class I alleles (e.g., `HLA-A*24:02`, `HLA-B*07:02`, `HLA-C*07:02`), which will be used in downstream neoantigen prediction.

---

## Step 2 – Neoantigen Prediction using DIPAN

After mutation identification and HLA typing, follow the DIPAN workflow to infer candidate neoantigens.

DIPAN integrates:

- IPA-derived or mutation-derived candidate events
- Sample-specific HLA alleles inferred by OptiType

### DIPAN output format

A typical DIPAN result table contains columns such as:

- `SYMBOL`
- `Terminal_exon`
- `IPAtype`
- `IPUI`
- `HLA`
- `Peptide`
- `%Rank`

Example:

```text
SYMBOL   Terminal_exon               IPAtype     IPUI   HLA           Peptide     %Rank
PAQR3    chr4:78923401-78923856      Composite   0.528  HLA-A*24:02   RYFPGRYLF   0.001
OXCT1    chr5:41849951-41850029      Composite   0.075  HLA-B*07:02   KPREVRNTL   0.001
NCKAP1   chr2:182980872-182981243    Composite   0.084  HLA-C*07:02   YYFPFVPSF   0.002
```

Where:

- `SYMBOL`: gene symbol
- `Terminal_exon`: genomic coordinates of the terminal exon / IPA region
- `IPAtype`: IPA category (e.g., Composite)
- `IPUI`: IPA usage index (higher indicates stronger usage)
- `HLA`: predicted presenting HLA allele
- `Peptide`: predicted neoantigen peptide sequence
- `%Rank`: binding rank score (lower indicates stronger predicted binding)


# PART III – HLA–Peptide Structural Modeling

## Step 1 – AlphaFold2 Structure Prediction

Use:

```
Neoantigen_visualization/alphafold2.py
```

Modify line 35:

```python
query_sequence = "HLA_SEQUENCE:PEPTIDE_SEQUENCE"
```

### Important notes

- Separate HLA and peptide sequences using ":"  
- Ensure correct spelling: `peptide`  

### Example

```
MAVMAPRTLVLLLSGALALTQTWA:LLFGYPVYV
```

---

## Step 2 – Generate PDB Structure

Running the script produces:

```
predicted_complex.pdb
```

---

## Step 3 – Visualization with PyMOL

Use the provided `pymol_command`.

### Example

```python
load predicted_complex.pdb
color cyan, chain A
color red, chain B
show cartoon
```

This visualizes:

- HLA structure  
- Peptide binding conformation  
- Neoantigen presentation  

---
---

## Step 4 – Variant Structural Modeling (Optional)

In addition to neoantigen–HLA complexes, individual protein variants can also be structurally evaluated.

This step allows visualization of mutation-induced structural changes at the protein level.

### 1. Prepare Mutant Protein Sequence

Generate the mutated protein sequence based on:

- Variant annotation results (from `annotated_output.csv`)
- Amino acid change (e.g., A284T)

Create two sequences:

- Wild-type protein sequence
- Mutant protein sequence

---

### 2. Structure Prediction

Use AlphaFold2 (or ColabFold) to predict structures for:

- Wild-type protein
- Mutant protein

Example:

```bash
colabfold_batch wt_sequence.fasta output_wt/
colabfold_batch mutant_sequence.fasta output_mut/
```

This generates predicted structure files:

```
output_wt/model_1.pdb
output_mut/model_1.pdb
```

---

### 3. Structural Comparison in PyMOL

Load both structures:

```python
load wt_model_1.pdb, WT
load mut_model_1.pdb, MUT
align MUT, WT
```

Highlight mutation site:

```python
select mutation_site, resi 284
show sticks, mutation_site
color red, mutation_site
```

Optional:

```python
color cyan, WT
color orange, MUT
```

This allows:

- Visualization of local conformational changes  
- Inspection of side-chain orientation differences  
- Evaluation of structural plausibility  

---

## Complete Integrated Workflow

```
Single-cell mutation analysis
→ HLA typing (OptiType)
→ DIPAN neoantigen prediction
→ AlphaFold2 structural modeling (HLA–peptide)
→ PDB structure generation
→ PyMOL visualization
→ Variant structural modeling (optional)
```

---

## Summary

This repository integrates:

- Single-cell mutation enrichment analysis  
- HLA genotype inference  
- Neoantigen prediction  
- HLA–peptide structural modeling  
- Variant structural comparison  

It enables both statistical mutation analysis and structural validation of neoantigen presentation and mutation effects.
