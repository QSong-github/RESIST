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

### Example

```bash
OptiTypePipeline.py \
   -i tumor_R1.fastq tumor_R2.fastq \
   -r \
   -o optitype_output \
   --dna
```

### Output

```
optitype_result.tsv
```

Extract predicted HLA alleles (e.g., HLA-A*02:01).

---

## Step 2 – Neoantigen Prediction using DIPAN

After mutation identification and HLA typing:

- Generate mutant peptide sequences  
- Predict HLA binding affinity  
- Identify candidate neoantigens  

### Output

```
neoantigen_candidates.tsv
```

---

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

# Complete Integrated Workflow

```
Single-cell mutation analysis
→ HLA typing (OptiType)
→ DIPAN neoantigen prediction
→ AlphaFold2 structural modeling
→ PDB structure generation
→ PyMOL visualization
```

---

# Summary

This repository integrates:

- Single-cell mutation enrichment analysis  
- HLA genotype inference  
- Neoantigen prediction  
- Structural modeling of HLA–peptide complexes  

It enables both statistical mutation analysis and structural validation of neoantigen presentation.
