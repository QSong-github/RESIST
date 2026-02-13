# RESIST – Neoantigen and Mutation Analysis Workflow

This repository contains a complete workflow for:

- Single-cell mutation enrichment analysis  
- HLA typing  
- Neoantigen prediction  
- HLA–peptide structural modeling  
- PyMOL visualization  
## External Dependencies

RESIST builds upon the following established community tools:

### Alternative Polyadenylation

- **scUTRquant**  
  https://github.com/Mayrlab/scUTRquant  

---

### Variant Detection

- **cellSNP-lite**  
  https://github.com/single-cell-genetics/cellsnp-lite  

---

### HLA Typing

- **OptiType**  
  https://github.com/FRED-2/OptiType  

---

### Neoantigen Prediction

- **DIPAN**  
  https://github.com/YY-TMU/DIPAN  

---

### Protein Structure Prediction

- **AlphaFold2**  
  https://github.com/deepmind/alphafold  

- **ColabFold** (lightweight AlphaFold implementation)  
  https://github.com/sokrypton/ColabFold  

---

### Structure Visualization

- **PyMOL**  
  https://pymol.org/

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
# PART IV – Alternative Polyadenylation (APA) Analysis

## Overview

The APA module in RESIST quantifies alternative polyadenylation events from single-cell RNA-seq data and evaluates differential 3′ UTR usage between groups.

APA analysis is performed using:

- **scUTRquant**  
  https://github.com/Mayrlab/scUTRquant  

Downstream processing computes relative expression (RE), differential APA events, and cell type–specific APA shifts.

---

## Step 1 – Identify Library Kit Version

Before running scUTRquant, determine the 10x Genomics library kit version.

Run:

```bash
python find_kit_version.py list.txt
```

### Input

`list.txt`

A text file containing one sample directory per line.

Example:

```
/path/to/sample1
/path/to/sample2
```

### Output

`kit.csv`

This file records inferred kit versions for each sample.

---

## Step 2 – Configure scUTRquant

Fill in the following files:

- `config.yaml`
- `sample_sheet.csv`

Using:

- Kit information from `kit.csv`
- Sample metadata
- Reference transcriptome path
- Output directory

Example configuration files are provided in the repository.

---

## Step 3 – Run scUTRquant

After configuration:

```bash
cd scUTRquant
snakemake --use-conda --configfile examples/config.yaml
```

This step performs:

- Proximal and distal poly(A) site quantification  
- Transcript-level abundance estimation  
- APA isoform assignment  

---

## Step 4 – Locate Output Files

After successful execution, results are generated under the configured output directory (e.g., `data/`).

The primary output file used for downstream analysis is:

```
GSE261898.txs.Rds
```

This file contains transcript-level quantification results.

---

## Step 5 – Downstream APA Processing and Visualization

Use the provided R script:

```
APA_processing_graph.R
```

Run in R:

```r
source("APA_processing_graph.R")
```

This script performs:

- Relative expression (RE) calculation  
- Differential APA analysis  
- Cell type–specific APA comparison  
- APA shift visualization  

---

## APA Output

The APA module generates:

- RE matrices  
- Differential APA statistics  
- APA shift visualizations  
- Processed APA summary tables  

---

## Integrated APA Workflow

```
Processed 10x data
→ find_kit_version.py
→ config.yaml + sample_sheet.csv
→ scUTRquant run
→ GSE261898.txs.Rds
→ APA_processing_graph.R
→ APA results and figures
```


