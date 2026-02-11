# RESIST

RESIST provides an integrated analysis framework for investigating multi-layer regulatory mechanisms underlying cancer drug resistance from single-cell sequencing data.

This repository contains analysis modules for:

- Alternative Polyadenylation (APA) analysis  
- Variant detection and mutation profiling  
- Neoantigen prediction  

The framework operates on processed 10x single-cell data (e.g., BAM files, barcode lists, count matrices) and builds upon established community tools.

---

## Implemented Analysis Modules

### 1. Alternative Polyadenylation (APA) Analysis

APA analysis is performed using:

- **scUTRquant**  
  https://github.com/mortazavilab/scUTRquant  

scUTRquant is used to quantify proximal and distal poly(A) site usage from single-cell RNA-seq data. Downstream analyses in RESIST compute relative expression (RE), differential APA, and cell type–specific APA shifts.

---

### 2. Variant Detection

Single-cell SNP detection is performed using:

- **cellSNP-lite**  
  https://github.com/single-cell-genetics/cellsnp-lite  

cellSNP-lite is used to generate allele depth (AD) and total depth (DP) matrices from aligned BAM files. These matrices are subsequently used for mutation profiling and downstream integration analyses.

---

### 3. HLA Typing

HLA typing is performed using:

- **OptiType**  
  https://github.com/FRED-2/OptiType  

OptiType is used to infer HLA genotypes from RNA-seq data for downstream neoantigen prediction.

---

### 4. Neoantigen Prediction

Neoantigen inference is performed using:

- **DIPAN**  
  (Insert GitHub link if public)

DIPAN integrates mutation information and HLA typing results to predict candidate neoantigens for resistant and sensitive groups.

---

## Required Inputs

Before running the variant or neoantigen modules, users should prepare:

- Sorted and indexed BAM files  
- Cell barcode lists  
- Reference VCF file (for variant calling)  
- HLA typing results (for neoantigen prediction)

Alignment, BAM preprocessing, and dataset-specific filtering depend on sequencing platform and computing environment, and are therefore not included in this repository.

---

## Output

The framework generates:

- AD/DP variant matrices  
- Differential APA results  
- Candidate neoantigen lists  
- Integrated multi-layer regulatory summaries  

---

## Notes

RESIST provides analysis modules and reference commands for key steps.  
Users are responsible for adapting upstream preprocessing (e.g., alignment, BAM merging) to their specific datasets and computational environments.
