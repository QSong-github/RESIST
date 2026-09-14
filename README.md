<h1 align="center">RESIST</h1>

<p align="center">
  <b>A single-cell and spatial omics resource for studying regulatory mechanisms<br/>
  underlying cancer drug resistance</b>
</p>

<p align="center">
  <a href="https://resist.website/">Explore the database</a> ·
  <a href="TUTORIAL.md">Example-data tutorial</a> ·
  <a href="docs/outputs.md">Output catalog</a> ·
  <a href="#external-tools">External tools</a> ·
  <a href="docs/hpc.md">HPC guide</a>
</p>

## Overview

Cancer drug resistance involves changes in cellular composition, transcriptional
programs, regulatory interactions, and immunogenomic features. RESIST brings
together single-cell RNA sequencing and spatial transcriptomics to investigate
these complementary aspects of treatment response.

The resource comprises **81 scRNA-seq datasets and 7 spatial transcriptomics
datasets**, covering **522 single-cell samples, 96 spatial samples, 13 cancer
types, and 59 drug regimens**. This repository contains analysis code organized
into four modules, with configuration files, reference-data instructions, and a
worked example. Processed data and results can be explored through the
[RESIST database](https://resist.website/).

## Workflow

<p align="center">
  <a href="docs/figures/Fig1.pdf"><img src="docs/figures/resist_overview.png"
       alt="RESIST workflow: data collection and curation; characterization, transcriptional, regulatory, and immunogenomic analyses; and an interactive database with downloadable results."
       width="1100"/></a>
</p>

**Figure 1. RESIST workflow.** Curated single-cell and spatial datasets support
four complementary analysis modules, connecting cellular characterization with
transcriptional, regulatory, and immunogenomic analyses. The web interface
supports interactive browsing, visualization, and access to data and results.
The figure summarizes the resource; individual analyses require different inputs
and are selected according to data availability.
[Download Figure 1 (PDF)](docs/figures/Fig1.pdf).

## Analysis modules

| Module | Biological focus | Key analyses |
|---|---|---|
| **[A.&nbsp;Characterization](docs/module-a.md)** | Cellular organization | UMAP; cluster composition;<br>CellChat; spatial distribution |
| **[B.&nbsp;Transcriptional](docs/module-b.md)** | Resistance-associated expression | DEGs; pathway enrichment;<br>ITH; EMT; drug connectivity |
| **[C.&nbsp;Regulatory](docs/module-c.md)** | Regulatory associations | Variant-analysis components;<br>RBP and miRNA enrichment |
| **[D.&nbsp;Immunogenomic](docs/module-d.md)** | RNA processing and immunity | APA; HLA typing;<br>peptide–MHC structure |

See the module guides for inputs and methods, and
[Implementation scope](#implementation-scope) for component availability.

## Quick start

On the HPC, [install dependencies and prepare the example](TUTORIAL.md), then run
these commands from the repository root. **Run B before C.** D requires
[a separate APA cohort](TUTORIAL.md#step-7--module-d-with-a-separately-prepared-apa-cohort).

### A. Characterization

```bash
bash scripts/run_A.sh --config config/example.yaml
```

**Output:** UMAP and cluster-composition figure (PDF/PNG) in `data/results/example/A/`.

### B. Transcriptional analysis

```bash
bash scripts/run_B.sh --config config/example.yaml
```

**Output:** Full/significant DEG tables (CSV), tumor-cell volcano and
intratumor-heterogeneity plots (PDF/PNG) in `data/results/example/B/`.

### C. Regulatory analysis

```bash
bash scripts/run_C.sh --config config/example.yaml
```

**Output:** RBP-enrichment table (CSV) and eligible circle plots (PDF/PNG) in
`data/results/example/C/`.

### D. Immunogenomic analysis

```bash
bash scripts/run_D.sh --config config/config.yaml --apa-config config/apa_config.yaml
```

**Output:** APA matrices (RDS), summary/mapping/differential tables (CSV), and
eligible ECDF/violin plots (PDF/PNG) in `data/results/D/APA/<cohort_id>/`.

Detailed steps and output definitions: [tutorial](TUTORIAL.md) ·
[output catalog](docs/outputs.md).

## Repository layout

```text
RESIST/
├── README.md                 resource overview and quick start
├── TUTORIAL.md               step-by-step example-data walkthrough
├── config/                   paths, sample sheets, dependencies, SLURM template
├── scripts/                  run_A.sh … run_D.sh, A–D analyses, shared helpers
├── data/                     compressed references and small input templates
└── docs/                     workflow figure, module guides, and technical notes
```

Input, expanded-reference, and result directories are created on the HPC as
needed. See the [repository guide](docs/repository-guide.md) for the script
organization and configuration conventions.

## Data and configuration

`config/config.yaml` defines the main data, reference, and result locations:

```yaml
paths:
  data: data/input
  spatial_data: data/spatial
  ref_data: data/reference
  results: data/results
```

Relative paths resolve against the repository root; absolute paths allow inputs
and references to remain on shared HPC storage. Use `config/example.yaml` for the
worked example and the separate APA configuration/sample sheet for Module D.

The [input-data guide](docs/input-data.md) describes the Seurat metadata and
branch-specific inputs. The [reference-data guide](docs/references.md) distinguishes
the compressed bundle from external resources such as LINCS and miRDB. Download
large references only for analyses that need them.

## External tools

RESIST uses the following key analysis tools. A–D refer to the modules above;
some tools support optional branches beyond the quick start. **Upstream** and
**separate** workflows run outside the module launchers.

| Tool | Used by | Role in the workflow | Official source |
|---|---|---|---|
| **Seurat** | Shared R environment; A, B, D analyses | Single-cell and spatial objects, existing embeddings, differential expression, and annotation handling | [Documentation](https://satijalab.org/seurat/) |
| **scType** | Upstream cell-type annotation | Assign cell types before RESIST analysis; the launchers read saved annotations | [GitHub](https://github.com/IanevskiAleksandr/sc-type) |
| **CellChat** | A · cell–cell communication | Infer and compare ligand–receptor communication networks between response conditions | [GitHub](https://github.com/jinworks/CellChat) |
| **clusterProfiler** | B · pathway enrichment | GO, KEGG, and gene-set enrichment of differentially expressed genes | [Bioconductor](https://bioconductor.org/packages/clusterProfiler/) |
| **msigdbr** | B · pathway enrichment | Retrieve MSigDB Hallmark gene sets for the enrichment analysis | [Documentation](https://igordot.github.io/msigdbr/) |
| **GSVA** | B · EMT scoring | Calculate EMT gene-set scores from expression data | [Bioconductor](https://bioconductor.org/packages/GSVA/) |
| **cmapR** | B · drug-signature analysis | Handle the prepared LINCS annotated-matrix object; RESIST calculates connectivity scores | [GitHub](https://github.com/cmap/cmapR) |
| **biomaRt** | C · miRNA enrichment | Map gene identifiers used to connect DEG lists with miRNA target predictions | [Bioconductor](https://bioconductor.org/packages/biomaRt/) |
| **cellSNP-lite** | C · variant analysis | Count reference and alternative alleles at selected SNP loci from BAM files and cell barcodes | [GitHub](https://github.com/single-cell-genetics/cellsnp-lite) |
| **scUTRquant** | Upstream quantification for D | Quantify single-cell 3′ UTR isoforms; RESIST consumes the prepared TXS outputs | [GitHub](https://github.com/Mayrlab/scUTRquant) |
| **OptiType** | D · separate HLA typing | Predict HLA-A, HLA-B, and HLA-C alleles; the wrapper expects an `OptiTypePipeline.py` installation | [GitHub](https://github.com/FRED-2/OptiType) |
| **DIPAN** | D · separate candidate generation | Identify intronic-polyadenylation-derived neoantigen candidates from RNA-seq; no automated RESIST wrapper is included | [GitHub](https://github.com/YY-TMU/DIPAN) |
| **ColabFold / AlphaFold2** | D · separate structural prediction | Model candidate peptide–MHC complexes in a dedicated Colab or GPU environment | [ColabFold](https://github.com/sokrypton/ColabFold), [AlphaFold2](https://github.com/google-deepmind/alphafold) |
| **PyMOL** | D · separate structural visualization | Inspect predicted complexes and render peptide–MHC interaction figures | [PyMOL](https://pymol.org/) |

For new Seurat inputs, store reviewed scType labels in `cell_type` or `celltype`,
following the [input-data guide](docs/input-data.md). The shared example is
already annotated, so its quick start does not require running scType again.

Full dependencies, including supporting libraries, are listed in the
[environment specification](config/setup/environment.yml),
[R package installer](config/setup/install_r_packages.R), and
[Python requirements](config/setup/requirements.txt). Follow the module guides
for installation and the [reference-data guide](docs/references.md) for HPC data
preparation. Cite the tools used in your analysis using their official guidance.

## Implementation scope

The example begins with an annotated, processed Seurat object. Module A displays its
existing annotations and embeddings. Communication, spatial, variant, and
immunogenomic analyses require additional data and external tools.

The **VCF-to-allele-frequency conversion** and the **TF motif-enrichment
implementation** are not included in the supplied code. Other branches retain
method- and cohort-specific assumptions documented in the module guides. Full
workflow execution and numerical validation remain to be completed on the HPC;
the [validation record](docs/validation.md) states the checks performed so far.

Use the [methods and limitations](docs/methods-and-limitations.md) when interpreting
results or adapting the pipeline to new studies. Cell-level contrasts and
enrichment results support exploratory associations; candidate mechanisms,
compounds, and neoantigens require appropriate experimental validation.

## Reproducibility

The R installer uses a fixed CellChat source revision and records resolved package
versions and software checks. Each module runner records selected steps,
configuration and script hashes, installed package versions, step logs, and
the files written during the run under
`<results>/logs/`. Failures stop execution, and steps that produce no new output
are identified explicitly. Use a separate results path when comparing analyses
with different inputs or settings. The [HPC guide](docs/hpc.md) covers submission
and environment recording, including [CellChat installation troubleshooting](docs/hpc.md#cellchat-installation).

## Citation and acknowledgements

Please cite the RESIST resource, the original datasets, and the analysis tools
used in your work. Resource citation metadata is provided in
[CITATION.cff](CITATION.cff). We acknowledge the investigators who generated the
underlying datasets and the developers of the community software used by RESIST.
