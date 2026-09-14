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
  <img src="docs/figures/resist_overview.png"
       alt="RESIST workflow: data collection and curation; characterization, transcriptional, regulatory, and immunogenomic analyses; and an interactive database with downloadable results."
       width="1100"/>
</p>

**Figure 1. RESIST workflow.** Curated single-cell and spatial datasets support
four complementary analysis modules, connecting cellular characterization with
transcriptional, regulatory, and immunogenomic analyses. The web interface
supports interactive browsing, visualization, and access to data and results.
The figure summarizes the resource; individual analyses require different inputs
and are selected according to data availability.

## Analysis modules

| Module | Biological focus | Analyses and products |
|-----|---|---|
| **[A · Characterization](docs/module-a.md)** | Cellular composition and organization | Existing UMAP embeddings, cluster composition by response, CellChat communication patterns, and spatial distributions |
| **[B · Transcriptional](docs/module-b.md)** | Expression programs associated with resistance | Differentially expressed genes, GO/KEGG/Hallmark enrichment, intratumor heterogeneity, EMT scores, and LINCS signature connectivity |
| **[C · Regulatory](docs/module-c.md)** | Candidate regulatory associations | Variant-analysis components and enrichment of RNA-binding-protein and miRNA target sets |
| **[D · Immunogenomic](docs/module-d.md)** | RNA-processing and immune-presentation features | Alternative polyadenylation, HLA typing, and candidate peptide–MHC structural analysis |

The module guides describe each script's inputs, dependencies, and scope.
Availability of individual components, including the TF motif panel shown in
Figure 1, is explained under [Implementation scope](#implementation-scope).

## Quick start

Run the analyses on an HPC system with suitable storage and memory. Begin from
the repository root. If using a downloaded package ZIP, extract it on the HPC,
enter its directory, and skip the two Git checkout commands below:

```bash
git clone https://github.com/QSong-github/RESIST.git RESIST
cd RESIST
conda env create -f config/setup/environment.yml
conda activate resist
Rscript config/setup/install_r_packages.R
```

Follow [tutorial Steps 2–3](TUTORIAL.md#step-2--download-the-example-directly-onto-the-hpc)
to download **`GSE104987_seurat_afterAnno.RDS`** directly onto the HPC, unpack the
bundled references, and inspect the input. The example profile writes to
`data/results/example/`. Its A–C walkthrough does not require the large LINCS
database or raw sequencing files.

**Each command below runs the analyses described in its section.** A, B, and C
use the shared example; C requires the differential-expression table produced by
B. D requires a separately prepared APA cohort. Additional analyses within each
module are described in the linked module guides.

### A. Characterization

```bash
bash run_A.sh --config config/example.yaml
```

**Produces one four-panel overview figure** containing:

- Three views of the existing UMAP, colored by cell cluster, response condition,
  and annotated cell type.
- A grouped bar chart comparing cluster composition between sensitive and
  resistant cells.

The figure is saved as **PDF and PNG** in `data/results/example/A/`.
It uses the input object's existing coordinates and annotations. Communication
and spatial analyses require additional inputs; see the [Module A guide](docs/module-a.md).

### B. Transcriptional analysis

```bash
bash run_B.sh --config config/example.yaml
```

**Produces differential-expression results and a heterogeneity comparison:**

- **Full DEG table (CSV):** tested genes and their expression differences between
  resistant and sensitive cells, within each eligible cell type.
- **Significant DEG table (CSV):** genes passing the analysis's
  significance and fold-change criteria.
- **Tumor-cell volcano plot (PDF and PNG):** expression effect sizes and adjusted
  p-values.
- **Intratumor-heterogeneity box plot (PDF and PNG):** comparison of PCA-based
  heterogeneity scores between response conditions.

Files are saved in `data/results/example/B/`. Positive log fold changes indicate
higher expression in resistant cells. Keep the full DEG table for Module C.
Pathway, EMT, and drug-signature analyses are described in the
[Module B guide](docs/module-b.md).

### C. Regulatory analysis

**Run B first**, using the same example configuration and results directory.

```bash
bash run_C.sh --config config/example.yaml
```

**Produces RNA-binding-protein (RBP) target-enrichment results:**

- **RBP-enrichment table (CSV):** overlap of up- and down-regulated tumor genes
  with the bundled RBP target sets, including enrichment statistics.
- **RBP circle plots (PDF and PNG):** summaries of eligible enriched target sets
  for each expression direction.

Files are saved in `data/results/example/C/`. This command reads the full DEG
table generated by B. Circle plots are produced only when significance and
count thresholds are met, so an enrichment table can exist without a plot.
The miRNA and variant branches have separate requirements in the
[Module C guide](docs/module-c.md).

### D. Immunogenomic analysis

Prepare a separate APA cohort following
[tutorial Step 7](TUTORIAL.md#step-7--module-d-with-a-separately-prepared-apa-cohort).
It requires scUTRquant TXS outputs, filtered 10x matrices, an annotation Seurat
object, and the matching GTF; the shared GSE104987 example does not supply these
inputs. Then run:

```bash
bash run_D.sh --config config/config.yaml --apa-config config/apa_config.yaml
```

**Produces alternative-polyadenylation (APA) summaries and comparisons:**

- **Gene-by-cell relative-expression matrices (RDS)** and **per-cell mean
  relative-expression summaries (CSV)**.
- **Barcode-to-annotation mappings and mapping diagnostics (CSV)** for checking
  how quantified cells match the annotation object.
- **Differential APA results (CSV)** for eligible comparisons.
- **ECDF and violin plots (PDF and PNG)** for eligible response groups and cell
  types.

Files are saved in `data/results/D/APA/<cohort_id>/`. Review mapping diagnostics
before interpreting the comparisons. HLA typing and peptide–MHC structural
analysis use separate commands and inputs in the [Module D guide](docs/module-d.md).

### Check inputs and select additional analyses

Use `--check` to inspect prerequisites without executing an analysis, and
`--list` to inspect the available analyses and expected products:

```bash
bash run_A.sh --config config/example.yaml --check
bash run_B.sh --list
```

The [tutorial](TUTORIAL.md) provides output checkpoints and extension commands.
The module guides explain individual analysis selections with `--steps`; the
[output catalog](docs/outputs.md) lists their filenames and eligibility conditions.

## Repository layout

```text
RESIST/
├── README.md                  resource overview and quick start
├── TUTORIAL.md                step-by-step example-data walkthrough
├── run_A.sh … run_D.sh        module entry points
├── config/                   paths, sample sheets, dependencies, SLURM template
├── scripts/                  A–D analyses, shared R functions, utilities
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

Each module runner records selected steps, configuration and script hashes,
package versions, step logs, and the files written during the run under
`<results>/logs/`. Failures stop execution, and steps that produce no new output
are identified explicitly. Use a separate results path when comparing analyses
with different inputs or settings. The [HPC guide](docs/hpc.md) covers submission
and environment recording.

## Citation and acknowledgements

Please cite the RESIST resource, the original datasets, and the analysis tools
used in your work. Resource citation metadata is provided in
[CITATION.cff](CITATION.cff). We acknowledge the investigators who generated the
underlying datasets and the developers of the community software used by RESIST.
