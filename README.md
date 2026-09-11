# RESIST · v3_02

**Single-cell and spatial omics analyses of cancer drug resistance**

RESIST organizes analyses into four biological modules: cellular characterization
(A), transcriptional programs (B), regulatory associations (C), and immunogenomic
features (D). The resource is available at [resist.website](https://resist.website/).

This release simplifies the supplied v2-2 package from **12 to 4 top-level folders**,
adds a worked tutorial, and makes each command's inputs and outputs explicit.
It is prepared for **subsequent HPC testing**; the full workflows have not been
executed or validated as part of this release.

**Start with [TUTORIAL.md](TUTORIAL.md)** for the real GSE104987 example, including
download, input inspection, execution order, expected filenames, and interpretation.
Download data and expand reference archives on the HPC, where they will be used.

## Quick start — A, B, C, D

Run these commands from the extracted `RESIST_v3_02` directory **on the HPC**.
Complete [tutorial Steps 1–3](TUTORIAL.md#step-1--prepare-the-hpc-environment)
first: activate the environment, download the example, and unpack the small
bundled reference archive. The selected example profile writes results to
`data/results/example/`.

### A. Characterization

```bash
bash run_A.sh --config config/example.yaml
```

Runs **A1**. Produces a four-panel UMAP/cluster-composition figure in **PDF and PNG**,
using the existing embeddings and annotations.

### B. Transcriptional analysis

```bash
bash run_B.sh --config config/example.yaml
```

Runs **B1 and B4**. Produces **full and significant DEG CSV tables**, a **tumor-cell
volcano plot**, and an **intratumor-heterogeneity box plot**; figures are PDF/PNG.
The contrast is **resistant minus sensitive**. Run this before C.

### C. Regulatory analysis

```bash
bash run_C.sh --config config/example.yaml
```

Runs **C7 and C8** using B1's full DEG table and the bundled RBP reference.
Produces an **RBP-enrichment CSV** and, when the inherited significance/count
criteria are met, **RBP circle plots** in PDF/PNG. No significant circle plot is
a possible outcome; inspect the table and log.

### D. Immunogenomic analysis — separate inputs required

The GSE104987 Seurat example does **not** contain the scUTRquant, 10x matrix, and
annotation/GTF inputs needed for APA. Prepare those inputs using
[tutorial Step 7](TUTORIAL.md#step-7--module-d-with-a-separately-prepared-apa-cohort),
then run:

```bash
bash run_D.sh --config config/config.yaml --apa-config config/apa_config.yaml
```

Runs **D3**. Produces **gene-by-cell relative-expression RDS matrices**, **per-cell
summary CSVs**, **barcode-mapping diagnostics**, **differential APA tables**, and
**ECDF/violin figures**, subject to eligible data. D1/D2 preparation, D4 HLA typing,
and D5/D6 structural work have their own inputs and commands in the module guide.

These four entry points are the only quick-start execution route. Their defaults
are deliberately small and explicit; they do not imply that every analysis in a
module can be run from the same example file.

## Additional analyses

Use `--list` to see available steps and outputs without loading data. Use `--check`
on the HPC to inspect prerequisites without executing analyses. Use `--steps` to
select additional work explicitly, for example:

```bash
bash run_B.sh --list
bash run_B.sh --config config/example.yaml --steps B5 --check
bash run_B.sh --config config/example.yaml --steps B5
```

| Module | Additional analyses | Guide |
|---|---|---|
| A | CellChat communication; study-specific spatial plots | [Module A](docs/module-a.md) |
| B | GO/KEGG/Hallmark enrichment; EMT; LINCS connectivity | [Module B](docs/module-b.md) |
| C | miRNA enrichment; variant calling/annotation; missing C2 and TF implementation documented | [Module C](docs/module-c.md) |
| D | APA preparation; HLA typing; candidate peptide–MHC modeling | [Module D](docs/module-d.md) |

## Four folders, one place to start

```text
RESIST_v3_02/
├── README.md                 overview and A–D quick start
├── TUTORIAL.md               detailed real-data walkthrough
├── run_A.sh … run_D.sh       consistent module entry points
├── config/                  paths, sample sheets, dependencies, SLURM template
├── scripts/                 A/, B/, C/, D/, shared code, and utilities
├── data/                    compressed references and small example templates
└── docs/                    module guides, output catalog, HPC and release notes
```

`data/example/`, `data/input/`, `data/reference/`, and `data/results/` are created
on the HPC as needed. Large files can live outside this package: set absolute paths
in `config/config.yaml`. Relative data paths resolve against the package root.
The original analysis filenames and A1–D6 numbering are retained for traceability.

## Results and reproducibility

[The output catalog](docs/outputs.md) lists products and their prerequisites.
Each runner records the selected steps, configuration copy and hash, R/package
versions, step logs, and new/updated output filenames and SHA-256 hashes under
`<results>/logs/`. A failed step stops the runner. A step that writes no output
is reported explicitly.

Read [HPC instructions](docs/hpc.md), [reference-data instructions](docs/references.md),
and [scientific scope and limitations](docs/methods-and-limitations.md) before
extending the example to new cohorts. The inherited methods remain exploratory
associations; this reorganization does not establish causal resistance mechanisms
or validate candidate drugs, variants, or neoantigens.

## Citation and release record

Cite the RESIST resource, the original datasets, and the tools used in your run.
[CITATION.cff](CITATION.cff) preserves the supplied resource metadata; publication
and author details should be completed by the project authors before publication.
See [release notes](docs/release-notes.md), the [file migration map](docs/migration.tsv),
and the [validation record](docs/validation.md) for what changed and what remains
to test on the HPC.
