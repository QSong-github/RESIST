# Package update notes

The README quick start now explains what each A–D command produces, with file
formats, result directories, and input dependencies. Script identifiers remain
in the technical guides and command-line step catalog.

## Documentation changes

- Replaced numbered-script introductions with descriptions of UMAP/composition,
  differential expression/heterogeneity, RBP enrichment, and APA outputs.
- Removed script IDs throughout the public README, including its optional-step
  example and implementation-scope text. Missing components remain disclosed by
  their function, with technical details in the module guides.
- Clarified that the quick start covers selected analyses within each module.
  Additional module capabilities have separate inputs and commands.
- Stated that C uses the full DEG table produced by B and that D needs an
  independently prepared APA cohort.
- Kept plot eligibility and barcode-mapping checks visible beside the outputs.
- Made the tutorial's introductory route use analysis names; detailed script
  references remain available for traceability.
- Removed release-number prose from active guides. The package identity is
  recorded in `VERSION`, `CITATION.cff`, and `config/steps.json`, which agree.

## Package contents and compatibility

The complete package retains four top-level folders: `config`, `scripts`, `data`,
and `docs`. It includes the supplied workflow figure, example-data tutorial,
all analysis scripts, configuration templates, shared code, and the compressed
reference bundle. No analysis scripts, default selections, output filenames,
statistical definitions, or reference contents changed in this update.

Software, reference-database, model, and sequencing-chemistry identifiers retain
their correct values. Historical source filenames remain in provenance records
and source-origin comments. Earlier documentation remains in the clearly labeled
[source-history archive](source-history.tar.gz).

The [migration map](migration.tsv) accounts for all 188 original supplied files,
including 119 unchanged reference files inside the compressed archive. Current
target hashes and the root `SHA256SUMS.txt` are regenerated for the package.

## Validation boundary

Checks cover syntax, launcher interfaces, release metadata, README descriptions,
local links, source traceability, and packaging integrity. Scientific workflows,
environment installation, scheduler jobs, and large external downloads remain
deferred to the HPC. See the [validation record](validation.md) and
[HPC acceptance guidance](hpc.md#first-run-acceptance).
