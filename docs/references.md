# Reference data and storage

Expand and download reference data **on the HPC**, not on a storage-limited laptop.
The package keeps the original supplied references in one compressed archive;
large external databases are not included or downloaded automatically.

## Bundled reference archive

`data/reference_bundle.tar.gz` is approximately 16.8 MB compressed and 75.3 MB
expanded. It contains all 119 supplied reference-data files. This page replaces
the original reference README. A per-file migration/checksum record is in
[migration.tsv](migration.tsv).

```bash
# Package-local storage, matching config/example.yaml and config/config.yaml:
tar -xzf data/reference_bundle.tar.gz -C data
```

For shared storage:

```bash
mkdir -p /your/HPC/shared/resist
tar -xzf data/reference_bundle.tar.gz -C /your/HPC/shared/resist
```

Then set `paths.ref_data: /your/HPC/shared/resist/reference` in your own config.
The archive expands into a `reference/` directory. Do not unpack it repeatedly
into each analysis cohort or duplicate it across project copies.

| Configuration key | Supplied content | Used by |
|---|---|---|
| `rbp_target_dir` | `All_RBP_TargetGene/`, one table per human RBP | C7 |
| `emt_gmt_human` | Human Hallmark EMT v2024.1 GMT | B5 |
| `biomart_human`, `biomart_mouse` | Cached RDS objects | C9 |
| `gene_map_human`, `gene_map_mouse` | Gene-identifier mapping CSVs | C9 |
| Other supplied references | Mouse EMT GMT, epithelial markers, mouse RBP list, drug metadata | Preserved for provenance; not all are used by active code |

Cached objects are not guaranteed to be compatible with every future package
version. Preserve the checksum and reference provenance used for your HPC run.

## Optional external references — download only for selected analyses

| Analysis | Missing resource | Storage/compute implication |
|---|---|---|
| B2/B3 | GO/KEGG/Hallmark package/database resources | Package versions and network/database access matter |
| B6/B7 | GSE70138 LINCS Level 4 data, then a locally serialized `GSE70138_LINCS_Level4.rds` | Multi-GB source and large in-memory object; high-memory HPC job |
| C9/C10 | `miRDB_v6.0_prediction_result.txt` | Separate prediction table; the source estimates about 57 MB compressed |
| C1 | Reference SNP VCF matching the alignment genome | Supplied by the analysis owner |
| D3 | scUTRquant target GTF and sample-level TXS outputs | Depends on the target and cohort |
| D4/D5 | HLA/tool-specific references or model resources | External tools and dedicated environments |

**miRDB:** obtain the specified v6.0 prediction table from
[the miRDB download page](https://mirdb.org/download.html), decompress it on the
HPC, and place it under `paths.ref_data` using the exact filename in the config.
Confirm file format and version before running C9.

**LINCS:** consult [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)
and obtain the Level 4 file corresponding to the source pipeline's
`n78980x22268_2015-06-30` signature matrix. The source instructions inconsistently
pairs a `.gct` filename with a `parse_gctx()` example. Confirm whether the downloaded
file is text GCT or HDF5 GCTX and use the parser appropriate to that format and your
installed cmapR version. Do not merely change its extension. Inspect the matrix
orientation and row/column annotation required by B6, serialize the parsed object
as `GSE70138_LINCS_Level4.rds`, and record source/checksum/parser versions. This
release does not claim that conversion has been validated.

The large LINCS download and conversion are intentionally outside the default
example route. The absence of this file should affect B6, not prevent A1/B1/B4/C7.
