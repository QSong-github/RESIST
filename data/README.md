# Data for HPC use

This folder ships only the compressed RESIST reference bundle and small templates.
Download the real example directly to the HPC and follow [TUTORIAL.md](../TUTORIAL.md).

- `reference_bundle.tar.gz`: unpack on the HPC with `tar -xzf data/reference_bundle.tar.gz -C data`.
- `templates/`: sequencing metadata, scUTRquant example configuration, and structure query examples; review and replace placeholders before use.
- `example/`, `input/`, `spatial/`, `reference/`, and `results/`: runtime locations created as needed, not populated with large data in the package.

Large references can live outside the package via absolute paths in the config.
See [reference instructions](../docs/references.md).
