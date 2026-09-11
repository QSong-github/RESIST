# HPC preparation and submission

Transfer this compact release, install dependencies in an HPC environment, and
place large data/reference files in project or shared storage. Nothing in the
quick start requires storing those large files on a laptop.

## Prepare once

1. Extract the release in persistent project space.
2. Create and activate the R/Python environment using tutorial Step 1.
3. Download the example and unpack references on the HPC.
4. Use `config/example.yaml` for the worked example; copy `config/config.yaml`
   for another cohort and edit only that cohort's paths/settings.
5. Run `--check` inside the environment and on an appropriate compute node.

The supplied software specifications are not a tested lockfile. Capture the
resolved installation for the run:

```bash
conda env export --no-builds > /your/HPC/project/resist-environment-resolved.yml
Rscript -e 'writeLines(capture.output(sessionInfo()), "/your/HPC/project/resist-sessionInfo.txt")'
```

## Submit one module at a time

From the package root, with the environment activated:

```bash
sbatch --export=ALL,RESIST_HOME="$PWD" config/slurm/submit_module.sbatch \
  A --config config/example.yaml
sbatch --export=ALL,RESIST_HOME="$PWD" config/slurm/submit_module.sbatch \
  B --config config/example.yaml
```

Submit C **after B succeeds**. To enforce that relationship in SLURM:

```bash
job_b=$(sbatch --parsable --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch B --config config/example.yaml)
sbatch --dependency="afterok:${job_b%%;*}" --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch C --config config/example.yaml
```

The example submission template requests one node, four CPUs, 16 GB RAM, and four
hours as initial placeholders, not measured requirements. Edit account/partition
settings for your cluster. It inherits the activated environment and does not
hard-code an institution-specific R module. Check that both `Rscript` and `python3`
are available inside the job.

B4 memory grows quadratically with tumor-cell count. B6 needs a much larger
allocation; the v2-2 comments describe tens of GB and previously used a 90 GB
submission request. An example override is:

```bash
sbatch --mem=90G --time=12:00:00 --export=ALL,RESIST_HOME="$PWD" \
  config/slurm/submit_module.sbatch B --config config/example.yaml --steps B6,B7
```

Measure actual peak memory and elapsed time with the scheduler and adjust the
request. The `--cpus-per-task` allocation does not itself make every R step
parallel. Submit D3 only once its independent cohort inputs are prepared; GPU
structural prediction should use the institution's separate GPU workflow.

## Paths, logs, and repeated runs

Every root runner can be called from another working directory by absolute path.
Use an absolute `--config` path in that case. The wrapper passes the package root
to the scripts; configured relative data paths still resolve against that root.
Direct R entry points also retain their root-discovery bootstrap, but direct
invocation bypasses the launcher's prerequisite checks and provenance logs.

Each runner logs to `<results>/logs/` and writes analysis products to
`<results>/A`, `/B`, `/C`, or `/D`. Configure a new results root for a new analysis
comparison. Same-name biological output files can be overwritten by repeated
steps; only the run log directory is always unique. Downstream steps scan their
configured results directories, so stale tables should not be mixed between
cohorts.

## First-run acceptance

For the GSE104987 example, check the observed input metadata, inspect A's UMAP,
verify B's CSV schemas and effect-size direction, and inspect C's adjusted-p-value
method and plotting gates. Then enable optional steps individually. A preflight
pass is an input/dependency check, not numerical validation. Record any package
compatibility or scientific-method changes made during HPC testing before using
results in a manuscript.
