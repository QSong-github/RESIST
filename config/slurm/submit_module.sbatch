#!/usr/bin/env bash
#SBATCH --job-name=resist
#SBATCH --output=resist-%j.out
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Initial resource placeholders; adjust for the cohort and step.
# Activate the project environment before sbatch, or add your site-specific setup.
set -euo pipefail
: "${RESIST_HOME:?Set RESIST_HOME to the extracted RESIST package root.}"
[[ $# -ge 1 ]] || { echo "Usage: sbatch ... submit_module.sbatch A|B|C|D [runner options]" >&2; exit 2; }
module_id="$1"
shift
case "$module_id" in A|B|C|D) ;; *) echo "Module must be A, B, C, or D" >&2; exit 2 ;; esac
cd "$RESIST_HOME"
command -v Rscript >/dev/null
command -v python3 >/dev/null
exec bash "$RESIST_HOME/run_${module_id}.sh" "$@"
