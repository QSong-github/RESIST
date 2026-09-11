#!/usr/bin/env python3
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C6_batch_annotate.py
# -----------------------------------------------------------------------------
# Purpose  : Walk a results tree and run C4 followed by C5 for every group directory
#            that already holds a fisher_results.tsv, skipping those that are done.
#            Convenience driver for cohorts with many sample groups.
# Inputs   : <tree>/**/result/group*/fisher_results.tsv
# Outputs  : <tree>/**/result/group*/input.csv and out.csv
# Usage    : python C6_batch_annotate.py --format_script C4_select_significant_variants.py \
#            --mutation_script C5_annotate_variants.py
# Origin   : Variant_detection/annotation.py
# Notes    : Scans the current working directory; cd into the results tree first.
# =============================================================================

import os
import argparse
import subprocess

def run_scripts(format_script, mutation_script):
    current_dir = os.getcwd()

    for root, dirs, files in os.walk(current_dir):
        if not root.endswith("result"):
            continue
        for d in dirs:
            if not d.startswith("group"):
                continue
            group_dir = os.path.join(root, d)
            fisher_file = os.path.join(group_dir, "fisher_results.tsv")
            input_csv = os.path.join(group_dir, "input.csv")
            output_csv = os.path.join(group_dir, "out.csv")

            if os.path.exists(output_csv):
                print(f"Skipping {output_csv}, already exists.")
                continue

            if os.path.exists(fisher_file):
                print(f"Processing: {fisher_file}")

                # Run format.py
                cmd1 = ["python", format_script, "--input", fisher_file, "--output", input_csv]
                subprocess.run(cmd1, check=True)

                # Run mutation_annotation.py
                cmd2 = ["python", mutation_script, input_csv, output_csv]
                subprocess.run(cmd2, check=True)

                print(f"Finished: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--format_script", required=True, help="Path to format.py script")
    parser.add_argument("--mutation_script", required=True, help="Path to mutation_annotation.py script")
    args = parser.parse_args()

    run_scripts(args.format_script, args.mutation_script)
