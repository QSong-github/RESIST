#!/usr/bin/env python3
# =============================================================================
# RESIST | Module C - Regulatory mechanisms
# C4_select_significant_variants.py
# -----------------------------------------------------------------------------
# Purpose  : Reduce the Fisher results of C3 to the loci worth annotating - nominally
#            significant, with a finite non-zero odds ratio - and emit the genomic
#            coordinates in the two-column form the annotation step expects.
# Inputs   : fisher_results.tsv (from C3)
# Outputs  : input.csv - columns Chr, Pos
# Usage    : python C4_select_significant_variants.py --input fisher_results.tsv --output input.csv
# Origin   : Variant_detection/format.py
# =============================================================================

import pandas as pd
import numpy as np
import argparse

def filter_and_extract(input_file, output_file):
    # Read input TSV file
    df = pd.read_csv(input_file, sep="\t")

    # Apply filtering conditions
    df_filtered = df[
        (df["p_value"] < 0.05) &
        (df["odds_ratio"].notna()) &
        (df["odds_ratio"] != 0) &
        (np.isfinite(df["odds_ratio"]))
    ]

    # Extract chromosome and position (integer format)
    df_filtered[["Chr", "Pos"]] = df_filtered["ID"].str.extract(r"chr(\d+|X|Y)_(\d+)")
    df_filtered["Pos"] = df_filtered["Pos"].astype(int)  # Ensure Pos is integer

    # Save as CSV (avoid scientific notation)
    df_filtered[["Chr", "Pos"]].to_csv(output_file, index=False)
    print(f"Output completed: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract Chr and Pos based on filtering criteria")
    parser.add_argument("--input", required=True, help="Input file path (TSV)")
    parser.add_argument("--output", required=True, help="Output file path (CSV)")
    args = parser.parse_args()

    filter_and_extract(args.input, args.output)

