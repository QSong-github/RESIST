#!/usr/bin/env python3
# =============================================================================
# RESIST | Module D - Immunogenomic features of resistance
# D1_detect_kit_version.py
# -----------------------------------------------------------------------------
# Purpose  : Read the Cell Ranger web summary of each sample and infer the 10x
#            chemistry (v2 or v3). scUTRquant needs the kit version to place poly(A)
#            sites correctly, so this runs before any APA quantification.
# Inputs   : list.txt - one sample directory or web_summary.html path per line
# Outputs  : kit.csv - input_path, kit_version
# Usage    : python D1_detect_kit_version.py --txt list.txt --out kit.csv
# Origin   : APA_analysis/find_kit_version.py
# Notes    : The original README showed `find_kit_version.py list.txt`; the
#            script has always required the `--txt` flag. Corrected here.
# =============================================================================

import re
import argparse
from pathlib import Path

def detect_10x_version(html_path: Path):
    """
    Extracts the Cell Ranger Chemistry version from web_summary.html.
    Returns '10xv2', '10xv3', or error status.
    """
    if not html_path.exists():
        return "FILE_NOT_FOUND"
    
    try:
        # Read the HTML content
        txt = html_path.read_text(errors="ignore")
        
        # Priority 1: Search for JSON structured data ["Chemistry","xxxxx"]
        m = re.search(r'\["Chemistry"\s*,\s*"([^"]+)"\]', txt)
        if m:
            chem = m.group(1)
        else:
            # Priority 2: Fallback search for 'Chemistry' string followed by quotes
            m2 = re.search(r'Chemistry[^"\']*["\']([^"\']+)["\']', txt)
            chem = m2.group(1) if m2 else "UNKNOWN"

        chem_low = chem.lower()
        
        if "v3" in chem_low:
            return "10xv3"
        if "v2" in chem_low:
            return "10xv2"
            
        return f"OTHER({chem})"
    except Exception as e:
        return f"ERROR_READING"

def main():
    parser = argparse.ArgumentParser(description="Extract 10x Chemistry version from a list of paths")
    parser.add_argument("--txt", required=True, help="Input txt file containing paths to web_summary.html or sample folders")
    parser.add_argument("--out", help="Optional output CSV file path")
    args = parser.parse_args()

    results = []
    
    # Load paths from input file
    with open(args.txt, 'r') as f:
        paths = [line.strip() for line in f if line.strip()]

    # Print header for CLI output
    print(f"{'Source_Path':<60} | {'Version':<12}")
    print("-" * 75)

    for p in paths:
        path_obj = Path(p)
        
        # Smart detection: if path is a directory, look for web_summary.html inside it or in outs/
        target_html = path_obj
        if path_obj.is_dir():
            if (path_obj / "web_summary.html").exists():
                target_html = path_obj / "web_summary.html"
            elif (path_obj / "outs" / "web_summary.html").exists():
                target_html = path_obj / "outs" / "web_summary.html"

        version = detect_10x_version(target_html)
        print(f"{str(path_obj)[-58:] if len(str(path_obj)) > 58 else str(path_obj):<60} | {version:<12}")
        results.append((p, version))

    # Write to CSV if --out is provided
    if args.out:
        with open(args.out, 'w') as f:
            f.write("input_path,kit_version\n")
            for p, v in results:
                f.write(f'"{p}","{v}"\n')
        print(f"\n[Done] Results saved to: {args.out}")

if __name__ == "__main__":
    main()