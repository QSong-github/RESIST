# Column descriptions:
# Chr: Chromosome number (input)
# Pos: Genomic position (input)
# rsid: rs identifier of the variation (unique identifier)
# Alleles: The type of alleles (e.g., C>G,T indicates reference allele C with variant alleles G and T)
# Chromosome (GRCh38): Chromosomal location based on GRCh38 genome build (format: ChromosomeNumber:Coordinate (BuildVersion))
# Chromosome (GRCh37): Chromosomal location based on GRCh37 genome build (format: ChromosomeNumber:Coordinate (BuildVersion))
# Gene: Gene name associated with the variation
# Functional Consequence: Functional impact of the variation (e.g., missense variant)
# Clinical significance (SNP): Clinical significance assessment from dbSNP database
# Identifiers (ClinVar): Variation identifiers from ClinVar (including description, ID, and accession)
# Type and length: Type and length of the variation (e.g., single nucleotide variant, 1bp)
# Protein change: Amino acid change at protein level (e.g., A284T indicates alanine to threonine at position 284)
# Molecular consequence (ClinVar): Molecular consequence of the variation recorded in ClinVar
# Gene (ClinVar): Gene name associated with the variation in ClinVar
# HI score: Haploinsufficiency score (gene dosage sensitivity assessment from ClinGen)
# TS score: Triplosensitivity score (gene dosage sensitivity assessment from ClinGen)
# Condition (Germline): Germline diseases/phenotypes associated with the variation
# Classification (Germline): Clinical classification of the variation for specific diseases in ClinVar
# Usage: python mutation_annotation.py input.csv output.csv

import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse, parse_qs
import pandas as pd
import re
import sys
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# ------------------------------
# Step 1: Get rsID from Chr and Pos
# ------------------------------
def get_rsid_from_chr_pos(chrom, pos):
    url = f"https://www.ncbi.nlm.nih.gov/snp/?term={pos}%5BPOSITION%5D+AND+{chrom}%5BCHR%5D"
    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
    except Exception as e:
        return "N/A"

    try:
        title_link = soup.find("p", class_="title").find("a", href=lambda h: h and "/snp/rs" in h)
        if title_link and title_link.string and title_link.string.startswith("rs"):
            return title_link.string.strip()
        
        rsid_span = soup.find("span", class_="rsid")
        if rsid_span and rsid_span.string and rsid_span.string.startswith("rs"):
            return rsid_span.string.strip()
        
        page_text = soup.get_text()
        match = re.search(r"rs\d+", page_text)
        if match:
            return match.group()
    except:
        pass
    
    return "N/A"

# ------------------------------
# Step 2: Get SNP info from rsID
# ------------------------------
def get_snp_info(rsid):
    url = f'https://www.ncbi.nlm.nih.gov/snp/?term={rsid}'
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
    except:
        return None

    result = {
        "rsid": rsid,
        "Alleles": "N/A",
        "Chromosome (GRCh38)": "N/A",
        "Chromosome (GRCh37)": "N/A",
        "Gene": "N/A",
        "Functional Consequence": "N/A",
        "Clinical significance (SNP)": "N/A"
    }

    try:
        if (alleles_dt := soup.find("dt", string=lambda text: text and "Alleles" in text.strip())) and (alleles_dd := alleles_dt.find_next("dd")):
            result["Alleles"] = alleles_dd.get_text().split("[")[0].strip()
    except:
        pass

    try:
        if (chrom_dt := soup.find("dt", string=lambda text: text and "Chromosome" in text.strip())) and (chrom_dd := chrom_dt.find_next("dd")):
            chrom_text = re.sub(r'\s+', ' ', chrom_dd.get_text()).strip()
            if grch38 := re.search(r'(\d+:\d+)\s*\(GRCh38\)', chrom_text):
                result["Chromosome (GRCh38)"] = f"{grch38.group(1)} (GRCh38)"
            if grch37 := re.search(r'(\d+:\d+)\s*\(GRCh37\)', chrom_text):
                result["Chromosome (GRCh37)"] = f"{grch37.group(1)} (GRCh37)"
    except:
        pass

    try:
        if (gene_dt := soup.find("dt", string=lambda text: text and "Gene" in text.strip())) and (gene_dd := gene_dt.find_next("dd")):
            result["Gene"] = gene_dd.get_text().split("(")[0].strip()
    except:
        pass

    try:
        if (func_dt := soup.find("dt", string=lambda text: text and "Functional Consequence" in text.strip())) and (func_dd := func_dt.find_next("dd")):
            result["Functional Consequence"] = func_dd.get_text().strip()
    except:
        pass

    try:
        clin_dt = soup.find("dt", string=lambda text: text and "Clinical significance" in text.strip())
        if not clin_dt:
            clin_a = soup.find("a", string=lambda text: text and "Clinical significance" in text.strip())
            if clin_a:
                clin_dt = clin_a.find_parent("dt")
        if clin_dt and (clin_dd := clin_dt.find_next("dd")):
            result["Clinical significance (SNP)"] = clin_dd.get_text().strip()
    except:
        pass

    return result

# ------------------------------
# Step 3: Get ClinVar info from rsID
# ------------------------------
def get_clinvar_info(rsid):
    url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, "html.parser")
    except:
        return None

    results = []
    identifiers = "N/A"
    try:
        if (identifiers_dd := soup.find("dt", string="Identifiers").find_next("dd")):
            var_desc = identifiers_dd.find("p").get_text(strip=True) if identifiers_dd.find("p") else "N/A"
            var_id = [line for line in identifiers_dd.stripped_strings if "Variation ID" in line][0].split(":")[1].strip() if identifiers_dd else "N/A"
            accession = [line for line in identifiers_dd.stripped_strings if "Accession" in line][0].split(":")[1].strip() if identifiers_dd else "N/A"
            identifiers = f"{var_desc}; Variation ID: {var_id}; Accession: {accession}"
    except:
        pass

    type_length = "N/A"
    try:
        if (tl_dd := soup.find("dt", string="Type and length").find_next("dd")):
            type_length = tl_dd.get_text(strip=True)
    except:
        pass

    protein_change = "N/A"
    try:
        if (pc_dd := soup.find("dt", string="Protein change").find_next("dd")):
            protein_change = pc_dd.get_text(strip=True)
    except:
        pass

    molecular_consequence = "N/A"
    try:
        if (hgvs_table := soup.find("table", class_="hgvstable")):
            if (missense_td := hgvs_table.find("td", string=lambda text: text and "missense" in text.lower())):
                molecular_consequence = missense_td.get_text(strip=True)
    except:
        pass

    gene = "N/A"
    try:
        if (variation_links := soup.find_all('a', {'data-ga-label': 'variation viewer for gene'})):
            for link in variation_links:
                if (href := link.get('href', '')) and (query_params := parse_qs(urlparse(href).query)) and 'q' in query_params:
                    gene = query_params['q'][0].strip()
                    if gene:
                        break
    except:
        pass

    hi_score = "N/A"
    try:
        if (hi_td := soup.find("th", string="HI score").find_next("td")):
            hi_score = hi_td.get_text(strip=True)
    except:
        pass

    ts_score = "N/A"
    try:
        if (ts_td := soup.find("th", string="TS score").find_next("td")):
            ts_score = ts_td.get_text(strip=True)
    except:
        pass

    try:
        if (conditions_table := soup.find("table", class_="conditions-germline-list")):
            for row in conditions_table.find("tbody").find_all("tr"):
                condition = "N/A"
                if (condition_td := row.find("td", class_="interpreted-conditions")):
                    condition = condition_td.get_text(strip=True)
                
                classification = "N/A"
                if len(row.find_all("td")) > 1 and (class_td := row.find_all("td")[1]):
                    classification = class_td.get_text(strip=True)
                
                results.append({
                    "rsid": rsid,
                    "Identifiers (ClinVar)": identifiers,
                    "Type and length": type_length,
                    "Protein change": protein_change,
                    "Molecular consequence (ClinVar)": molecular_consequence,
                    "Gene (ClinVar)": gene,
                    "HI score": hi_score,
                    "TS score": ts_score,
                    "Condition (Germline)": condition,
                    "Classification (Germline)": classification
                })
        else:
            results.append({
                "rsid": rsid,
                "Identifiers (ClinVar)": identifiers,
                "Type and length": type_length,
                "Protein change": protein_change,
                "Molecular consequence (ClinVar)": molecular_consequence,
                "Gene (ClinVar)": gene,
                "HI score": hi_score,
                "TS score": ts_score,
                "Condition (Germline)": "N/A",
                "Classification (Germline)": "N/A"
            })
    except:
        results.append({
            "rsid": rsid,
            "Identifiers (ClinVar)": identifiers,
            "Type and length": type_length,
            "Protein change": protein_change,
            "Molecular consequence (ClinVar)": molecular_consequence,
            "Gene (ClinVar)": gene,
            "HI score": hi_score,
            "TS score": ts_score,
            "Condition (Germline)": "N/A",
            "Classification (Germline)": "N/A"
        })

    return results

# ------------------------------
# Single record processing
# ------------------------------
def process_record(chrom, pos):
    results = []
    
    # Get rsID
    rsid = get_rsid_from_chr_pos(chrom, pos)
    
    # Prepare base data with input coordinates
    base_data = {"Chr": chrom, "Pos": pos, "rsid": rsid}
    
    # If no rsID found
    if rsid == "N/A":
        return [base_data]
    
    # Get SNP info
    snp_data = get_snp_info(rsid)
    if not snp_data:
        return [base_data]
    
    # Get ClinVar info
    clinvar_data_list = get_clinvar_info(rsid)
    if not clinvar_data_list:
        merged = {**base_data,** snp_data}
        return [merged]
    
    # Merge all data
    for clinvar_data in clinvar_data_list:
        merged = {**base_data,** snp_data, **clinvar_data}
        results.append(merged)
    
    return results

# ------------------------------
# Main workflow with parallel processing
# ------------------------------
def main():
    if len(sys.argv) != 3:
        print("Usage: python get_all_parallel.py input_chr_pos.csv output_annotations.csv", file=sys.stderr)
        print("Input file must contain 'Chr' and 'Pos' columns", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Define column order
    cols = ["Chr", "Pos", "rsid", "Alleles", "Chromosome (GRCh38)", "Chromosome (GRCh37)", 
            "Gene", "Functional Consequence", "Clinical significance (SNP)", 
            "Identifiers (ClinVar)", "Type and length", "Protein change", 
            "Molecular consequence (ClinVar)", "Gene (ClinVar)", "HI score", 
            "TS score", "Condition (Germline)", "Classification (Germline)"]

    # Read input data
    try:
        df = pd.read_csv(input_file, sep=None, engine="python")
        if "Chr" not in df.columns or "Pos" not in df.columns:
            print("Input file must contain 'Chr' and 'Pos' columns", file=sys.stderr)
            sys.exit(1)
        records = list(zip(df["Chr"], df["Pos"]))
    except Exception as e:
        print(f"Failed to read input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Process in parallel with progress bar
    all_results = []
    max_workers = 5  # Adjust based on your network capacity
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        futures = {executor.submit(process_record, chrom, pos): (chrom, pos) for chrom, pos in records}
        
        # Process results as they complete
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing records"):
            chrom, pos = futures[future]
            try:
                results = future.result()
                all_results.extend(results)
            except Exception as e:
                print(f"Error processing Chr{chrom}:{pos}: {e}", file=sys.stderr)
                all_results.append({"Chr": chrom, "Pos": pos, "rsid": "Error"})

    # Write output
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=cols, restval="N/A")
        writer.writeheader()
        for result in all_results:
            # Ensure all columns are present
            row = {col: result.get(col, "N/A") for col in cols}
            writer.writerow(row)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()

