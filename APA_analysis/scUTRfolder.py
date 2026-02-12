#!/usr/bin/env python3
import os
import re
import csv
from pathlib import Path
import argparse


# -------------------------
# 精确解析 Chemistry → v2 / v3
# -------------------------
def detect_10x_version(html_path: Path):
    """从 web_summary.html 精确解析 Cell Ranger Chemistry."""
    txt = html_path.read_text(errors="ignore")

    # 优先解析 ["Chemistry","xxxxx"] 结构
    m = re.search(r'\["Chemistry"\s*,\s*"([^"]+)"\]', txt)
    if m:
        chem = m.group(1)
    else:
        # fallback：尝试解析 'Chemistry' 后的字段
        m2 = re.search(r'Chemistry[^"\']*["\']([^"\']+)["\']', txt)
        chem = m2.group(1) if m2 else ""

    chem_low = chem.lower()

    if "v3" in chem_low:
        return "10xv3"
    if "v2" in chem_low:
        return "10xv2"

    print(f"[WARN] 无法从 Chemistry 判断版本，默认 v3。Chemistry='{chem}'")
    return "10xv3"


# -------------------------
# 工具函数：提取 GSE / SRR
# -------------------------
def extract_gse(path: str):
    m = re.search(r"(GSE\d+)", path)
    return m.group(1) if m else "UNKNOWN"


def extract_srr(path: str):
    m = re.search(r"(SRR\d+)", path)
    return m.group(1) if m else "UNKNOWN"


# -------------------------
# 主样本处理函数
# -------------------------
def process_one(matrix_path: str, out_root: Path):
    matrix_path = Path(matrix_path).resolve()

    gse = extract_gse(str(matrix_path))
    srr = extract_srr(str(matrix_path))
    sample_id = f"{gse}-{srr}"     # 你要求的命名格式

    outs_dir = matrix_path.parent          # outs/
    sample_root = outs_dir.parent          # SRRxxxxxx/

    bam = sample_root / "outs" / "possorted_genome_bam.bam"
    annots = sample_root / "outs" / "annots.csv"
    summary = sample_root / "outs" / "web_summary.html"

    if not summary.exists():
        print(f"[ERR] Chemistry 文件缺失: {summary}, 跳过")
        return

    tech = detect_10x_version(summary)

    # 输出目录
    outdir = out_root / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    # -------------------------
    # 写 sample_sheet.csv
    # -------------------------
    sample_sheet_path = outdir / "sample_sheet.csv"
    with sample_sheet_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "file_type", "files"])
        w.writerow([sample_id, "bam", str(bam)])

    # -------------------------
    # 写 config.yaml
    # -------------------------
    config_text = f"""dataset_name: "{gse}"
sample_file: "{sample_sheet_path}"
cell_annots: "{annots}"
target: "utrome_hg38_v1"

tech: "{tech}"
strand: "--fr-stranded"

min_umis: 1000
"""

    (outdir / "config.yaml").write_text(config_text)

    print(f"[OK] {sample_id}: 写入 config.yaml & sample_sheet.csv")


# -------------------------
# 主程序入口
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="为 scUTRquant 批量生成输入配置")
    ap.add_argument("--txt", required=True, help="含 bc_matrix 路径的 txt 文件")
    ap.add_argument("--out", required=True, help="输出目录")
    args = ap.parse_args()

    out_root = Path(args.out).resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    paths = []
    with open(args.txt) as f:
        for line in f:
            line = line.strip()
            if line:
                paths.append(line)

    for p in paths:
        process_one(p, out_root)


if __name__ == "__main__":
    main()
