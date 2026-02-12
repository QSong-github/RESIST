#!/bin/bash
set -euo pipefail

META_CSV="meta_condition_clean_4.csv"
REF_VCF="/home/liangjialu/blue_qsong1/liangjialu/data/variant-test/reference/Homo_sapiens_assembly38.dbsnp138.vcf"

echo "===== Step 1: 从 CSV 生成 sensitive/resistant barcodes ====="

rm -f sensitive_barcodes.txt resistant_barcodes.txt

awk -F',' '
NR>1 {
    cb=$3;                              # 第 3 列是 CB
    gsub(/[\r ]/, "", cb);              # 清理
    gsub("-1$", "", cb);                # 保险：去掉 -1（避免匹配错误）

    if ($2=="sensitive")  print cb >> "sensitive_barcodes.txt";
    if ($2=="resistant")  print cb >> "resistant_barcodes.txt";
}
' "$META_CSV"

echo "Sensitive barcodes:  $(wc -l < sensitive_barcodes.txt)"
echo "Resistant barcodes:  $(wc -l < resistant_barcodes.txt)"

echo "Files generated:"
ls -lh sensitive_barcodes.txt resistant_barcodes.txt


echo "===== Step 2: 使用 sensitive_merged.sorted.bam / resistant_merged.sorted.bam 跑 SNP ====="

mkdir -p vcf_sensitive vcf_resistant

echo "--- Running cellSNP for sensitive group ---"
cellsnp-lite \
    -s sensitive_merged.sorted.bam \
    -b sensitive_barcodes.txt \
    -O vcf_sensitive \
    -R "$REF_VCF" \
    --cellTAG CB \
    --UMItag UB \
    --gzip \
    --genotype \
    --minMAF 0.1 \
    --minCOUNT 20 \
    -p 8

echo "--- Running cellSNP for resistant group ---"
cellsnp-lite \
    -s resistant_merged.sorted.bam \
    -b resistant_barcodes.txt \
    -O vcf_resistant \
    -R "$REF_VCF" \
    --cellTAG CB \
    --UMItag UB \
    --gzip \
    --genotype \
    --minMAF 0.1 \
    --minCOUNT 20 \
    -p 8

echo "===== 完成啦！！ SNP calling 全部结束 🎉 ====="
echo "输出： vcf_sensitive/, vcf_resistant/"
