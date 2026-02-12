#!/bin/bash

# ===== 0. 输入：你的 parent_path =====
PARENT_DIR="/path/to/parent_path"      # 改成你的父目录路径
TRANSCRIPTOME="/home/liangjialu/blue_qsong1/liangjialu/data/reference/refdata-gex-GRCh38-2024-A"

log(){ echo "[`date '+%H:%M:%S'`] $1"; }

log "Scanning parent directory: $PARENT_DIR"
echo

# ===== 1. 遍历 parent_path 下所有 GSE* 目录 =====
for GSE_DIR in "$PARENT_DIR"/GSE*; do
    [ -d "$GSE_DIR" ] || continue
    GSE_NAME=$(basename "$GSE_DIR")

    log "=== Processing GSE project: $GSE_NAME ==="

    # 创建 result 目录，例如 GSE108397-result
    RESULT_DIR="${GSE_DIR}-result"
    mkdir -p "$RESULT_DIR"
    log "Created result directory: $RESULT_DIR"

    # test_fastqs 目录路径
    TEST_FASTQ_DIR="$GSE_DIR/test_fastqs"
    if [ ! -d "$TEST_FASTQ_DIR" ]; then
        log "  No test_fastqs folder found → SKIP $GSE_NAME"
        echo
        continue
    fi

    # ===== 2. 遍历 test_fastqs 下的每个 SRR 目录 =====
    for SRR_DIR in "$TEST_FASTQ_DIR"/SRR*; do
        [ -d "$SRR_DIR" ] || continue
        SRR=$(basename "$SRR_DIR")
        log "  → Found SRR: $SRR"

        # 获取 FASTQ 数量
        fastqs=( "$SRR_DIR"/*.fastq.gz )
        n=${#fastqs[@]}

        log "     FASTQ count = $n"

        # 默认使用原目录
        FASTQ_INPUT="$SRR_DIR"

        # =====================================
        # Case 1 — 正常 2 个 fastq，直接用
        # =====================================
        if (( n == 2 )); then
            log "     2 fastqs → normal 10x → use original directory"
        fi

        # =====================================
        # Case 2 — 3 个 fastq，需要修复
        # =====================================
        if (( n == 3 )); then
            log "     3 fastqs → fixing naming (_3 → R2, _2 → R1)"

            FIX_DIR="$RESULT_DIR/${SRR}-fixed"
            mkdir -p "$FIX_DIR"

            declare -A mapping
            for f in "${fastqs[@]}"; do
                bn=$(basename "$f")
                [[ "$bn" == *_3.fastq.gz ]] && mapping[R2]="$f"
                [[ "$bn" == *_2.fastq.gz ]] && mapping[R1]="$f"
            done

            if [ -z "${mapping[R1]}" ] || [ -z "${mapping[R2]}" ]; then
                log "     ERROR: Cannot map R1/R2 → SKIP"
                echo
                continue
            fi

            cp "${mapping[R1]}" "$FIX_DIR/${SRR}_R1.fastq.gz"
            cp "${mapping[R2]}" "$FIX_DIR/${SRR}_R2.fastq.gz"

            FASTQ_INPUT="$FIX_DIR"
            log "     Fixed fastqs stored in: $FIX_DIR"
        fi

        # =====================================
        # Case 3 — 非 2/3 个 fastq，跳过
        # =====================================
        if (( n != 2 && n != 3 )); then
            log "     $n fastqs → Not 10x (SMART-seq/inDrop?) → SKIP"
            echo
            continue
        fi

        # =====================================
        # 执行 cellranger count
        # =====================================
        log "     Running cellranger count for $SRR"

        (
          cd "$RESULT_DIR" || exit 1
          cellranger count \
              --id="$SRR" \
              --transcriptome="$TRANSCRIPTOME" \
              --fastqs="$FASTQ_INPUT" \
              --sample="$SRR" \
              --localcores=12 \
              --localmem=128 \
              --create-bam true
        )

        echo
    done
    echo

done

log "✔ ALL DONE."
