import pandas as pd
import numpy as np
import argparse

def filter_and_extract(input_file, output_file):
    # 读取输入文件
    df = pd.read_csv(input_file, sep="\t")

    # 筛选
    df_filtered = df[
        (df["p_value"] < 0.05) &
        (df["odds_ratio"].notna()) &
        (df["odds_ratio"] != 0) &
        (np.isfinite(df["odds_ratio"]))
    ]

    # 提取 Chr 和 Pos（整数形式）
    df_filtered[["Chr", "Pos"]] = df_filtered["ID"].str.extract(r"chr(\d+|X|Y)_(\d+)")
    df_filtered["Pos"] = df_filtered["Pos"].astype(int)  # 保证 Pos 是整数

    # 保存为 CSV（避免科学记数法）
    df_filtered[["Chr", "Pos"]].to_csv(output_file, index=False)
    print(f" 输出完成：{output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="提取满足条件的Chr和Pos")
    parser.add_argument("--input", required=True, help="输入文件路径（TSV）")
    parser.add_argument("--output", required=True, help="输出文件路径（CSV）")
    args = parser.parse_args()

    filter_and_extract(args.input, args.output)
