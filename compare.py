import pandas as pd

# 讀取檔案
df1 = pd.read_csv("shapeit5/shapeit_chrx/complete/chrX_3_head.tsv", sep="\t")
df2 = pd.read_csv("shapeit5/shapeit_chrx/complete/chrX_4_head.tsv", sep="\t")

# 找出所有差異的列
diff_rows_mask = df1.ne(df2).any(axis=1)
df1_diff = df1[diff_rows_mask].reset_index(drop=True)
df2_diff = df2[diff_rows_mask].reset_index(drop=True)
comparison = df1_diff.ne(df2_diff)

# 建立標記 DataFrame，格式 (前)1/(後)1
diff_output = pd.DataFrame()
for col in df1_diff.columns:
    # 只保留有差異的欄位
    if comparison[col].any():
        diff_output[col] = df1_diff[col].astype(str) + " / " + df2_diff[col].astype(str)

# 輸出
diff_output.to_csv("shapeit5/shapeit_chrx/complete/diff_format.tsv", sep="\t", index=False)
print("已輸出 diff_format.tsv")

with open("shapeit5/shapeit_chrx/complete/diff_columns.txt", "w") as f:
    for col in df1_diff.columns:
        # 找出此欄位有差異的列
        col_diff_mask = df1_diff[col] != df2_diff[col]
        if col_diff_mask.any():
            for v1, v2 in zip(df1_diff[col][col_diff_mask], df2_diff[col][col_diff_mask]):
                f.write(f"{v1} / {v2}\n")
            f.write("\n")

print("已輸出 diff_columns.txt")