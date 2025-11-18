import matplotlib.pyplot as plt
import pandas as pd

# 檔案路徑
qc_file = "sample_qc_metrics.xlsx"
sex_file = "/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam"

# 讀取 QC 結果
df_qc = pd.read_excel(qc_file)

# 讀取 psam，注意第一欄 '#IID'
df_sex = pd.read_csv(sex_file, delim_whitespace=True)

# 改欄位名稱，方便 merge
df_sex = df_sex.rename(columns={"#IID": "sample", "SEX": "sex"})

# 合併
merged = df_qc.merge(df_sex, on="sample", how="left")

# 把數字 sex 轉成文字
sex_map = {1: "Male", 2: "Female", 0: "Unknown"}
merged["sex"] = merged["sex"].map(sex_map).fillna("NA")

# 輸出
out_file = "sample_qc_metrics_sex.csv"
merged.to_csv(out_file, index=False)

print(f"✅ 已合併完成，輸出檔案：{out_file}")

df = merged

# 假設 df 已經有 'nVariants' 和 'sex' 欄位
plt.figure(figsize=(6,4))

# 拆分不同性別
sex_groups = df['sex'].dropna().unique()
data = [df.loc[df['sex'] == s, 'nVariants'] for s in sex_groups]

# 疊加直方圖 (不同性別顏色不同)
counts, bins, patches = plt.hist(
    data,
    bins=30,
    stacked=True,
    edgecolor="black",
    label=sex_groups
)

plt.xlabel("Number of variants")
plt.ylabel("Number of samples")
plt.title("Distribution of Variants per Sample by Sex")
plt.legend(title="Sex")
plt.tight_layout()
plt.savefig("sample_variants_by_sex.png", dpi=300)
plt.close()

# 輸出每個 bin 各性別的計數
# hist_df = pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:]})
# for i, s in enumerate(sex_groups):
#     hist_df[s] = counts[i].astype(int)

# hist_df.to_excel(excel_writer, sheet_name="nVariants_by_sex", index=False)
