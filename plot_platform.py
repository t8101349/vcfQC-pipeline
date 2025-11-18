import matplotlib.pyplot as plt
import pandas as pd
import sys

# Usage: python plot_platform.py sample_qc_metrics_platform.csv

if len(sys.argv) < 2:
    print("Usage: python plot_platform.py sample_qc_metrics_platform.csv")
    sys.exit(1)

file = sys.argv[1]
df = pd.read_csv(file)

# excel_writer = pd.ExcelWriter("sample_qc_summary.xlsx", engine="xlsxwriter")

# 假設 df 已經包含 'nVariants' 與 'platform'
plt.figure(figsize=(6,4))

# 拆分不同平台資料
platforms = df['platform'].dropna().unique()
data = [df.loc[df['platform'] == p, 'nVariants'] for p in platforms]

# 疊加直方圖
counts, bins, patches = plt.hist(
    data,
    bins=30,
    stacked=True,
    edgecolor="black",
    label=platforms
)

plt.xlabel("Number of variants")
plt.ylabel("Number of samples")
plt.title("Distribution of Variants per Sample by Platform")
plt.legend(title="Platform")
plt.tight_layout()
plt.savefig("sample_variants_by_platform.png", dpi=300)
plt.close()

# 輸出每個 bin 各平台的計數
# hist_df = pd.DataFrame(
#     {"bin_start": bins[:-1], "bin_end": bins[1:]}
# )
# for i, p in enumerate(platforms):
#     hist_df[p] = counts[i].astype(int)

# hist_df.to_excel(excel_writer, sheet_name="nVariants_by_platform", index=False)
