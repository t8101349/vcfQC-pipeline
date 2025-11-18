import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 檔案路徑
qc_file = "sample_qc_metrics_sex.csv"

# 讀取 QC 結果
df_qc = pd.read_csv(qc_file, sep=",")
df = df_qc[["sample", "nSingletons"]]

# 設定 bins (5000 一格，覆蓋所有數值範圍)
bin_size = 5000
max_value = df["nSingletons"].max()
bins = np.arange(0, max_value + bin_size, bin_size)

# 畫直方圖
plt.figure(figsize=(6,4))
counts, bins, patches = plt.hist(
    df["nSingletons"],
    bins=bins,
    edgecolor="black",
    color="skyblue"
)

plt.xlabel("nSingletons")
plt.ylabel("Number of Samples")
plt.title("Distribution of Singletons per Sample")
plt.tight_layout()
plt.savefig("nSingleton_hist.png", dpi=300)
plt.close()

# 匯出統計結果
hist_df = pd.DataFrame({
    "bin_start": bins[:-1],
    "bin_end": bins[1:],
    "count": counts.astype(int)
})
hist_df.to_csv("nSingleton_bins.csv", index=False)

print("✅ Done! 已輸出 nSingleton_hist.png 與 nSingleton_bins.csv")
