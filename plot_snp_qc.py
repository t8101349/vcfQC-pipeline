import pandas as pd
import matplotlib.pyplot as plt
import re
import sys
import numpy as np

# 使用方法: python plot_snp_qc.py vcfqc_3_snp_info.tsv

# 讀取輸入檔
if len(sys.argv) < 2:
    infile = "snp_info.tsv"
else:
    infile = sys.argv[1]

df = pd.read_csv(infile, sep="\t")
df["MAF"] = df["AF"].apply(lambda x: x if x <= 0.5 else 1 - x)


excel_writer = pd.ExcelWriter("snp_qc_summary.xlsx", engine="xlsxwriter")

# AF histogram
counts, bins, _ = plt.hist(df["MAF"], bins=50, color="salmon", edgecolor="black")
plt.xlabel("Allele Frequency (AF)")
plt.ylabel("Number of SNPs")
plt.title("Allele Frequency Distribution")
plt.savefig("AF_hist.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="MAF", index=False)


# HWE histogram (-log10)
plt.hist(df["HWE"], bins=50, color="lightblue", edgecolor="black", log=True)
plt.xlabel("-log10(HWE P-value)")
plt.ylabel("Number of SNPs")
plt.title("HWE P-value Distribution")
plt.savefig("HWE_hist.png", dpi=300)
plt.close()
counts, bins = np.histogram(df["HWE"], bins=50)
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts}).to_excel(excel_writer, sheet_name="HWE", index=False)

# (1b) HWE zoom-in (1e-6 ~ 0.02)
subset_hwe = df[(df["HWE"] >= 1e-6) & (df["HWE"] <= 0.02)]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset_hwe["HWE"], bins=30, color="teal", edgecolor="black")
plt.xlabel("HWE P-value (1e-6 ~ 0.02)")
plt.ylabel("Number of SNPs")
plt.title("HWE P-value Distribution (Zoom-in)")
plt.savefig("HWE_hist_1e-6_0.02.png", dpi=300)
plt.close()


# 存到 Excel
pd.DataFrame({
    "bin_start": bins[:-1],
    "bin_end": bins[1:], 
    "count": counts.astype(int)
}).to_excel(excel_writer, sheet_name="HWE_1e-6_0.02", index=False)

# (2b) HWE zoom-in (0.9 ~ 1)
subset_hwe = df[(df["HWE"] >= 0.9) & (df["HWE"] <= 1)]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset_hwe["HWE"], bins=30, color="teal", edgecolor="black")
plt.xlabel("HWE P-value (0.9 ~ 1.0)")
plt.ylabel("Number of SNPs")
plt.title("HWE P-value Distribution (Zoom-in)")
plt.savefig("HWE_hist_0.9_1.0.png", dpi=300)
plt.close()


# 存到 Excel
pd.DataFrame({
    "bin_start": bins[:-1],
    "bin_end": bins[1:], 
    "count": counts.astype(int)
}).to_excel(excel_writer, sheet_name="HWE_0.9_1.0", index=False)


# 定義區間
bins = [0, 1e-6, 5e-6, 1e-5, 5e-4, 0.001, 0.01, 0.05, 1.01]
labels = [
    "0–1E-6",
    "1E-6–5E-6",
    "5E-6–1E-5",
    "1E-5–5E-4",
    "5E-4–0.001",
    "0.001–0.01",
    "0.01–0.05",
    "0.05–1"
]

# 分箱並計數
df["HWE_bin"] = pd.cut(df["HWE"], bins=bins, labels=labels, right=False)
counts = df["HWE_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="HWE_bins", index=False)


# Missing rate histogram (0.0–0.1 區間)
subset = df[df["F_MISSING"] <= 0.1]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset["F_MISSING"], bins=20, color="lightgreen", edgecolor="black")
plt.xlabel("Missing rate (F_MISSING, 0.0–0.1)")
plt.ylabel("Number of SNPs")
plt.title("SNP Missing Rate Distribution (0.0–0.1)")
plt.savefig("F_MISSING_hist_0-0.1.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="F_MISSING_0-0.1", index=False)

# (2b) Missing rate zoom-in (0–0.01)
subset2 = df[df["F_MISSING"] <= 0.01]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset2["F_MISSING"], bins=20, color="orange", edgecolor="black")
plt.xlabel("Missing rate (F_MISSING, 0.0–0.01)")
plt.ylabel("Number of SNPs")
plt.title("SNP Missing Rate Distribution (0.0–0.01)")
plt.savefig("F_MISSING_hist_0-0.01.png", dpi=300)
plt.close()

pd.DataFrame({
    "bin_start": bins[:-1],
    "bin_end": bins[1:], 
    "count": counts.astype(int)
}).to_excel(excel_writer, sheet_name="F_MISSING_0-0.01", index=False)

# 定義區間
bins = [0, 0.001, 0.01, 0.05, 0.11]
labels = ["0-0.001", "0.001-0.01", "0.01-0.05", "0.05-0.1"]

# 分箱並計數
df["F_MISS_bin"] = pd.cut(df["F_MISSING"], bins=bins, labels=labels, right=False)
counts = df["F_MISS_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="F_MISSING_bins", index=False)


# F_MISSING vs HWE scatter
plt.scatter(df["F_MISSING"], df["HWE"], alpha=0.5, s=10)
plt.xlabel("F_MISSING")
plt.ylabel("HWE P-value")
plt.yscale("log")
plt.title("F_MISSING vs HWE")
plt.savefig("Missing_vs_HWE.png", dpi=300)
plt.close()


# (3b) MAF zoom-in (0–0.05)
subset_maf = df[df["MAF"] <= 0.05]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset_maf["MAF"], bins=20, color="green", edgecolor="black")
plt.xlabel("Allele Frequency (AF)")
plt.ylabel("Number of SNPs")
plt.title("Allele Frequency Distribution (Zoom-in)(0–0.05)")
plt.xticks(fontsize=8)
plt.savefig("MAF_hist_0-0.05.png", dpi=300)
plt.close()

subset_maf = df[df["MAF"] <= 0.002]
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(subset_maf["MAF"], bins=20, color="green", edgecolor="black")
plt.xlabel("Allele Frequency (AF)")
plt.ylabel("Number of SNPs")
plt.title("Allele Frequency Distribution (Zoom-in)(0–0.002)")
plt.xticks(fontsize=8)
plt.savefig("MAF_hist_0-0.002.png", dpi=300)
plt.close()

# 存 Excel (細分 bin 統計)
pd.DataFrame({
    "bin_start": bins[:-1],
    "bin_end": bins[1:], 
    "count": counts.astype(int)
}).to_excel(excel_writer, sheet_name="MAF_0-0.02", index=False)

# 定義 MAF 區間
maf_bins = [0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.51]
maf_labels = ["0–0.0001", "0.0001–0.001", "0.001–0.01", "0.01–0.05", "0.05–0.1", "0.1–0.5"]

# 分箱計數
df["MAF_bin_zoom"] = pd.cut(df["MAF"], bins=maf_bins, labels=maf_labels, right=False)
counts_zoom = df["MAF_bin_zoom"].value_counts().reindex(maf_labels, fill_value=0)

# 存 Excel (四區間統計)
pd.DataFrame({
    "range": maf_labels,
    "count": counts_zoom.values
}).to_excel(excel_writer, sheet_name="MAF_bins_0-0.1", index=False)


# 存 Excel
excel_writer.close()

print("Plots saved: AF_hist.png, HWE_hist.png, Missing_vs_HWE.png")
