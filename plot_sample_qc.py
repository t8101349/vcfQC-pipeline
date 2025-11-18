import pandas as pd
import matplotlib.pyplot as plt
import re
import sys
import numpy as np

# Usage: python plot_sample_qc.py stats_file.stats

if len(sys.argv) < 2:
    print("Usage: python plot_sample_qc.py stats_file.stats")
    sys.exit(1)

stats_file = sys.argv[1]
samples = []

# 第一遍讀 PSC 行收集基本 info
with open(stats_file) as f:
    for line in f:
        if line.startswith("PSC"):
            fields = line.strip().split("\t")
            nTransversions = int(fields[7])
            samples.append({
                "sample": fields[2],
                "nRefHom": int(fields[3]),
                "nNonRefHom": int(fields[4]),
                "nHets": int(fields[5]),
                "nTransitions": int(fields[6]),
                "nTransversions": nTransversions,
                "nIndels": int(fields[8]),
                "avg_depth": float(fields[9]),
                "nSingletons": int(fields[10]),
                "nHapRef": int(fields[11]),
                "nHapAlt": int(fields[12]),
                "nMissing": int(fields[13]),
                "titv": float(fields[6])/nTransversions if nTransversions != 0 else float('nan')
            })

# 第二遍讀 TSTV 行更新 titv
with open(stats_file) as f:
    for line in f:
        if line.startswith("TSTV"):
            _, sample, ts, tv, ts_tv_ratio = line.strip().split()[:5]
            ts, tv, ts_tv_ratio = int(ts), int(tv), float(ts_tv_ratio)
            for s in samples:
                if s["sample"] == sample:
                    s["ts"] = ts
                    s["tv"] = tv
                    s["titv"] = ts_tv_ratio


df = pd.DataFrame(samples)
# 計算 missing rate
df["missing_rate"] = df["nMissing"] / (df["nRefHom"] + df["nNonRefHom"] + df["nHets"] + df["nMissing"])
# 計算 nVariants
df["nVariants"] = df["nRefHom"] + df["nNonRefHom"] + df["nHets"]
# 計算雜合率
df["het_rate"] = df["nHets"] / df["nVariants"]

# === Save to Excel ===
df.to_excel("sample_qc_metrics.xlsx", index=False)

excel_writer = pd.ExcelWriter("sample_qc_summary.xlsx", engine="xlsxwriter")
# === Plotting ===
plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(df["missing_rate"], bins=30, color="skyblue", edgecolor="black")
#plt.hist(df["missing_rate"], bins=30, color="skyblue", edgecolor="black")
plt.xlabel("Missing rate per sample")
plt.ylabel("Number of samples")
plt.title("Distribution of Missing Rate")
plt.savefig("sample_missing_rate.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="missing_rate", index=False)

plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(df["titv"].dropna(), bins=30, color="orange", edgecolor="black")
plt.xlabel("Ti/Tv ratio")
plt.ylabel("Number of samples")
plt.title("Distribution of Ti/Tv ratio")
plt.savefig("sample_titv.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="titv", index=False)

plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(df["nVariants"], bins=30, color="green", edgecolor="black")
plt.xlabel("Number of variants")
plt.ylabel("Number of samples")
plt.title("Distribution of Variants per Sample")
plt.savefig("sample_variants.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="nVariants", index=False)

plt.figure(figsize=(6,4))
counts, bins, _ = plt.hist(df["het_rate"], bins=30, color="purple", edgecolor="black")
plt.xlabel("Heterozygosity rate")
plt.ylabel("Number of samples")
plt.title("Distribution of Heterozygosity Rate")
plt.savefig("sample_het_rate.png", dpi=300)
plt.close()
pd.DataFrame({"bin_start": bins[:-1], "bin_end": bins[1:], "count": counts.astype(int)}).to_excel(excel_writer, sheet_name="het_rate", index=False)

# missing_rate
# 定義區間
bins = [0, 0.005, 0.01, 0.02, 0.05]
labels = [
    "0–0.005",
    "0.005–0.01",
    "0.01–0.02",
    "0.02–0.05",
]

# 分箱並計數
df["missing_rate_bin"] = pd.cut(df["missing_rate"], bins=bins, labels=labels, right=False)
counts = df["missing_rate_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="missing_rate_bin", index=False)

# titv
# 定義區間
bins = [2.07, 2.08, 2.085, 2.09, 2.10]
labels = [
    "2.07–2.08",
    "2.08–2.085",
    "2.085–2.09",
    "2.09–2.10",
]

# 分箱並計數
df["titv_bin"] = pd.cut(df["titv"], bins=bins, labels=labels, right=False)
counts = df["titv_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="titv_bin", index=False)

# nVariants
# 定義區間
bins = [40800000, 41200000, 41600000, 42000000, 42400000, 42800000, 43200000, 43600000]
labels = [
    "40800000–41200000",
    "41200000–41600000",
    "41600000–42000000",
    "42000000–42400000",
    "42400000–42800000",
    "42800000–43200000",
    "43200000–43600000",
]

# 分箱並計數
df["nVariants_bin"] = pd.cut(df["nVariants"], bins=bins, labels=labels, right=False)
counts = df["nVariants_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="nVariants_bin", index=False)

# heterozygosity_rate
# 定義區間
bins = [0.04, 0.042, 0.044, 0.046, 0.048]
labels = [
    "0.040–0.042",
    "0.042–0.044",
    "0.044–0.046",
    "0.046–0.048",
]

# 分箱並計數
df["het_rate_bin"] = pd.cut(df["het_rate"], bins=bins, labels=labels, right=False)
counts = df["het_rate_bin"].value_counts().reindex(labels, fill_value=0)

# 存到 Excel
pd.DataFrame({
    "range": labels,
    "count": counts.values
}).to_excel(excel_writer, sheet_name="het_rate_bin", index=False)


excel_writer.close()

print("✅ Done! Results saved to sample_qc_metrics.xlsx, sample_qc_summary.xlsx and PNG figures.")
