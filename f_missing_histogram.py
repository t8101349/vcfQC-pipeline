import gzip
import matplotlib.pyplot as plt
import numpy as np

vcf_file = "/home/Weber/vcfQC/TWB_Drogen_vcfqc_2-2.vcf.gz"
f_missing = []

# 讀取 VCF，抓 INFO 欄位的 F_MISSING
with gzip.open(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        info = fields[7]  # INFO 欄位
        # 找 F_MISSING
        for entry in info.split(';'):
            if entry.startswith('F_MISSING='):
                try:
                    f_missing.append(float(entry.split('=')[1]))
                except:
                    pass
                break

print(f"總 SNP 數: {len(f_missing)}")

# 統計不同 cutoff
cutoffs = [0.1, 0.05, 0.01]
for c in cutoffs:
    n = sum([x < c for x in f_missing])
    print(f"Missing rate < {c}: {n} SNPs")

# 畫 histogram
plt.hist(f_missing, bins=50, color='skyblue', edgecolor='black')
for c, col in zip(cutoffs, ['red','orange','green']):
    plt.axvline(c, color=col, linestyle='--', linewidth=2)
plt.xlabel("Missing Rate (F_MISSING)")
plt.ylabel("Number of SNPs")
plt.title("SNP Missing Rate Distribution")
# 設定 x 軸 0.01 為一刻度
plt.xticks(np.arange(0, 1.01, 0.01))  # 從 0 到 1，每 0.01 一格
plt.legend([f"<{c}" for c in cutoffs])

plt.savefig("/home/Weber/vcfQC/f_missing_distribution.png", dpi=300)
plt.show()
print("Histogram 已存檔：/home/Weber/vcfQC/f_missing_distribution.png")