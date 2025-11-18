import pandas as pd
import re
import sys
import numpy as np

# Usage: python platform_mix.py sample_qc_metrics.xlsx /CMUH_server/DATA7/SNParray/TWB2023/4.lab_info.csv

if len(sys.argv) < 3:
    print("Usage: python platform_mix.py sample_qc_metrics.xlsx /CMUH_server/DATA7/SNParray/TWB2023/4.lab_info.csv")
    sys.exit(1)

platform_file = sys.argv[2]
sample_file = sys.argv[1]

# 讀檔
lab_info = pd.read_csv(platform_file)
sample_qc = pd.read_excel(sample_file)
lab_info = lab_info[["Sample_Name","Platform"]]

# 對應 Sample_Name (lab_info) 與 sample (sample_qc)
merged = sample_qc.merge(
    lab_info,
    left_on="sample", right_on="Sample_Name",
    how="left"
)
merged = merged.drop(columns=["Sample_Name"]).rename(columns={"Platform": "platform"})

out_file = "sample_qc_metrics_platform.csv"
merged.to_csv(out_file, index=False)

print(f"✅ 合併完成，輸出檔案：{out_file}")