import re

file = "shapeit5/shapeit_chrx/complete/diff_columns.txt"

rows_start_dot = []
count_00 = 0
count_11 = 0
count_change = 0

with open(file, "r") as f:
    for line in f:
        line = line.strip()

        # 跳過空行
        if not line:
            continue

        # 是否以 "." 開頭
        if line.startswith("."):
            rows_start_dot.append(line)

            # 檢查結尾是否為 0|0 或 1|1
            if re.search(r"0\|0$", line):
                count_00 += 1
            elif re.search(r"1\|1$", line):
                count_11 += 1

        if line.startswith("1|1") or line.startswith("0|0"):
            # 檢查結尾是否改變
            if re.search(r"0\|1$", line) or re.search(r"1\|0$", line):
                count_change += 1

# 顯示結果
print("=== 以 '.' 開頭的行 ===")
for r in rows_start_dot:
    print(r)

print("\n=== 統計 ===")
print(f"缺值改為 0|0：{count_00}")
print(f"缺值改為 1|1：{count_11}")
print(f"count_change:{count_change}")

"""
with open("dot_rows_filtered.txt", "w") as f:
    for r in rows_start_dot:
        f.write(r + "\n")

print("\n已輸出 dot_rows_filtered.txt")
"""