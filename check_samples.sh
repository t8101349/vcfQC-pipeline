#!/bin/bash
# check_samples.sh autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <autosome.vcf.gz> <chrX_PAR.vcf.gz> <chrX_nonPAR.vcf.gz>"
    exit 1
fi

autosome=$1
par=$2
nonpar=$3

echo "==> Extracting sample lists..."
bcftools query -l ${autosome} > autosome.samples.txt
bcftools query -l ${par} > par.samples.txt
bcftools query -l ${nonpar} > nonpar.samples.txt

echo "==> Comparing sample lists..."
diff -q autosome.samples.txt par.samples.txt || echo "❌ autosome vs PAR 不一致"
diff -q autosome.samples.txt nonpar.samples.txt || echo "❌ autosome vs nonPAR 不一致"
diff -q par.samples.txt nonpar.samples.txt || echo "❌ PAR vs nonPAR 不一致"

if cmp -s autosome.samples.txt par.samples.txt && cmp -s autosome.samples.txt nonpar.samples.txt; then
    echo "✅ 三個 VCF 的 sample 名單完全一致，可以直接用 bcftools concat"
else
    echo "⚠️  sample 名單不一致，請用 bcftools merge"
fi
