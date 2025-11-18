#!/bin/bash
# Usage: ./run_vcfqc.sh autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz  vcfqc_final_merged.vcf.gz

#合併
bcftools concat \
  -Oz -o vcfqc_final_merged.vcf.gz \
  autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz

tabix -p vcf vcfqc_final_merged.vcf.gz
bcftools sort vcfqc_final_merged.vcf.gz -Oz -o vcfqc_final_merged.sorted.vcf.gz
tabix -p vcf vcfqc_final_merged.sorted.vcf.gz


echo "==> Running bcftools stats..."
bcftools stats -s - vcfqc_final_merged.vcf.gz > merged.stats
echo "==> Generating plots with plot-vcfstats..."
plot-vcfstats -p merged_plots merged.stats

echo "==> Done!"
