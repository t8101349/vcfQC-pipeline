#!/bin/bash

# Usage: bash extract_snp_info.sh autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz combine_vcfqc_3_snp_info.tsv

autosome_vcf=$1
par_vcf=$2
nonpar_vcf=$3
output_tsv=$4

# 輸出表頭
echo -e "CHROM\tPOS\tID\tREF\tALT\tAC\tAF\tHWE\tF_MISSING\tREGION" > ${output_tsv}

# 定義 function，加入 REGION 標籤
extract_info () {
  local vcf=$1
  local region=$2
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AF\t%INFO/HWE\t%INFO/F_MISSING\t${region}\n" \
    ${vcf} >> ${output_tsv}
}

# 依序處理三個 VCF
extract_info ${autosome_vcf} "autosome"
extract_info ${par_vcf} "chrX_PAR"
extract_info ${nonpar_vcf} "chrX_nonPAR"
