
input_vcf=/CMUH_server/DATA11/Research/dragen_run_joint/VCF/TWBK_1484Cases/TWBK_1484Cases.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $input_vcf \
  | awk -F'\t' '{ n=split($4, a, ","); if(n>1) print $1,$2,n }' OFS='\t' > multiallelic_sites.txt


bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $input_vcf \
  | awk -F'\t' '{ n=split($4,a,","); if(n>1) print n }' \
  | sort | uniq -c

# wc -l multiallelic_sites.txt
# awk '{sum += $3} END {print sum}' multiallelic_sites.txt
