#!/usr/bin/env bash
set -euo pipefail
# sexcheck_fail_compute.sh autosome_singletonremoved.vcf.gz sex_check/chrX_sexcheck_failed_samples.txt

# ----- edit these -----
VCF=$1
SAMPLES=$2  # sex_check/chrX_sexcheck_failed_samples.txt
# VCF="autosome_singletonremoved.vcf.gz"   # your autosomal vcf (chr1-22)
# SAMPLE="NGS20140610G"                  # sample of interest
# ----------------------

REGIONS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

while read SAMPLE; do
  if [[ -z "$SAMPLE" ]]; then
    continue
  fi
  # 1) total autosomal variant sites
  total_sites=$(bcftools view -r ${REGIONS} ${VCF} | bcftools view -H | wc -l)

  # 2) called (non-missing) genotypes for SAMPLE on autosomes
  called=$(bcftools view -s "${SAMPLE}" -r ${REGIONS} ${VCF} | \
    bcftools query -f '[%GT\n]' | \
    grep -vE '^\./\.|^\.\|\.|^\.$' | wc -l)

  # 3) heterozygous genotype count for SAMPLE on autosomes
  hets=$(bcftools view -s "${SAMPLE}" -r ${REGIONS} ${VCF} | \
    bcftools query -f '[%GT\n]' | \
    grep -E '^(0[\/|]1|1[\/|]0)$' | wc -l)

  # 4) number of singletons carried by SAMPLE on autosomes
  singletons=$(bcftools view -r ${REGIONS} -i 'INFO/AC=1' ${VCF} | \
    bcftools query -s "${SAMPLE}" -f '[%GT\n]' | \
    grep -E '^(0[\/|]1|1[\/|]0|1[\/|]1)$' | wc -l)

  # 5) compute rates
  call_rate="NA"
  het_rate="NA"
  if [ "${total_sites}" -gt 0 ]; then
    call_rate=$(awk -v c=${called} -v t=${total_sites} 'BEGIN{printf "%.6f", (t>0 ? c/t : 0)}')
  fi
  if [ "${called}" -gt 0 ]; then
    het_rate=$(awk -v h=${hets} -v c=${called} 'BEGIN{printf "%.6f", h/c}')
  fi

  # 6) output
  echo -e "${SAMPLE}\t${VCF}\t${total_sites}\t${called}\t${call_rate}\t${hets}\t${het_rate}\t${singletons}"

done < ${SAMPLES}