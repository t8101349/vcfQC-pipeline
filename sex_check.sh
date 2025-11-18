# =========================================
# Describe: Step3-Sex QC
# Parameter: $1 > Raw bfile pathway
#            $2 > Sex info. file pathway
#            $3 > Output pathway
# Author: Rogen
# Date: 2023.03.28
# Last Update: 2024.04.18
# Update Info:
# 2024.04.18 - add para. for input sex info.
# =========================================
#!/bin/bash
set -euo pipefail

# module load R/3.4.0
# source ~/R-3.4.0-ownlib.bash
# 03_sex.check.sh ${bfile_name} ${sex_file} sex_check

vcf=$1                 # chrX_merged_singletonremoved.vcf.gz
sexfile=$2             # TWBK_1484_allvariants.psam
outdir=$3              # sex_check
scrdir='/CMUH_server/DATA6/Rogen/Impute_pipeline/script'
female_keep='/home/Weber/vcfQC/female_samples_keep.txt'


mkdir -p $outdir

# 將 sex info 轉成 PLINK2 格式 (IID SEX)
awk 'NR>1{print $1, $2}' $sexfile > $outdir/sexinfo_plink2.txt

# 1️⃣ 轉 PLINK2 PGEN 格式，並自動拆 PAR
plink2 --vcf $vcf \
    --split-par hg38 \
    --update-sex $outdir/sexinfo_plink2.txt \
    --make-pgen \
    --out $outdir/chrX_plink2

# 2️⃣ LD pruning on chrX (non-PAR + PAR separately)
plink2 --pfile $outdir/chrX_plink2 \
    --chr 23 \
    --geno 0.02 \
    --maf 0.05 \
    --indep-pairwise 200 100 0.1 \
    --out $outdir/chrX_ldprune

# 3️⃣ Sex check
plink2 \
  --pfile sex_check/chrX_plink2 \
  --extract sex_check/chrX_ldprune.prune.in \
  --check-sex max-female-xf=0.25 min-male-xf=0.75 \
  --out sex_check/chrX_sexcheck


# 4️⃣ 可選：輸出 R script 畫圖
Rscript 03_sex.dist.R sex_check chrX_sexcheck 0.25 0.75


