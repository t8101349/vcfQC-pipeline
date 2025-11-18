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

# module load R/3.4.0
# source ~/R-3.4.0-ownlib.bash
# 03_sex.check.sh ${bfile_name} ${sex_file} sex_check

bfile=$1
sexinfo=$2
outdir=$3
scrdir='/CMUH_server/DATA6/Rogen/Impute_pipeline/script'
female_keep='/home/Weber/vcfQC/female_samples_keep.txt'

mkdir -p $outdir
# cd $outdir

# Perform LD pruning: chrX
name=$(basename $bfile)


# Perform LD pruning: chrX
echo "[INFO] LD pruning on chrX..."

plink \
  --bfile $bfile \
  --keep ${female_keep} \
  --chr 23 \
  --geno 0.02 \
  --maf 0.05 \
  --snps-only just-acgt \
  --indep-pairwise 200 100 0.1 \
  --out $outdir/${name}-ldpr-chrx

if [ ! -s $outdir/${name}-ldpr-chrx.prune.in ]; then
  echo "[WARN] No SNPs passed LD pruning on chrX."
  exit 1
fi

# Check sex
echo "[INFO] Running sex check..."
plink \
  --bfile $bfile \
  --update-sex $sexinfo 1 2 \
  --extract $outdir/${name}-ldpr-chrx.prune.in \
  --check-sex 0.25 0.75 \
  --out $outdir/${name}_SNPSex_0.25_0.75



# Rscript $scrdir/03_sex.dist.R $outdir $name 0.2 0.8
Rscript $scrdir/03_sex.dist.R $outdir $name 0.25 0.75

