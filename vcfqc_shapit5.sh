# sex check
sex_file=/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam
bfile_name=chrX_singletonremoved
plink2 --vcf chrX_merged_singletonremoved.vcf.gz \
  --update-sex ${sex_file} \
  --split-par hg38 \
  --make-pgen \
  --threads 64 \
  --out ${bfile_name}
plink2 --vcf chrX_merged_singletonremoved.vcf.gz \
  --update-sex ${sex_file} \
  --split-par hg38 \
  --make-bed \
  --threads 64 \
  --out ${bfile_name}

awk '{if($1=="PAR1" || $1=="PAR2") $1="23"; print $0}' chrX_singletonremoved.bim > chrX_singletonremoved_fixed.bim
mv chrX_singletonremoved_fixed.bim chrX_singletonremoved.bim
awk '{print "0",$1}' female_samples.txt > female_samples_keep.txt

bfile=/home/Weber/vcfQC/chrX_merged_singletonremoved.vcf.gz
bash sex_check.sh ${bfile} ${sex_file} sex_check


# bash sexcheck_fail_compute.sh autosome_singletonremoved.vcf.gz sex_check/chrX_sexcheck_failed_samples.txt
# 直接找sample_qc_metrics_sex.csv就好

# bcftools view -S <(bcftools query -l autosome_singletonremoved.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
#   -Oz -o cleaned_autosome_singletonremoved.vcf.gz autosome_singletonremoved.vcf.gz
# tabix -p vcf cleaned_autosome_singletonremoved.vcf.gz

# bcftools view -S <(bcftools query -l chrX_merged_singletonremoved.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
#   -Oz -o cleaned_chrX_merged_singletonremoved.vcf.gz chrX_merged_singletonremoved.vcf.gz
# tabix -p vcf cleaned_chrX_merged_singletonremoved.vcf.gz

bcftools view -S <(bcftools query -l autosome_singletonremoved2.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
  -Oz -o cleaned_autosome_singletonremoved2.vcf.gz autosome_singletonremoved2.vcf.gz
tabix -p vcf cleaned_autosome_singletonremoved2.vcf.gz

bcftools view -S <(bcftools query -l chrX_merged_singletonremoved2.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
  -Oz -o cleaned_chrX_merged_singletonremoved2.vcf.gz chrX_merged_singletonremoved2.vcf.gz
tabix -p vcf cleaned_chrX_merged_singletonremoved2.vcf.gz

# PCA


# index查看每條染色體的變異數
for chr in chr{1..22} chrX chrY chrM; do
  echo -ne "$chr\t"
  bcftools view -r $chr TWB_Drogen_vcfqc.vcf.gz | bcftools view -H | wc -l
done
for chr in chr{1..22} chrX; do
  echo -ne "$chr\t"
  bcftools view -r $chr TWB_Drogen_vcfqc_3.vcf.gz | bcftools view -H | wc -l
done

# autosome
before=$(bcftools view -H autosome.vcf.gz | wc -l)
after=$(bcftools view -H autosome_singletonremoved.vcf.gz | wc -l)
after2=$(bcftools view -H autosome_singletonremoved2.vcf.gz | wc -l)
echo -e "autosome\t$before\t$after\t$after2" >> variant_counts.tsv

# chrX
before=$(bcftools view -H chrX_merged.sorted.vcf.gz | wc -l)
after=$(bcftools view -H chrX_merged_singletonremoved.vcf.gz | wc -l)
after2=$(bcftools view -H chrX_merged_singletonremoved2.vcf.gz | wc -l)
echo -e "chrX\t$before\t$after\t$after2" >> variant_counts.tsv


# bcftools annotate --remove FORMAT/LAF cleaned_autosome_singletonremoved.vcf.gz -Oz -o clean_autosome_singletonremoved_NOLAF.vcf.gz
# tabix -p vcf clean_autosome_singletonremoved_NOLAF.vcf.gz

# bcftools annotate --remove FORMAT/LAF cleaned_chrX_merged_singletonremoved.vcf.gz -Oz -o cleaned_chrX_merged_singletonremoved_NOLAF.vcf.gz
# tabix -p vcf cleaned_chrX_merged_singletonremoved_NOLAF.vcf.gz


# phasing開始
# 分割autosome
mkdir -p phasing
# INPUT=clean_autosome_singletonremoved_NOLAF.vcf.gz
INPUT=cleaned_autosome_singletonremoved2.vcf.gz

for CHR in {1..22} X; do
    echo "Extracting chr${CHR}..."
    bcftools view -r chr${CHR} ${INPUT} -Oz -o phasing/chr${CHR}.vcf.gz
    tabix -p vcf phasing/chr${CHR}.vcf.gz
done

# region切分
for CHR in {1..22}; do
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    zcat ${MAP} | awk -v chr=${CHR} 'NR > 1 {
            # 初始化
            if (NR == 2) {
                start_bp = $1
                start_cM = $3
                next
            }
            
            # 每累積超過 20 cM 就輸出
            if ($3 - start_cM >= 20) {
                print "chr" chr ":" start_bp "-" last_bp
                start_bp = $1 + 1
                start_cM = $3
            }
            
            last_bp = $1
        }
        END {
            # 若最後一段未輸出則補上
            if (last_bp > start_bp)
                print "chr" chr ":" start_bp "-" last_bp
        }' > shapeit5/phasing/chr${CHR}.chunks.txt
done
# chr15需手動調整region大小(稀疏密集差太多)


# 可能需要切開男女sample 記憶體爆炸
# grep -v '^NGS20140610G[[:space:]]*' female_samples.txt > female_samples.new.txt

# # female sample subset
# for CHR in {1..22}; do
#     echo "Extracting female chr${CHR}..."
#     bcftools view -S female_samples.new.txt phasing/chr${CHR}.vcf.gz -Oz -o phasing/female.chr${CHR}.vcf.gz
#     tabix -p vcf phasing/female.chr${CHR}.vcf.gz
# done

# # male sample subset
# for CHR in {1..22}; do
#     echo "Extracting male chr${CHR}..."
#     bcftools view -S ^female_samples.new.txt phasing/chr${CHR}.vcf.gz -Oz -o phasing/male.chr${CHR}.vcf.gz
#     tabix -p vcf phasing/male.chr${CHR}.vcf.gz
# done


# reference
cut -f1 /CMUH_server/home/Rogen/Resource/1000Genome/population/1000Genome_Asia_Samples.tsv | tail -n +2 > EAS_samples.list
bcftools query -l /CMUH_server/home/Rogen/Resource/1000Genome/dataset/vcf_add_id/add_id/ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz \
  | grep -F -f EAS_samples.list > EAS_samples_inVCF.list

for CHR in {1..22} X; do
    bcftools view --threads 4 -S EAS_samples_inVCF.list -Oz -o 1000G_phased/1k_genome_EAS_phased_chr${CHR}.vcf.gz /CMUH_server/home/Rogen/Resource/1000Genome/dataset/vcf_add_id/add_id/ALL.chr${CHR}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz 
    bcftools index -c --threads 4 1000G_phased/1k_genome_EAS_phased_chr${CHR}.vcf.gz
done

# bcftools norm -m-any -Oz -o target.norm.vcf.gz target.vcf.gz
# bcftools norm -m-any -Oz -o ref.norm.vcf.gz reference.vcf.gz

mkdir -p shapeit5/logs
mkdir -p shapeit5/shapeit_chr
for CHR in {1..6}; do
    echo "Phasing chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H phasing/chr${CHR}.vcf.gz | wc -l)
        echo " SNP count: $lines "
        SHAPEIT5_phase_common \
          --input phasing/chr${CHR}.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output /home/Weber/vcfQC/shapeit5/shapeit_chr/chr${CHR}_${REGION}.phased.bcf \
          --thread 5 \
          --log shapeit5/logs/chr${CHR}_${REGION}.log
    done < <(cat shapeit5/phasing/chr${CHR}.chunks.txt)
done

for CHR in {7..13}; do
    echo "Phasing chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H phasing/chr${CHR}.vcf.gz | wc -l)
        echo " SNP count: $lines "
        SHAPEIT5_phase_common \
          --input phasing/chr${CHR}.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output /home/Weber/vcfQC/shapeit5/shapeit_chr/chr${CHR}_${REGION}.phased.bcf \
          --thread 5 \
          --log shapeit5/logs/chr${CHR}_${REGION}.log
    done < <(cat shapeit5/phasing/chr${CHR}.chunks.txt)
done

for CHR in {13..22}; do
    echo "Phasing chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H phasing/chr${CHR}.vcf.gz | wc -l)
        echo " SNP count: $lines "
        SHAPEIT5_phase_common \
          --input phasing/chr${CHR}.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output /home/Weber/vcfQC/shapeit5/shapeit_chr/chr${CHR}_${REGION}.phased.bcf \
          --thread 5 \
          --log shapeit5/logs/chr${CHR}_${REGION}.log
    done < <(cat shapeit5/phasing/chr${CHR}.chunks.txt)
done


mkdir -p shapeit5/shapeit_chr/complete


for CHR in {1..22}; do
    echo "Merging chr${CHR}..."

    # 生成 chunk 檔案清單，依 phasing/chr${CHR}.chunks.txt 排序
    chunks=$(awk -v chr=${CHR} '{print "shapeit5/shapeit_chr/chr" chr "_" $1 ".phased.bcf"}' shapeit5/phasing/chr${CHR}.chunks.txt)

    # concat 並直接排序，輸出到 complete
    bcftools concat -Ob $chunks | bcftools sort -Ob -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf

    # 建立索引
    bcftools index -c shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf
done

'''
# sample 確認
for CHR in {1..22}; do
    bcftools view \
      -S target_samples.txt \
      -Ob -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf \
      shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf
    
    bcftools index -c shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf

    echo -n "chr${CHR}: "
    bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf | wc -l
done
'''
for CHR in {1..22}; do
    bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
    bcftools query -l shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
done

for CHR in {1..22}; do
    bcftools view -Oz -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz \
        shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf
    bcftools index -t shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz
done

# chrx
mkdir -p shapeit5/shapeit_chrx
mkdir -p shapeit5/shapeit_chrx/chrX
mkdir -p shapeit5/shapeit_chrx/chrX_par1
mkdir -p shapeit5/shapeit_chrx/chrX_par2
mkdir -p shapeit5/shapeit_chrx/complete
CHR=X

# nonpar divide male and female
# female phasing
echo "divide female chrx..."
bcftools view -S female_samples.txt chrX_merged_singletonremoved2.vcf.gz \
    -Oz -o chrX_merged_singletonremoved2_female.vcf.gz
tabix -p vcf chrX_merged_singletonremoved2_female.vcf.gz

bcftools view -S <(bcftools query -l chrX_merged_singletonremoved2_female.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
  -Oz -o cleaned_chrX_merged_singletonremoved2_female.vcf.gz chrX_merged_singletonremoved2_female.vcf.gz
tabix -p vcf cleaned_chrX_merged_singletonremoved2_female.vcf.gz

bcftools view -S <(bcftools query -l chrX_merged_singletonremoved2.vcf.gz | grep -v -f sex_check/chrX_sexcheck_failed_samples.txt) \
  -Oz -o cleaned_chrX_merged_singletonremoved2.vcf.gz chrX_merged_singletonremoved2.vcf.gz
tabix -p vcf cleaned_chrX_merged_singletonremoved2.vcf.gz

echo "Before:" $(bcftools query -l chrX_merged_singletonremoved2.vcf.gz | wc -l)
echo "Before:" $(bcftools query -l chrX_merged_singletonremoved2_female.vcf.gz | wc -l)
echo "After:" $(bcftools query -l cleaned_chrX_merged_singletonremoved2_female.vcf.gz | wc -l)

MAP_DIR=/home/Weber/hg38_chr
for REGION in chrX chrX_par1 chrX_par2; do
    MAP=${MAP_DIR}/chr${REGION}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Processing ${REGION}..."
    
    zcat ${MAP} | awk -v region=${REGION} '
        NR > 1 {
            # 初始化
            if (NR == 2) {
                start_bp = $1
                start_cM = $3
                next
            }
            
            # 每累積超過 20 cM 就輸出
            if ($3 - start_cM >= 20) {
                end_bp = $1
                print "chrX" ":" start_bp "-" end_bp
                start_bp = $1 + 1
                start_cM = $3
            }
            
            last_bp = $1
        }
        END {
            # 若最後一段未輸出則補上
            if (last_bp > start_bp)
                print "chrX" ":" start_bp "-" last_bp
        }
    ' > shapeit5/phasing/${REGION}.chunks.txt
done


MAP_DIR=/home/Weber/hg38_chr
# par
for REGION_TYPE in chrX_par1 chrX_par2; do
    MAP=/home/Weber/hg38_chr/chr${REGION_TYPE}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Phasing ${REGION_TYPE}..."
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l)
        echo " SNP count: $lines "
        SHAPEIT5_phase_common \
          --input cleaned_chrX_merged_singletonremoved2.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output /home/Weber/vcfQC/shapeit5/shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_${REGION}.phased.bcf \
          --thread 10 \
          --log shapeit5/logs/chrX_${REGION}.log
    done < shapeit5/phasing/${REGION_TYPE}.chunks.txt
done

# cleaned_chrX_merged_singletonremoved2 的 par2 沒有snp
#nonpar
for REGION_TYPE in chrX; do
    MAP=/home/Weber/hg38_chr/chr${REGION_TYPE}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Phasing ${REGION_TYPE}..."
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H cleaned_chrX_merged_singletonremoved2_female.vcf.gz | wc -l)
        echo " SNP count: $lines "
        SHAPEIT5_phase_common \
          --input cleaned_chrX_merged_singletonremoved2_female.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output shapeit5/shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_${REGION}.phased.bcf \
          --thread 10 \
          --log shapeit5/logs/chrX_${REGION}.log
    done < shapeit5/phasing/${REGION_TYPE}.chunks.txt
done


for REGION_TYPE in chrX chrX_par1 chrX_par2; do
    echo "Concatenating in ${REGION_TYPE}..."
    bcftools concat -Ob -o shapeit5/shapeit_chrx/complete/${REGION_TYPE}_1.bcf \
  $(awk -v grp="$REGION_TYPE" '{print "shapeit5/shapeit_chrx/" grp "/" grp "_" $1 ".phased.bcf"}' shapeit5/phasing/${REGION_TYPE}.chunks.txt)
    bcftools index -c shapeit5/shapeit_chrx/complete/${REGION_TYPE}_1.bcf
done


# nonpar_male
bcftools view -S male_samples.txt cleaned_chrX_merged_singletonremoved2.vcf.gz \
  | sed 's/\/|/|/g' | sed 's/\//|/g' \
  | bgzip -c > cleaned_chrX_merged_singletonremoved2_male.vcf.gz
tabix -p vcf cleaned_chrX_merged_singletonremoved2_male.vcf.gz

bcftools view -r chrX:2785350-155699751 cleaned_chrX_merged_singletonremoved2_male.vcf.gz \
    | bgzip -c > shapeit5/shapeit_chrx/complete/chrX_male_nonpar.vcf.gz
tabix -p vcf shapeit5/shapeit_chrx/complete/chrX_male_nonpar.vcf.gz

# 修正
{ bcftools view -h shapeit5/shapeit_chrx/complete/chrX_male_nonpar.vcf.gz;
  bcftools view -H shapeit5/shapeit_chrx/complete/chrX_male_nonpar.vcf.gz \
  | awk 'BEGIN{OFS="\t"}{
       # 前 9 欄不變
       for(i=1;i<=9;i++) printf "%s\t",$i;

       # 從第10欄開始是 sample
       for(i=10;i<=NF;i++){
         # 拆出 GT 與其他欄位
         split($i, a, ":")
         gt=a[1]

         # 修改 GT
         if(gt=="0|0" || gt=="0/0") a[1]="0"
         else if(gt=="1|1" || gt=="1/1") a[1]="1"
         else if(gt=="0|1" || gt=="1|0") a[1]="1"
         else if(gt==".|." || gt=="./.") a[1]="."

         # 重組 sample
         printf "%s", a[1]
         for(j=2;j<=length(a);j++) printf ":%s", a[j]
         printf "%s", (i==NF?"\n":"\t")
       }
  }'
} | bgzip -c > shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed.vcf.gz
tabix -p vcf shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed.vcf.gz

# fixploidy
bcftools +fixploidy shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed.vcf.gz -Ob -o shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed_diploid.vcf.gz
bcftools index -t shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed_diploid.vcf.gz

# nonpar_male --haploids
for REGION_TYPE in chrX; do
    MAP=/home/Weber/hg38_chr/chr${REGION_TYPE}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Phasing ${REGION_TYPE}..."
    SHAPEIT5_phase_common \
      --input shapeit5/shapeit_chrx/complete/chrX_male_nonpar_fixed_diploid.vcf.gz \
      --haploids male_samples.txt\
      --map ${MAP} \
      --region chrX \
      --output shapeit5/shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_nonpar_male_haploids.phased.vcf.gz \
      --thread 10 \
      --log shapeit5/logs/chrX_nonpar_male_haploids.log
done
bcftools index -t shapeit5/shapeit_chrx/complete/chrX_nonpar_male_haploids.phased.vcf.gz


# merge nonpar
bcftools merge -Oz -o shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz \
  shapeit5/shapeit_chrx/complete/chrX_1.bcf \
  shapeit5/shapeit_chrx/complete/chrX_nonpar_male_haploids.phased.vcf.gz
bcftools index -c shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz

bcftools query -l shapeit5/shapeit_chrx/complete/chrX_par1_1.bcf > chrX_order.txt
bcftools view -S chrX_order.txt -Oz -o shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz
bcftools index -c shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz


# concat
bcftools concat -Ob -o shapeit5/shapeit_chr/complete/merged_chrX.target.bcf \
  shapeit5/shapeit_chrx/complete/chrX_par1_1.bcf \
  shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz \
  shapeit5/shapeit_chrx/complete/chrX_par2_1.bcf
bcftools index -c shapeit5/shapeit_chr/complete/merged_chrX.target.bcf


for CHR in X; do
    bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf | wc -l
    bcftools query -l shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf | wc -l
done
#CHECK
bcftools view -H cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l
bcftools query -l cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l

for CHR in X; do
    bcftools view -Oz -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz \
        shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf
    bcftools index -t shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz
done





"""
# debug
CHR=15
echo "Phasing chr${CHR}..."
MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
for groupfile in group_phasing/${CHR}_samples_group_*; do
    [[ $groupfile == *.csi ]] && continue

    group=$(basename "$groupfile" .vcf.gz)
    echo "GROUP $group..."
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H group_phasing/${group}.vcf.gz | wc -l)
        echo " SNP count: $lines "
        outfile="shapeit_chr4/chr${group}_${REGION}.phased.vcf.gz"
        if [ ! -f "$outfile" ]; then
            echo "  → Running shapeit4 for ${REGION}"
            shapeit4 \
              --input group_phasing/${group}.vcf.gz \
              --map ${MAP} \
              --region ${REGION} \
              --output "$outfile" \
              --thread 5 \
              --pbwt-depth 1 \
              --window 3.5 \
              --sequencing \
              --log shapeit_chr4/logs/${group}_${REGION}.log
        else
            echo "  → Skipping ${REGION} (already done)"
        fi
    done < <(cat phasing/chr${CHR}.chunks.txt)
done

echo "Concatenating groups in chr${CHR}..."
for groupfile in group_phasing/${CHR}_samples_group_*; do
    [[ $groupfile == *.csi ]] && continue
    group=$(basename "$groupfile" .vcf.gz)
    echo "concat ${group}"
    bcftools concat -Oz -o shapeit_chr4/merged_${group}.vcf.gz \
$(awk -v grp="$group" '{print "shapeit_chr4/chr" grp "_" $1 ".phased.vcf.gz"}' phasing/chr${CHR}.chunks.txt)
    bcftools index -c shapeit_chr4/merged_${group}.vcf.gz
done


echo "Merging chr${CHR}..."
bcftools merge --force-samples -Oz -o shapeit_chr4/merged_chr${CHR}.vcf.gz \
    shapeit_chr4/merged_${CHR}_samples_group_*.vcf.gz
bcftools index -c shapeit_chr4/merged_chr${CHR}.vcf.gz

bcftools query -l shapeit_chr4/merged_chr${CHR}.target.vcf.gz | wc -l


for CHR in {1..22}; do
    bcftools view \
      -S target_samples.txt \
      -Oz -o temp/merged_chr${CHR}.target.vcf.gz \
      temp/merged_chr${CHR}.vcf.gz
    
    bcftools index -c temp/merged_chr${CHR}.target.vcf.gz

    echo -n "chr${CHR}: "
    bcftools view -H temp/merged_chr${CHR}.target.vcf.gz | wc -l
done


# debug'

MAP=/home/Weber/hg38_chr/chr7.chrprefix.b38.gmap.gz
shapeit4 \
  --input group_phasing/7_samples_group_ak.vcf.gz \
  --map ${MAP} \
  --region chr7:77477480-80428848 \
  --output shapeit_chr4/chr7_samples_group_ab_chr7:77477480-80428848.phased.vcf.gz \
  --thread 2 \
  --pbwt-depth 1 \
  --window 2.0 \
  --sequencing \
  --log shapeit_chr4/logs/chr7_samples_group_ab_chr7:77477480-80428848.log

bcftools view -H shapeit_chr4/chr7_samples_group_ab_chr7:77477480-80428848_debug.phased.vcf.gz | wc -l
bcftools view -H shapeit_chr4/chr7_samples_group_ab_chr7:77477480-80428848.phased.vcf.gz | wc -l
bcftools view -H shapeit_chr4/chr7_samples_group_ak_chr7:77477480-80428848_debug.phased.vcf.gz | wc -l
bcftools view -H shapeit_chr4/chr7_samples_group_ak_chr7:77477480-80428848.phased.vcf.gz | wc -l
MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
chr7:77477480-80428848
shapeit_chr4/group/merged_7_samples_group_ak.vcf.gz
shapeit_chr4/chr7_samples_group_ak_chr7:77477480-80428848.phased.vcf.gz

for CHR in 7; do
    echo "Phasing chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    for group in 7_samples_group_ab; do
        echo "GROUP $group..."
        while read REGION; do
            echo "  REGION ${REGION}"
            lines=$(bcftools view -r ${REGION} -H group_phasing/${group}.vcf.gz | wc -l)
            echo " SNP count: $lines "
            shapeit4 \
              --input group_phasing/${group}.vcf.gz \
              --map ${MAP} \
              --region ${REGION} \
              --output shapeit_chr4/chr${group}_${REGION}.phased.vcf.gz \
              --thread 2 \
              --pbwt-depth 1 \
              --window 2.0 \
              --sequencing \
              --log shapeit_chr4/logs/${group}_${REGION}.log
        done < <(cat phasing/chr${CHR}.chunks.txt)
    done
done
# Step1: concat 每個 group 的 chunks
for CHR in 7; do
    echo "Concatenating groups in chr${CHR}..."
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue
        group=$(basename "$groupfile" .vcf.gz)
        echo "concat ${group}"
        bcftools concat -Oz -o shapeit_chr4/group/merged_${group}.vcf.gz \
  $(awk -v grp="$group" '{print "shapeit_chr4/chr" grp "_" $1 ".phased.vcf.gz"}' phasing/chr${CHR}.chunks.txt)
        bcftools index -c shapeit_chr4/group/merged_${group}.vcf.gz
    done
done


# Step2: merge 每條 chr 的所有 group
for CHR in 7; do
    echo "Merging chr${CHR}..."
    bcftools merge --force-samples -Oz -o shapeit_chr4/complete/merged_chr${CHR}.vcf.gz \
        shapeit_chr4/group/merged_${CHR}_samples_group_*.vcf.gz
    bcftools index -c shapeit_chr4/complete/merged_chr${CHR}.vcf.gz
done

# sample 確認
for CHR in 7; do
    bcftools view \
      -S target_samples.txt \
      -Oz -o shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz \
      shapeit_chr4/complete/merged_chr${CHR}.vcf.gz
    
    bcftools index -c shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz

    echo -n "chr${CHR}: "
    bcftools view -H shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz | wc -l
done


###parallel

#!/bin/bash

THREAD=10

for CHR in {1..2}; do
    echo "==> Phasing chr${CHR} with GNU parallel..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    INPUT=phasing/chr${CHR}.vcf.gz
    CHUNKS=shapeit5/phasing/chr${CHR}.chunks.txt
    OUTDIR=shapeit5/shapeit_chr
    LOGDIR=shapeit5/logs

    mkdir -p ${OUTDIR}/test ${LOGDIR}

    export CHR MAP INPUT OUTDIR LOGDIR THREAD

    # 使用 parallel 執行每個 REGION
    cat ${CHUNKS} | grep -v '^$' | parallel -j ${THREAD} "
        REGION={}
        echo '  REGION ${REGION}'
        lines=\$(bcftools view -r \${REGION} -H ${INPUT} | wc -l)
        echo '  SNP count: '\$lines

        if [ \$lines -gt 0 ]; then
            SHAPEIT5_phase_common \
              --input ${INPUT} \
              --map ${MAP} \
              --region \${REGION} \
              --output ${OUTDIR}/test/chr${CHR}_\${REGION}.phased.bcf \
              --thread 5 \
              --log ${LOGDIR}/chr${CHR}_\${REGION}.log
        else
            echo '⚠️  REGION \${REGION} has 0 SNPs, skipping...'
        fi
    "
done


mkdir -p shapeit5/shapeit_chr/complete

# 1️⃣ 平行 concat 每條 chr 的 chunk
# 將每個 chr 的工作封裝成函數
merge_chr() {
    CHR=$1
    echo "Merging chr${CHR}..."

    # 生成 chunk 檔案清單，依 phasing/chr${CHR}.chunks.txt 排序
    chunks=$(awk -v chr=${CHR} '{print "shapeit5/shapeit_chr/chr" chr "_" $1 ".phased.bcf"}' shapeit5/phasing/chr${CHR}.chunks.txt)

    # concat 並排序，輸出到 complete
    bcftools concat -Ob $chunks | bcftools sort -Ob -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf

    # 建立索引
    bcftools index -c shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf
}

export -f merge_chr

# 用 parallel 處理多個 chr，這裡同時跑 5 個 chr
parallel -j 5 merge_chr ::: {1..22}
'

# 2️⃣ 平行篩 target samples 並 index
seq 1 22 | parallel -j 5 '
CHR={}
echo -n "chr${CHR}: "
bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
bcftools query -l shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
bcftools view -Oz -o shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz \
    shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf
bcftools index -t shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.vcf.gz
'
"""