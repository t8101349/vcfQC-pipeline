# check input_vcf, REF, REF.fai existed
# produce REF.fai :
# bcftools faidx /CMUH_server/DataShare/REFERENCE/SEQ/UCSC/hg38/hg38.fa

#!/bin/bash

Help() {
    echo "
使用方法：
  bash vcfqc_shapeit5.sh --input input.vcf.gz --output out_name --ref ref.fa --sex_file xxx.psam

必要參數：
  -i / --input        輸入 VCF
  -o / --output       輸出命名

可選參數：
  -r / --ref          參考基因組 (預設 hg38)
  -s / --sex_file     PSAM sex file (目前未用)
"
    exit 0
}

# === 參數解析 ===
re='^(--help|-h)$'
if [[ $1 =~ $re ]]; then
    Help
else
    while [ "$#" -gt 0 ]; do
        case "$1" in
            -r|--ref) ref="$2"; shift 2;;
            -i|--input) input="$2"; shift 2;;
            -o|--output) output="$2"; shift 2;;
            -s|--sex_file) sex_file="$2"; shift 2;;
            *) echo "unknown option: $1" >&2; exit 1;;
        esac
    done

    # === 檢查必要參數 ===
    if [[ -z "$input" || -z "$output" ]]; then
        echo "❌ 必要參數缺失，請確認 --input --output 是否有指定。" >&2
        exit 1
    fi
    if [[ -z "$ref" ]]; then
        ref="/CMUH_server/DataShare/REFERENCE/SEQ/UCSC/hg38/hg38.fa"
    fi
    if [[ -z "$sex_file" ]]; then
        sex_file="/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam"
    fi
fi

out_vcf1=/home/Weber/vcfQC/TWB_Drogen_vcfqc.vcf.gz

bcftools norm -Ou -m -any ${input} |\
  bcftools norm -Ou -f ${ref} |\
  bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' -O z -o ${out_vcf1}

# 建立索引
tabix -p vcf ${out_vcf1}

# chrx.vcf
bcftools view -r chrX \
    ${out_vcf1} -Oz -o chrX.vcf.gz
tabix -p vcf chrX.vcf.gz
# 取 PAR 區域 (PAR1 + PAR2)
bcftools view -r chrX:10001-2781479,chrX:155701383-156030895 chrX.vcf.gz -Oz -o chrX_PAR.vcf.gz
tabix -p vcf chrX_PAR.vcf.gz
# 取 nonPAR 區域 (整條 chrX 扣掉 PAR1, PAR2)
bcftools view -r chrX chrX.vcf.gz | \
    bcftools view -T ^<(echo -e "chrX\t10001\t2781479\nchrX\t155701383\t156030895") -Oz -o chrX_nonPAR.vcf.gz
tabix -p vcf chrX_nonPAR.vcf.gz

# out_vcf2=/home/Weber/vcfQC/TWB_Drogen_vcfqc_2.vcf.gz
# bcftools filter -e 'FMT/GQ<20 | FMT/GT="mis"' -S . -Oz -o ${out_vcf2} ${out_vcf1}

out_vcf2=/home/Weber/vcfQC/TWB_Drogen_vcfqc_2-2.vcf.gz
bcftools +fill-tags ${out_vcf1} -Oz -o withDP.vcf.gz -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'
bcftools filter -e 'FMT/GQ<20 | FMT/DP<5 | FMT/GT="mis"' -S . withDP.vcf.gz | \
  bcftools +fill-tags -Oz -o ${out_vcf2} -- -t F_MISSING,HWE,AF

out_vcf_par1=chrX_PAR.vcf.gz
out_vcf_par2=chrX_PAR_2.vcf.gz
bcftools +fill-tags ${out_vcf_par1} -Oz -o PARwithDP.vcf.gz -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'
bcftools filter -e 'FMT/GQ<20 | FMT/DP<5 | FMT/GT="mis"' -S . PARwithDP.vcf.gz | \
  bcftools +fill-tags -Oz -o ${out_vcf_par2} -- -t F_MISSING,HWE,AF
tabix -p vcf ${out_vcf_par2}

out_vcf_nopar1=chrX_nonPAR.vcf.gz
out_vcf_nopar1_2=chrX_nonPAR_1_2.vcf.gz
bcftools +fill-tags ${out_vcf_nopar1} \
    -Oz -o noPARwithDP.vcf.gz \
    -- -t 'FORMAT/DP:1=int(smpl_sum(FORMAT/AD))'
tabix -p vcf noPARwithDP.vcf.gz

# 過濾 + 填 F_MISSING, AF
bcftools filter -e 'FMT/GQ<20 | FMT/DP<5 | FMT/GT="mis"' -S . noPARwithDP.vcf.gz \
| bcftools +fill-tags -Oz -o ${out_vcf_nopar1_2} -- -t F_MISSING,AF
tabix -p vcf ${out_vcf_nopar1_2}

# 去除male
awk 'NR>1 && $2==2 {print $1}' ${sex_file} > female_samples.txt
bcftools view -S female_samples.txt ${out_vcf_nopar1_2} \
    -Oz -o noPAR_flt_female.vcf.gz
tabix -p vcf noPAR_flt_female.vcf.gz

# female-only 算 HWE
out_vcf_nopar2=chrX_nonPAR_2.vcf.gz
bcftools +fill-tags noPAR_flt_female.vcf.gz \
    -Oz -o ${out_vcf_nopar2} \
    -- -t HWE
tabix -p vcf ${out_vcf_nopar2}


# 先把 male-only VCF 抽出來
# sex_file=/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam
# awk 'NR>1 && $2 == 1 {print $1}' ${sex_file} > male_samples.txt

# bcftools view -S male_samples.txt ${out_vcf_nopar1_2} -Oz -o noPAR_flt_male.vcf.gz
# tabix -p vcf noPAR_flt_male.vcf.gz

# 合併 female (HWE) + male #合併會導致sample清單最後不一致
# out_vcf_nopar2_2=chrX_nonPAR_merged.vcf.gz
# bcftools merge -m all -Oz -o ${out_vcf_nopar2_2} ${out_vcf_nopar2} noPAR_flt_male.vcf.gz
# tabix -p vcf ${out_vcf_nopar2_2}

out_vcf_nopar2=chrX_nonPAR_2.vcf.gz
out_vcf_nopar2_2=chrX_nonPAR_merged.vcf.gz
# 抽出 HWE tag
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HWE\n' ${out_vcf_nopar2} \
  | bgzip -c > female_HWE.tsv.gz
tabix -s1 -b2 -e2 female_HWE.tsv.gz

# 把 HWE annotate 回 full nonPAR
bcftools annotate -a female_HWE.tsv.gz -h <(echo '##INFO=<ID=HWE,Number=1,Type=Float,Description="Hardy-Weinberg Equilibrium from females only">') \
  -c CHROM,POS,REF,ALT,INFO/HWE ${out_vcf_nopar1_2} -Oz -o ${out_vcf_nopar2_2}
tabix -p vcf ${out_vcf_nopar2_2}

# f_missing.py

out_vcf2=/home/Weber/vcfQC/TWB_Drogen_vcfqc_2-2.vcf.gz
out_vcf3=/home/Weber/vcfQC/TWB_Drogen_vcfqc_3.vcf.gz
bcftools view -i "INFO/F_MISSING <= 0.1 & INFO/HWE >= 1e-6" \
  ${out_vcf2} -Oz -o ${out_vcf3}

out_vcf_par3=chrX_PAR_3.vcf.gz
bcftools view -i "INFO/F_MISSING <= 0.1 & INFO/HWE >= 1e-6" \
  ${out_vcf_par2} -Oz -o ${out_vcf_par3}
tabix -p vcf ${out_vcf_par3}

out_vcf_nopar3=chrX_nonPAR_3.vcf.gz
bcftools view -i "INFO/F_MISSING <= 0.1 & INFO/HWE >= 1e-6" \
  ${out_vcf_nopar2_2} -Oz -o ${out_vcf_nopar3}
tabix -p vcf ${out_vcf_nopar3}


# autosome.vcf
bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 \
    ${out_vcf3} -Oz -o autosome.vcf.gz
tabix -p vcf autosome.vcf.gz


bash extract_snp_info.sh autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz combine_vcfqc_3_snp_info.tsv

# plot_snp_qc.py
python plot_snp_qc.py combine_vcfqc_3_snp_info.tsv

# make stats_file
# out_vcf3=/home/Weber/vcfQC/TWB_Drogen_vcfqc_3.vcf.gz
# out_file=/home/Weber/vcfQC/stats_file.stats
# bcftools stats -s - ${out_vcf3} > ${out_file}
# plot-vcfstats -p ./stats_plots ${out_file}
# bcftools concat：用在 同樣sample、不同染色體或區段 的情況。
# bcftools merge：用來合併 同樣變異、不同sample 的 VCF。
bash check_samples.sh autosome.vcf.gz chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz
bash bcf_stats.sh  #singleton graph

python plot_sample_qc.py stats_file.stats

out_vcf4=/home/Weber/vcfQC/TWB_Drogen_vcfqc_4.vcf.gz
bcftools view -i "INFO/F_MISSING <= 0.01 & INFO/AF >= 0.05" vcfqc_final_merged.sorted.vcf.gz -O z -o ${out_vcf4}

# relative
out_vcf5=TWB_Drogen_vcfqc_3_relative.txt
vcftools --gzvcf autosome.vcf.gz \
  --relatedness2 \
  --out ${out_vcf5}

# machine type
# platform_file=/CMUH_server/DATA7/SNParray/TWB2023/4.lab_info.csv
# python platform_mix.py sample_qc_metrics.xlsx /CMUH_server/DATA7/SNParray/TWB2023/4.lab_info.csv
# python plot_platform.py

# sex_file=/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam
python sex_mix.py

# singleton sample counts
python plot_singleton.py

# 計算 singleton 數量
bcftools view -i 'INFO/AC=1' vcfqc_final_merged.sorted.vcf.gz | bcftools view -H | wc -l

# 移除 singleton SNP
bcftools view -e 'INFO/AC=1' -Oz -o autosome_singletonremoved.vcf.gz autosome.vcf.gz
tabix -p vcf autosome_singletonremoved.vcf.gz
# X移除 singleton SNP
bcftools concat \
  -Oz -o chrX_merged.vcf.gz \
  chrX_PAR_3.vcf.gz chrX_nonPAR_3.vcf.gz
bcftools sort chrX_merged.vcf.gz -Oz -o chrX_merged.sorted.vcf.gz
tabix -p vcf chrX_merged.sorted.vcf.gz
bcftools view -e 'INFO/AC=1' -Oz -o chrX_merged_singletonremoved.vcf.gz chrX_merged.sorted.vcf.gz
tabix -p vcf chrX_merged_singletonremoved.vcf.gz

# 移除 singleton SNP
bcftools view -e 'INFO/AC<2' -Oz -o autosome_singletonremoved2.vcf.gz autosome.vcf.gz
tabix -p vcf autosome_singletonremoved2.vcf.gz
# X移除 singleton SNP
bcftools view -e 'INFO/AC<2' -Oz -o chrX_merged_singletonremoved2.vcf.gz chrX_merged.sorted.vcf.gz
tabix -p vcf chrX_merged_singletonremoved2.vcf.gz

# bcftools query -f '%FILTER\n' autosome_singletonremoved.vcf.gz | grep -v '^PASS$' | wc -l
# bcftools query -f '%FILTER\n' chrX_merged_singletonremoved.vcf.gz | grep -v '^PASS$' | wc -l
bcftools query -f '%FILTER\n' autosome_singletonremoved2.vcf.gz | grep -v '^PASS$' | wc -l
bcftools query -f '%FILTER\n' chrX_merged_singletonremoved2.vcf.gz | grep -v '^PASS$' | wc -l


# sex check
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


# index查看每條染色體的變異數
for chr in chr{1..22} chrX chrY chrM; do
  echo -ne "$chr\t"
  bcftools view -r $chr ${out_vcf1} | bcftools view -H | wc -l
done
for chr in chr{1..22} chrX; do
  echo -ne "$chr\t"
  bcftools view -r $chr ${out_vcf3} | bcftools view -H | wc -l
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
#paralell
for CHR in {1..22}; do
    echo "==> Phasing chr${CHR} with GNU parallel..."

    MAP="/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz"
    INPUTVCF="phasing/chr${CHR}.vcf.gz"
    CHUNKS="shapeit5/phasing/chr${CHR}.chunks.txt"
    OUTDIR="shapeit5/shapeit_chr"
    LOGDIR="shapeit5/logs"
    THREAD=5

    export CHR MAP INPUT OUTDIR LOGDIR THREAD

    parallel -j 10 --env CHR --env MAP --env INPUT --env OUTDIR --env LOGDIR --env THREAD '
        REGION={}

        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H "${INPUTVCF}" | wc -l)
        echo "  SNP count: $lines"

        if [ "$lines" -gt 0 ]; then
            SHAPEIT5_phase_common \
                --input "${INPUTVCF}" \
                --map "${MAP}" \
                --region "${REGION}" \
                --output "${OUTDIR}/chr${CHR}_${REGION}.phased.bcf" \
                --thread "${THREAD}" \
                --log "${LOGDIR}/chr${CHR}_${REGION}.log"
        else
            echo "⚠️  REGION ${REGION} has 0 SNPs, skipping..."
        fi
    ' :::: "${CHUNKS}"

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

for CHR in {1..22}; do
    bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
    bcftools query -l shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf | wc -l
done


# chrx
mkdir -p shapeit5/shapeit_chrx
mkdir -p shapeit5/shapeit_chrx/chrX
mkdir -p shapeit5/shapeit_chrx/chrX_par1
mkdir -p shapeit5/shapeit_chrx/chrX_par2
mkdir -p shapeit5/shapeit_chrx/complete


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
         else if(gt=="0|1" || gt=="1|0") a[1]="."
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
      --output shapeit5/shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_nonpar_male_haploids.phased.bcf \
      --thread 10 \
      --log shapeit5/logs/chrX_nonpar_male_haploids.log
done


# merge nonpar
bcftools merge -Oz -o shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz \
  shapeit5/shapeit_chrx/complete/chrX_1.bcf \
  shapeit5/shapeit_chrx/chrX/chrX_nonpar_male_haploids.phased.bcf
bcftools index -c shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz

bcftools query -l shapeit5/shapeit_chrx/complete/chrX_par1_1.bcf > chrX_order.txt
bcftools view -S chrX_order.txt -Oz -o shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz shapeit5/shapeit_chrx/complete/chrX_3.vcf.gz
bcftools index -c shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz


# concat
bcftools concat -Ob -o shapeit5/shapeit_chr/complete/merged_chrX.target.vcf.gz \
  shapeit5/shapeit_chrx/complete/chrX_par1_1.bcf \
  shapeit5/shapeit_chrx/complete/chrX_3.sorted.vcf.gz \
  shapeit5/shapeit_chrx/complete/chrX_par2_1.bcf
bcftools index -c shapeit5/shapeit_chr/complete/merged_chrX.target.vcf.gz


for CHR in X; do
    bcftools view -H shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf | wc -l
    bcftools query -l shapeit5/shapeit_chr/complete/merged_chr${CHR}.target.bcf | wc -l
done

#CHECK
bcftools view -H cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l
bcftools query -l cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l

for CHR in {1..22} X; do
    if [[ "$CHR" == "X" ]]; then
        INFILE="shapeit5/shapeit_chr/complete/merged_chrX.target.bcf"
    else
        INFILE="shapeit5/shapeit_chr/complete/merged_chr${CHR}.bcf"
    fi

    OUTFILE="shapeit5/shapeit_chr/complete/${output}_merged_chr${CHR}.target.vcf.gz"

    echo "Processing chr${CHR} ..."
    bcftools view -Oz -o "${OUTFILE}" "${INFILE}"
    bcftools index -t "${OUTFILE}"
done

###
# shapeit5 --> phase_common(MAF>0.1%) > phase_common(related,haploids,MAF>0.1%) > phase_rare(scaffold)
###
