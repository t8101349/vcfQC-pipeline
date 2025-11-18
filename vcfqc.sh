# check input_vcf, REF, REF.fai existed
# produce REF.fai :
# bcftools faidx /CMUH_server/DataShare/REFERENCE/SEQ/UCSC/hg38/hg38.fa

input_vcf=/CMUH_server/DATA11/Research/dragen_run_joint/VCF/TWBK_1484Cases/TWBK_1484Cases.vcf.gz
out_vcf1=/home/Weber/vcfQC/TWB_Drogen_vcfqc.vcf.gz
REF=/CMUH_server/DataShare/REFERENCE/SEQ/UCSC/hg38/hg38.fa

bcftools norm -Ou -m -any ${input_vcf} |\
  bcftools norm -Ou -f ${REF} |\
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
sex_file=/CMUH_server/home2/liuTY/stroke_plan_2021/TWBK_1484_allvariants.psam
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
output=TWB_Drogen_vcfqc_3_relative.txt
vcftools --gzvcf autosome.vcf.gz \
  --relatedness2 \
  --out ${output}

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
            
            # 每累積超過 4 cM 就輸出
            if ($3 - start_cM >= 4) {
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
        }' > phasing/chr${CHR}.chunks.txt
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

# 依舊爆炸 100人一組(最小100)
bcftools query -l cleaned_chrX_merged_singletonremoved2.vcf.gz > sampleslist.txt
for CHR in {1..22} X; do
    inputvcf=phasing/chr${CHR}.vcf.gz
    split -l 100 sampleslist.txt group_phasing/samples_group_
    for f in group_phasing/samples_group_*; do
        group=$(basename $f)
        echo "Processing $group..."
        bcftools view -S $f -Oz -o group_phasing/${CHR}_${group}.vcf.gz ${inputvcf}
        bcftools index group_phasing/${CHR}_${group}.vcf.gz
    done
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

'''
# 100人chunk, (0.7-1.4)1.2cM
for CHR in {1..22}; do
    echo "Phasing female chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue

        group=$(basename "$groupfile" .vcf.gz)
        echo "GROUP $group..."
        while read REGION; do
            echo "  REGION ${REGION}"
            lines=$(bcftools view -r ${REGION} -H group_phasing/${group}.vcf.gz | wc -l)
            echo " SNP count: $lines "
            shapeit4 \
              --input group_phasing/${group}.vcf.gz \
              --map ${MAP} \
              --region ${REGION} \
              --output shapeit_chr/chr${group}_${REGION}.phased.vcf.gz \
              --thread 1 \
              --pbwt-depth 1 \
              --window 0.6 \
              --sequencing
        done < <(cat phasing/chr${CHR}.chunks.txt)
    done
done

# 100人chunk, (0.7-1.4)1.2cM reference
for CHR in {1..22}; do
    echo "Phasing female chr${CHR}..."
    MAP=/home/Weber/hg38_chr/chr${CHR}.chrprefix.b38.gmap.gz
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue

        group=$(basename "$groupfile" .vcf.gz)
        echo "GROUP $group..."
        while read REGION; do
            echo "  REGION ${REGION}"
            lines=$(bcftools view -r ${REGION} -H group_phasing/${group}.vcf.gz | wc -l)
            echo " SNP count: $lines "
            shapeit4 \
              --input group_phasing/${group}.vcf.gz \
              --map ${MAP} \
              --region ${REGION} \
              --output shapeit_chr2/chr${group}_${REGION}.phased.vcf.gz \
              --thread 1 \
              --pbwt-depth 1 \
              --window 0.6 \
              --reference /home/Weber/1000G_phased/1k_genome_EAS_phased_chr${CHR}.vcf.gz \
              --sequencing
        done < <(cat phasing/chr${CHR}.chunks.txt)
    done
done

mkdir -p shapeit_chr2/group
mkdir -p shapeit_chr2/complete
# Step1: concat 每個 group 的 chunks
for CHR in {1..22} X; do
    echo "Concatenating groups in chr${CHR}..."
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue
        group=$(basename "$groupfile" .vcf.gz)
        echo "concat ${group}"
        # bcftools concat -Oz -o shapeit_chr2/group/merged_${group}.vcf.gz shapeit_chr2/chr${group}_*.phased.vcf.gz
        bcftools concat -Oz -o shapeit_chr2/group/merged_${group}.vcf.gz \
  $(awk -v grp="$group" '{print "shapeit_chr2/" grp "_" $1 ".phased.vcf.gz"}' phasing/chr${CHR}.chunks.txt)
        bcftools index -c shapeit_chr2/group/merged_${group}.vcf.gz
    done
done


# Step2: merge 每條 chr 的所有 group
for CHR in {1..22} X; do
    echo "Merging chr${CHR}..."
    bcftools merge -Oz -o shapeit_chr2/complete/merged_chr${CHR}.vcf.gz \
        shapeit_chr2/group/merged_${CHR}_samples_group_*.vcf.gz
    bcftools index -c shapeit_chr2/complete/merged_chr${CHR}.vcf.gz
done

mkdir -p shapeit_chr/group
mkdir -p shapeit_chr/complete
for CHR in {1..22}; do
    echo "Concatenating groups in chr${CHR}..."
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue
        group=$(basename "$groupfile" .vcf.gz)
        echo "concat ${group}"
        # bcftools concat -Oz -o shapeit_chr/group/merged_${group}.vcf.gz shapeit_chr2/chr${group}_*.phased.vcf.gz
        bcftools concat -Oz -o shapeit_chr/group/merged_${group}.vcf.gz \
  $(awk -v grp="$group" '{print "shapeit_chr/" grp "_" $1 ".phased.vcf.gz"}' phasing/chr${CHR}.chunks.txt)
        bcftools index -c shapeit_chr/group/merged_${group}.vcf.gz
    done
done


# Step2: merge 每條 chr 的所有 group
for CHR in {1..22}; do
    echo "Merging chr${CHR}..."
    bcftools merge -Oz -o shapeit_chr/complete/merged_chr${CHR}.vcf.gz \
        shapeit_chr/group/merged_${CHR}_samples_group_*.vcf.gz
    bcftools index -c shapeit_chr/complete/merged_chr${CHR}.vcf.gz
done
'''
'''
for CHR in {1..22} X; do
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
            shapeit4 \
              --input group_phasing/${group}.vcf.gz \
              --map ${MAP} \
              --region ${REGION} \
              --output shapeit_chr3/chr${group}_${REGION}.phased.vcf.gz \
              --thread 1 \
              --pbwt-depth 1 \
              --window 2.0 \
              --reference /home/Weber/1000G_phased/1k_genome_EAS_phased_chr${CHR}.vcf.gz \
              --sequencing
        done < <(cat phasing/chr${CHR}.chunks.txt)
    done
done

mkdir -p shapeit_chr3/group
mkdir -p shapeit_chr3/complete
# Step1: concat 每個 group 的 chunks
for CHR in {1..22} X; do
    echo "Concatenating groups in chr${CHR}..."
    for groupfile in group_phasing/${CHR}_samples_group_*; do
        [[ $groupfile == *.csi ]] && continue
        group=$(basename "$groupfile" .vcf.gz)
        echo "concat ${group}"
        bcftools concat -Oz -o shapeit_chr3/group/merged_${group}.vcf.gz \
  $(awk -v grp="$group" '{print "shapeit_chr3/chr" grp "_" $1 ".phased.vcf.gz"}' phasing/chr${CHR}.chunks.txt)
        bcftools index -c shapeit_chr3/group/merged_${group}.vcf.gz
    done
done


# Step2: merge 每條 chr 的所有 group
for CHR in {1..22} X; do
    echo "Merging chr${CHR}..."
    bcftools merge --force-samples -Oz -o shapeit_chr3/complete/merged_chr${CHR}.vcf.gz \
        shapeit_chr3/group/merged_${CHR}_samples_group_*.vcf.gz
    bcftools index -c shapeit_chr3/complete/merged_chr${CHR}.vcf.gz
done

# sample 確認
for CHR in {1..22} X; do
    bcftools view \
      -S target_samples.txt \
      -Oz -o shapeit_chr3/complete/merged_chr${CHR}.target.vcf.gz \
      shapeit_chr3/complete/merged_chr${CHR}.vcf.gz
    
    bcftools index -c shapeit_chr3/complete/merged_chr${CHR}.target.vcf.gz

    echo -n "chr${CHR}: "
    bcftools view -H shapeit_chr3/complete/merged_chr${CHR}.target.vcf.gz | wc -l
done
'''
mkdir -p shapeit_chr4/logs
for CHR in {1..6}; do
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

for CHR in {7..13}; do
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

for CHR in {13..22}; do
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

mkdir -p shapeit_chr4/group
mkdir -p shapeit_chr4/complete
# Step1: concat 每個 group 的 chunks
for CHR in {1..22}; do
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
for CHR in {1..22}; do
    echo "Merging chr${CHR}..."
    bcftools merge --force-samples -Oz -o shapeit_chr4/complete/merged_chr${CHR}.vcf.gz \
        shapeit_chr4/group/merged_${CHR}_samples_group_*.vcf.gz
    bcftools index -c shapeit_chr4/complete/merged_chr${CHR}.vcf.gz
done

# sample 確認
for CHR in {1..22}; do
    bcftools view \
      -S target_samples.txt \
      -Oz -o shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz \
      shapeit_chr4/complete/merged_chr${CHR}.vcf.gz
    
    bcftools index -c shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz

    echo -n "chr${CHR}: "
    bcftools view -H shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz | wc -l
done
for CHR in {1..22}; do
    bcftools query -l shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz | wc -l
done

# chrx
mkdir -p shapeit_chrx/chrX
mkdir -p shapeit_chrx/chrX_par1
mkdir -p shapeit_chrx/chrX_par2
mkdir -p shapeit_chrx/complete
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
            
            # 每累積超過 9 cM 就輸出
            if ($3 - start_cM >= 8) {
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
    ' > phasing/${REGION}.chunks.txt
done


MAP_DIR=/home/Weber/hg38_chr

for REGION_TYPE in chrX_par1 chrX_par2; do
    MAP=/home/Weber/hg38_chr/chr${REGION_TYPE}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Phasing ${REGION_TYPE}..."
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l)
        echo " SNP count: $lines "
        shapeit4 \
          --input cleaned_chrX_merged_singletonremoved2.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_${REGION}.phased.vcf.gz \
          --thread 10 \
          --pbwt-depth 1 \
          --window 3.0 \
          --sequencing \
          --log shapeit_chr4/logs/${REGION}.log
    done < phasing/${REGION_TYPE}.chunks.txt
done

# cleaned_chrX_merged_singletonremoved2 的 par2 沒有snp

for REGION_TYPE in chrX; do
    MAP=/home/Weber/hg38_chr/chr${REGION_TYPE}.b38.gmap.gz.chrprefix.b38.gmap.gz
    echo "Phasing ${REGION_TYPE}..."
    while read REGION; do
        echo "  REGION ${REGION}"
        lines=$(bcftools view -r ${REGION} -H cleaned_chrX_merged_singletonremoved2_female.vcf.gz | wc -l)
        echo " SNP count: $lines "
        shapeit4 \
          --input cleaned_chrX_merged_singletonremoved2_female.vcf.gz \
          --map ${MAP} \
          --region ${REGION} \
          --output shapeit_chrx/${REGION_TYPE}/${REGION_TYPE}_${REGION}.phased.vcf.gz \
          --thread 10 \
          --pbwt-depth 1 \
          --window 3.0 \
          --sequencing \
          --log shapeit_chr4/logs/${REGION}.log
    done < phasing/${REGION_TYPE}.chunks.txt
done

for REGION_TYPE in chrX chrX_par1 chrX_par2; do
    echo "Concatenating in ${REGION_TYPE}..."
    bcftools concat -Oz -o shapeit_chrx/complete/${REGION_TYPE}_1.vcf.gz \
  $(awk -v grp="$REGION_TYPE" '{print "shapeit_chrx/" grp "/" grp "_" $1 ".phased.vcf.gz"}' phasing/${REGION_TYPE}.chunks.txt)
    bcftools index -c shapeit_chrx/complete/${REGION_TYPE}_1.vcf.gz
done


# nonpar_male
bcftools view -S male_samples.txt cleaned_chrX_merged_singletonremoved2.vcf.gz \
  | sed 's/\/|/|/g' | sed 's/\//|/g' \
  | bgzip -c > cleaned_chrX_merged_singletonremoved2_male.vcf.gz
tabix -p vcf cleaned_chrX_merged_singletonremoved2_male.vcf.gz

bcftools view -r chrX:2785350-155699751 cleaned_chrX_merged_singletonremoved2_male.vcf.gz \
| awk 'BEGIN{OFS="\t"}
       /^#/ {
         if($0 ~ /^#CHROM/) {
           # 重新印出 header
           for(i=1; i<=NF; i++) {
             if(i <= 8) printf "%s\t", $i;
             else if(i == 9) printf "FORMAT";
             else printf "\t%s", $i;
           }
           print "";
         } else print;
         next;
       }
       {
         # 若無 ID，補上 chrX:POS:REF:ALT
         if($3 == "." || $3 == "") $3=$1":"$2":"$4":"$5;
         # 印出前 8 欄與 INFO
         printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT", $1,$2,$3,$4,$5,$6,$7,$8;
         # sample 部分
         for(i=10; i<=NF; i++) {
           gt=$i;
           sub(/:.*/,"",gt);      # 移除冒號後資訊
           gsub(/[\/|]/,"|",gt);  # 統一成 |
           printf "\t%s", gt;
         }
         print "";
       }' \
 | bgzip -c > shapeit_chrx/complete/chrX_2.vcf.gz
tabix -p vcf shapeit_chrx/complete/chrX_2.vcf.gz

# merge nonpar
bcftools merge -Oz -o shapeit_chrx/complete/chrX_3.vcf.gz \
  shapeit_chrx/complete/chrX_1.vcf.gz \
  shapeit_chrx/complete/chrX_2.vcf.gz
bcftools index -c shapeit_chrx/complete/chrX_3.vcf.gz

bcftools query -l shapeit_chrx/complete/chrX_par1_1.vcf.gz > chrX_order.txt
bcftools view -S chrX_order.txt -Oz -o shapeit_chrx/complete/chrX_3.sorted.vcf.gz shapeit_chrx/complete/chrX_3.vcf.gz
bcftools index -c shapeit_chrx/complete/chrX_3.sorted.vcf.gz

# 修正
{ bcftools view -h shapeit_chrx/complete/chrX_3.sorted.vcf.gz;
  bcftools view -H shapeit_chrx/complete/chrX_3.sorted.vcf.gz \
  | awk 'BEGIN{OFS="\t"}
         {
           printf "%s", $1;
           for(i=2;i<=9;i++) printf "\t%s",$i;
           for(i=10;i<=NF;i++){
             gt=$i;
             if(gt=="0") gt="0|0";
             else if(gt=="1") gt="1|1";
             printf "\t%s", gt;
           }
           print "";
         }'
} | bgzip -c > shapeit_chrx/complete/chrX_3.phased.vcf.gz
bcftools index -c shapeit_chrx/complete/chrX_3.phased.vcf.gz


# concat
bcftools concat -Oz -o shapeit_chr4/complete/merged_chrX.target.vcf.gz \
  shapeit_chrx/complete/chrX_par1_1.vcf.gz \
  shapeit_chrx/complete/chrX_3.phased.vcf.gz 
bcftools index -c shapeit_chr4/complete/merged_chrX.target.vcf.gz

for CHR in X; do
    bcftools view -H shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz | wc -l
    bcftools query -l shapeit_chr4/complete/merged_chr${CHR}.target.vcf.gz | wc -l
done

bcftools view -H cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l
bcftools query -l cleaned_chrX_merged_singletonremoved2.vcf.gz | wc -l



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
"""