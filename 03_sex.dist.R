#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
wkdir <- args[1]
bfile <- args[2]
f_low_thres <- args[3]
f_high_thres <- args[4]

# bfile <- "mafcallrateMono_qc"
# f_low_thres <- 0.2
# f_high_thres <- 0.8

setwd(wkdir)
library(ggplot2)
library(dplyr)

# Read in data
imp_sex <- read.table(paste0(bfile,".sexcheck"), header=FALSE)
colnames(imp_sex) <- c("IID", "PEDSEX", "SNPSEX", "STATUS", "FSTAT")
head(imp_sex)
imp_sex$F <- as.numeric(as.character(imp_sex$FSTAT))


sex_info <- read.table("sexinfo_plink2.txt", header=FALSE, col.names=c("IID","PEDSEX"))
head(sex_info)
imp_sex <- merge(imp_sex, sex_info, by="IID", all.x=TRUE, suffixes=c("", "_sexinfo"))
head(imp_sex)

# 使用 sexinfo 來更新 PEDSEX
imp_sex$PEDSEX <- factor(imp_sex$PEDSEX_sexinfo, levels=c(0,1,2), labels=c("Unknown","Male","Female"))

# SNPSEX 也轉成 factor
imp_sex$SNPSEX <- factor(imp_sex$SNPSEX, levels=c(0,1,2), labels=c("Unknown","Male","Female"))

head(imp_sex)

# 用 PEDSEX vs SNPSEX 判斷 Fail
imp_sex$SEXCHECK_FAIL <- ifelse(
  is.na(imp_sex$SNPSEX) | (as.character(imp_sex$PEDSEX) != as.character(imp_sex$SNPSEX)),
  "FAIL", "PASS"
)

# 檢查結果
table(imp_sex$SEXCHECK_FAIL)

# 存完整表格
write.table(imp_sex, file="chrX_sexcheck_annotated.txt",
            quote=FALSE, row.names=FALSE, sep="\t")

# 輸出 fail 的 sample ID
fail_ids <- imp_sex$IID[imp_sex$SEXCHECK_FAIL == "FAIL"]
write.table(fail_ids, file="chrX_sexcheck_failed_samples.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#check the sample size
imp_sex %>% group_by(PEDSEX, SNPSEX) %>% tally() %>% data.frame()


# Plot F-stat distribution: histogram
p <- ggplot(imp_sex, aes(x=F, fill=SNPSEX)) +
  geom_histogram(alpha=0.6, binwidth=0.02, color="black") +
  geom_vline(xintercept=as.numeric(f_high_thres), linetype="dashed") +  #default cutoff: >0.8: male
  geom_vline(xintercept=as.numeric(f_low_thres), linetype="dashed") +  #default cutoff: <0.2: female
  theme_bw() +
  labs(x='F-statistic', y='Frequency', title='sample QC', fill='Imputed sex') +
  scale_fill_manual(
    values = c("grey30", "#56B4E9", "red"),
    breaks = c("Unknown","Male","Female"),
    labels = c("Unknown","Male","Female")
  )
ggsave(paste0(bfile,"_sexcheck_fstat_F",f_low_thres,"_M",f_high_thres,".tiff"), p, width=7.5, height=4)
 
 
# Plot: imputed vs. reported gender
tiff(paste0(bfile,"_sexcheck_fstat_F",f_low_thres,"_M",f_high_thres,"_ped.vs.snpsex.imputed.tiff"), width=7.5, height=5, units = "in", res = 300)
ggplot(imp_sex, aes(y=F,x=factor(SNPSEX), color=factor(SNPSEX))) +
  geom_jitter(alpha=0.7) +
  geom_hline(yintercept=as.numeric(f_high_thres), linetype="dashed") + #default cutoff: >0.8: male
  geom_hline(yintercept=as.numeric(f_low_thres), linetype="dashed") + #default cutoff: <0.2: female
  labs(x='Imputed gender',y='chrX F-statistic',
       title=paste0('sample QC'),
       color='Imputed gender') +
  scale_color_manual(
  values = c("gray", "#56B4E9", "red"),
  breaks = c("Unknown","Male","Female"),
  labels = c("Unknown","Male","Female")
) +
  theme_bw()
dev.off()


tiff(paste0(bfile,"_sexcheck_fstat_F",f_low_thres,"_M",f_high_thres,"_ped.vs.snpsex.reported.tiff"), width=7.5, height=5, units = "in", res = 300)
ggplot(imp_sex, aes(y=F,x=factor(SNPSEX), color=factor(PEDSEX))) +
  geom_jitter(alpha=0.7) +
  geom_hline(yintercept=as.numeric(f_high_thres), linetype="dashed") + #default cutoff: >0.8: male
  geom_hline(yintercept=as.numeric(f_low_thres), linetype="dashed") + #default cutoff: <0.2: female
  labs(x='Imputed gender',y='chrX F-statistic',
       title=paste0('sample QC'),
       color='Reported gender') +
  scale_color_manual(
  values = c("gray", "#56B4E9", "red"),
  breaks = c("Unknown","Male","Female"),
  labels = c("Unknown","Male","Female")
) +
  theme_bw()
dev.off()
 
 
# # Write out a list of sex discordant samples for removal: given a few differnet thresholds
# sex_mismatch_F02_M08_inds <- imp_sex[(imp_sex$STATUS=="PROBLEM"), c("FID","IID")]
 
# imp_sex$PEDSEX <- as.character(imp_sex$PEDSEX)
# imp_sex$SNPSEX <- as.character(imp_sex$SNPSEX)
# imp_sex$SNPSEX <- ifelse(imp_sex$F > 0.75, "Male", imp_sex$SNPSEX)
# imp_sex$SNPSEX <- ifelse(imp_sex$F < 0.25, "Female", imp_sex$SNPSEX)
# sex_mismatch_F025_M075_inds <- imp_sex[imp_sex$SNPSEX != imp_sex$PEDSEX, c("FID","IID")]
# #dim(sex_mismatch_F025_M075_inds)
 
# imp_sex$SNPSEX <- ifelse(imp_sex$F > 0.7, "Male", imp_sex$SNPSEX)
# imp_sex$SNPSEX <- ifelse(imp_sex$F < 0.3, "Female", imp_sex$SNPSEX)
# sex_mismatch_F03_M07_inds <- imp_sex[imp_sex$SNPSEX != imp_sex$PEDSEX, c("FID","IID")]
# #dim(sex_mismatch_F03_M07_inds)
 
# imp_sex$SNPSEX <- ifelse(imp_sex$F > 0.6, "Male", imp_sex$SNPSEX)
# imp_sex$SNPSEX <- ifelse(imp_sex$F < 0.4, "Female", imp_sex$SNPSEX)
# sex_mismatch_F04_M06_inds <- imp_sex[imp_sex$SNPSEX != imp_sex$PEDSEX, c("FID","IID")]
# #dim(sex_mismatch_F04_M06_inds)
 
# imp_sex$SNPSEX <- ifelse(imp_sex$F > 0.5, "Male", imp_sex$SNPSEX)
# imp_sex$SNPSEX <- ifelse(imp_sex$F < 0.5, "Female", imp_sex$SNPSEX)
# sex_mismatch_F05_M05_inds <- imp_sex[imp_sex$SNPSEX != imp_sex$PEDSEX, c("FID","IID")]
# dim(sex_mismatch_F05_M05_inds)
 
 
# write.table(sex_mismatch_F02_M08_inds, paste0(bfile,"_sex_mismatch_F02_M08.indlist"), quote=F, col.names=F, row.names=F)
# write.table(sex_mismatch_F025_M075_inds, paste0(bfile,"_sex_mismatch_F025_M075.indlist"), quote=F, col.names=F, row.names=F)
# write.table(sex_mismatch_F03_M07_inds, paste0(bfile,"_sex_mismatch_F03_M07.indlist"), quote=F, col.names=F, row.names=F)
# write.table(sex_mismatch_F04_M06_inds, paste0(bfile,"_sex_mismatch_F04_M06.indlist"), quote=F, col.names=F, row.names=F)
# write.table(sex_mismatch_F05_M05_inds, paste0(bfile,"_sex_mismatch_F05_M05.indlist"), quote=F, col.names=F, row.names=F)
