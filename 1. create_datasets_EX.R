#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script creates exposure datasets for epigenetic clock acceleration
# measures, which can later be used in two-sample MR analyses
# This script also estimates R2 and F-statistics for each of the exposures

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) ##删除目前所有的变量

# Set working directory 
setwd("your_working_directory") #设置当前工作目录

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools")
#批量安装和加载R包。判断环境中是否有这个包，如果没有，先安装再加载；如果有，直接加载。
#---------------------------------------------------------------------#
#                           Exposures                                  #----
#---------------------------------------------------------------------#
# Create function to extract GWAS significant SNPs and perform LD-clumping 
# Output in two-sample MR format
# NOTE: this takes really long if run in RStudio
exposure_func <- function(data_name, file, SNP) {   #自定义函数exposure_func
  x <- read.table(file = file, header = T)  # 读取Exposures GWAS txt文件，第一行为变量名
  x$pheno <- "DNA methylation ageing"      #给pheno整列赋值
  x <- format_data(x, type = "exposure",   #调用TwoSampleMR R包中format_date函数，将输入文件格式进行调整，以适用于MR分析
                   phenotype_col = "pheno",    #定义phenotype_col在输入文件中的对应列
                   snp_col = SNP,              #定义snp_col在输入文件中的对应列，SNP由参数传入，*必要
                   beta_col = "Effect",        #同上,The effect size. If the trait is binary then log(OR) should be used*必要
                   se_col = "SE",              #The standard error of the effect size*必要
                   eaf_col = "Freq1",          #The effect allele frequency
                   effect_allele_col = "A1",   #The allele of the SNP which has the effect marked in beta*必要
                   other_allele_col = "A2",     #The non-effect allele
                   pval_col = "P",              #The P-value for the SNP’s association with the exposure
                   samplesize_col = "N",        #Sample size for estimating the effect size
                   chr_col = "chr",             #Physical position of variant (chromosome)
                   pos_col = "bp")              #Physical position of variant (position)
  x$id.exposure <- data_name                              #在最后一列加入Clock Age标签
  x <- clump_data(x, clump_p1 = 5e-08, clump_p2 = 5e-08)  # 调用TwoSampleMR R包中clump_data函数，用PLINK clumping法, 识别和保留每个LD块中最重要的SNP（最低p值), 
  # p-value < 5*10-8, r2 <0.001(default), 该函数与OpenGWAS API进行交互，存储了千人基因组中5个群体（EUR, SAS, EAS, AFR, AMR）的LD数据。Defult = EUR(European reference)
}
# Apply function to raw epigenetic age acceleration datasets
GrimAge_exp_dat <- exposure_func("GrimAge","GrimAge_EUR_summary_statistics.txt", "rsID") #执行GrimAge
Hannum_exp_dat <- exposure_func("Hannum","Hannum_EUR_summary_statistics.txt", "roblrsID") #执行Hannum
IEAA_exp_dat <- exposure_func("IEAA","IEAA_EUR_summary_statistics.txt", "rsID")           #执行IEAA
PhenoAge_exp_dat <- exposure_func("PhenoAge","PhenoAge_EUR_summary_statistics.txt", "rsID") #执行PhenoAge

#combine exposures
exp_dat <- rbind(GrimAge_exp_dat, Hannum_exp_dat, IEAA_exp_dat, PhenoAge_exp_dat) #合并结果，HannumAge (9 SNPs)), Horvath Intrinsic Age (24 SNPs), PhenoAge (11 SNPs), and GrimAge (4 SNPs)

#save unique list of SNPs
#write.table(unique(exp_dat$SNP), "SNP_list.txt", row.names = F, col.names = F)

#---------------------------------------------------------------------#
#                          R2 and F-statistic                         #----显著性检验，用来检验结果精密度和偶然误差
#---------------------------------------------------------------------#

# Calculate R2 and F statistics for each exposure dataset（R2检验和F检验）
#method 1（R2检验和F检验的计算公式）
exp_dat$r2 <- (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure)) /
  (2 * (exp_dat$beta.exposure^2) * exp_dat$eaf.exposure * (1 - exp_dat$eaf.exposure) +
     2 * exp_dat$samplesize.exposure * exp_dat$eaf.exposure * 
     (1 - exp_dat$eaf.exposure) * exp_dat$se.exposure^2) 
#
#R2 = (2β2×MAF×(1-MAF))/(2β2×MAF×(1-MAF)+ 2 N × MAF × (1-MAF) ×SE2) where MAF = effect allele frequency, β = effect estimate of the SNP in the exposure GWAS,SE = standard error, N = sample size
#
exp_dat$F <- exp_dat$r2 * (exp_dat$samplesize.exposure - 2) / (1 - exp_dat$r2)
#
#F = (R2×(N-2))/(1-R2) 
#
#method 2
# exp_dat$F_stat <- exp_dat$beta.exposure^2 / exp_dat$se.exposure^2
# exp_dat$R2_stat <- exp_dat$F_stat/(exp_dat$samplesize.exposure-2+exp_dat$F_stat)

# Calculate total R2 for each exposure dataset （计算每个clock的R2之和）
r2_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  sum(x$r2, na.rm = T)
}

variance_GrimAge <- r2_func("GrimAge") # GrimAgeSNP R2之和为 0.47%
variance_Hannum <- r2_func("Hannum") # Hannum R2之和为1.48%
variance_IEAA <- r2_func("IEAA") # IEAASNP R2之和为4.41%
variance_PhenoAge <- r2_func("PhenoAge") #PhenoAgeSNP R2之和为1.86%

# Calculate minimum F-statistic for each exposure dataset （计算每个clock的F检验最小值）
Fmin_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  min(x$F, na.rm = T)
}

Fmin_GrimAge <- Fmin_func("GrimAge") # 31
Fmin_Hannum <- Fmin_func("Hannum") # 31
Fmin_IEAA <- Fmin_func("IEAA") # 31
Fmin_PhenoAge <- Fmin_func("PhenoAge") # 32 

# Calculate maximum F-statistic for each exposure dataset （计算每个clock的F检验最大值）
Fmax_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  max(x$F, na.rm = T)
}

Fmax_GrimAge <- Fmax_func("GrimAge") # 45
Fmax_Hannum <- Fmax_func("Hannum") # 99
Fmax_IEAA <- Fmax_func("IEAA") # 240
Fmax_PhenoAge <- Fmax_func("PhenoAge") # 89 

# Calculate median F-statistic for each exposure dataset （计算每个clock的F检验中位数）
Fmedian_func <- function(id)
{
  x <- exp_dat[which(exp_dat$id.exposure==id),]
  median(x$F, na.rm = T)
}

Fmedian_GrimAge <- Fmedian_func("GrimAge") # 36
Fmedian_Hannum <- Fmedian_func("Hannum") # 38
Fmedian_IEAA <- Fmedian_func("IEAA") # 47
Fmedian_PhenoAge <- Fmedian_func("PhenoAge") # 45 

#save
#write.table(exp_dat, "exp_data.txt", row.names = F)
