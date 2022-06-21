
#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 26 August 2021

# This script runs two-sample MR for epigenetic clock acceleration
# measures and cancer in consortiums

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools")

#---------------------------------------------------------------------#
#                            Available data                           #----
#---------------------------------------------------------------------#

ao <- available_outcomes() # 读取MR-Base平台上的GWAS数据


#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#
#read exposure
exp_dat <- read.table("exp_data.txt", header = T) #读取exposure数据

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#
# Extract outcomes from the IEU catalog (MR-Base) #读取outcome特定癌种数据中与exposureSNPmatch的数据
IEU_out_dat <- extract_outcome_data(
  snps = exp_dat$SNP,
  outcomes = c("ieu-a-1126", "ieu-b-85", "ieu-a-966", "ieu-a-1120", "ieu-a-1127", 
               "ieu-a-1128", "ieu-a-1121", "ieu-a-1122", "ieu-a-1123", "ieu-a-1124",
               "ieu-a-1125", "ieu-a-965", "ieu-a-967")
) 

# Read GECCO outcomes #读取GECCO数据
x <- read.csv("GECCO.csv", 
                header = T)
x$SNP <- exp_dat$SNP[match(x$Position, exp_dat$pos.exposure)]


out_func_GECCO <- function(name) #对GECCO数据进行调整format，读取outcome数据中与exposureSNPmatch的数据
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_var <- x %>% format_data(snps = exp_dat$SNP, 
                          type = "outcome", 
                          snp_col = "SNP", 
                          beta_col = "Effect", 
                          se_col = "StdErr", 
                          effect_allele_col = "Allele1", 
                          other_allele_col = "Allele2",
                          chr_col = "Chr",
                          pos_col = "Position",
                          eaf_col = "Freq1", 
                          pval_col = "P.value")
  outcome_var$outcome <- name
  return(outcome_var)
}

GECCO <- out_func_GECCO("Colorectal cancer")

IEU_out_dat <- smartbind(IEU_out_dat, GECCO) #合并
#---------------------------------------------------------------------#
#                            Harmonisation                             #----
#---------------------------------------------------------------------#
# Harmonise exposure and outcome datasets #调整正反链，使它们方向一致。
data_IEU <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = IEU_out_dat
)

#---------------------------------------------------------------------#
#                            Results                                   #----
#---------------------------------------------------------------------#
# Run two-sample MR 执行MR分析，用四种计算方法
results_IEU <- mr(data_IEU, method_list = c("mr_ivw", "mr_egger_regression", 
                                    "mr_weighted_median", "mr_weighted_mode"))

#---------------------------------------------------------------------#
#                            Save results                             #----
#---------------------------------------------------------------------#

#write.table(results_IEU, "results_IEU.txt", row.names = F)
