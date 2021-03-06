#######################################################################
#               EPIGENETIC CLOCKS AND MULTIPLE CANCERS                #
#######################################################################
# R version 4.0.2 (2020-06-22)
# Last modified: 17 August 2021

# This script runs two-sample MR for epigenetic clock acceleration
# measures and cancer in FinnGen

#---------------------------------------------------------------------#
#                            Housekeeping                             #----
#---------------------------------------------------------------------#

# Clear environment
rm(list=ls()) #Remove any existing objects in R 

# Set working directory 
setwd("your_working_directory") 

if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot", "openxlsx")

#---------------------------------------------------------------------#
#                            Read exposure                             #----
#---------------------------------------------------------------------#
# Read exposure dataset
exp_dat <- read.table("exp_data.txt", header = T)

#---------------------------------------------------------------------#
#                            Outcomes                                  #----
#---------------------------------------------------------------------#

# Upload outcome datasets including proxies (these are already in two-sample MR format)
breast_cancer_excall_proxies <- read.table("finngen_breast_exallc_replaced_proxies.txt", header = T)
ovarian_cancer_excall_proxies <- read.table("finngen_ovarian_exallc_replaced_proxies.txt", header = T)
prostate_cancer_excall_proxies <- read.table("finngen_prostate_exallc_replaced_proxies.txt", header = T)
lung_cancer_excall_proxies <- read.table("finngen_lung_exallc_replaced_proxies.txt", header = T)
colorectal_cancer_excall_proxies <- read.table("finngen_colorectal_exallc_replaced_proxies.txt", header = T)
#读取已经format好的finngen5个癌种位点信息，合并
FINNGEN_out_dat <- rbind(breast_cancer_excall_proxies, ovarian_cancer_excall_proxies,
                         prostate_cancer_excall_proxies, lung_cancer_excall_proxies,
                         colorectal_cancer_excall_proxies)

#---------------------------------------------------------------------#
#                            Harmonisation                             #----
#---------------------------------------------------------------------#
# Harmonise exposure and outcome datasets #根据allele和频率判断两个gwas结果中的SNP位点是否一致，不一致的进行去除，或进行调整来协同统一。将效应等位基因（effect allele）对齐。
data_FINNGEN <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = FINNGEN_out_dat
)
#---------------------------------------------------------------------#
#                            Results                                   #----
#---------------------------------------------------------------------#
#MR
results_FINNGEN <- mr(data_FINNGEN, method_list = c("mr_ivw", "mr_egger_regression", 
                                        "mr_weighted_median", "mr_weighted_mode")) 用四种方法进行计算MR结果

#---------------------------------------------------------------------#
#                            Save results                             #----
#---------------------------------------------------------------------#

#write.table(results_FINNGEN, "results_FINNGEN.txt", row.names = F)
