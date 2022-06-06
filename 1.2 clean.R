rm(list=ls())
setwd("E:/")
if (!require("pacman")) install.packages("pacman")
pacman::p_load("MRInstruments", "TwoSampleMR", "tidyverse", "dplyr", "ggpubr", "ggplot2", "ggforce", "data.table", "ggforestplot","gtools")
#exposure_func <- function(data_name, file, SNP) {
  x <- read.table(file = file, header = T)  
  x$pheno <- "DNA methylation ageing"     
  x <- format_data(x, type = "exposure",   
                   phenotype_col = "pheno",    
                   snp_col = "rsID",              
                   beta_col = "Effect",        
                   se_col = "SE",              
                   eaf_col = "Freq1",          
                   effect_allele_col = "A1",   
                   other_allele_col = "A2",     
                   pval_col = "P",              
                   samplesize_col = "N",        
                   chr_col = "chr",             
                   pos_col = "bp")              
  x$id.exposure <- data_name
  re <- split(x,x$chr.exposure)
  y <- clump_data(re[[22]], clump_p1 = 5e-08, clump_p2 = 5e-08)
  
  for (i in 1:22){
  x <- rbind(assign(paste0(file,"_chr",i,sep=""),clump_data(result[[i]], clump_p1 = 5e-08, clump_p2 = 5e-08))
  #write.table(result[[i]], file = paste0(file,"_chr",i,sep=""), quote = FALSE)
  #y <- clump_data(result[[i]], clump_p1 = 5e-08, clump_p2 = 5e-08)
  #write.table(y, file = paste0(file,"_chr",i,sep=""), quote = FALSE)
  }
