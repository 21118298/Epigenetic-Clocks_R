一、Finngen并没有在MR-Base platform，因此需要用‘LDlinkR’ R package to find proxies。

#finngen dataset过大，运行过慢。所以每个癌种按染色体拆分，只留下chr6，chr10，chr17三条相关染色体文件合并。

x <- read.table("finngen_LUNG_EXALLC.txt",col.names = c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","maf","maf_cases","maf_controls"), sep = "\t")
for (i in c(6,10,17)){                  
   result <- split(x, x$chrom)
   write.table(result[[i]], file = paste0("lung_chr",i,sep=""), quote = FALSE, row.names = FALSE)
}

结果意义：利用LDlinkR-LDproxy找到Exposore GrimAge4个SNP在所有Finngen5个癌种的GWASdataset都存在的proxy。(5个癌种breast,prostate, colorectal, ovarian and lung cancers)

结果如下1，3，5，7为original，2，4，6，8为proxy：
    RS_Number           Coord Alleles    MAF Distance Dprime     R2 Correlated_Alleles RegulomeDB Function   original
1   rs4065321  chr17:38143548   (C/T) 0.4642        0      1 1.0000            C=C,T=T          5     <NA>  rs4065321
4  rs12453764  chr17:38144187   (T/C) 0.4642      639      1 1.0000            C=T,T=C          7     <NA>  rs4065321
11  rs9386796  chr6:109618704   (C/T) 0.4284        0      1 1.0000            C=C,T=T         1f     <NA>  rs9386796
2   rs5005289  chr6:109617673   (C/G) 0.4284    -1031      1 1.0000            C=C,T=G          4     <NA>  rs9386796
12   rs887466   chr6:31143511   (G/A) 0.4076        0      1 1.0000            G=G,A=A         1b     <NA>   rs887466
21  rs1265155   chr6:31143694   (G/A) 0.3877      183      1 0.9203            G=G,A=A         2b     <NA>   rs887466
13 rs17094148 chr10:101280279   (A/G) 0.3022        0      1 1.0000            A=A,G=G          4     <NA> rs17094148
22 rs11190133 chr10:101278725   (C/T) 0.3022    -1554      1 1.0000            A=C,G=T          5     <NA> rs17094148

