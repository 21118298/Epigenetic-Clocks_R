一、Finngen并没有在MR-Base platform，因此需要用‘LDlinkR’ R package to find proxies。

#finngen dataset过大，运行过慢。所以每个癌种按染色体拆分，只留下chr6，chr10，chr17三条相关染色体文件合并。

x <- read.table("finngen_LUNG_EXALLC.txt",col.names = c("chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta","sebeta","maf","maf_cases","maf_controls"), sep = "\t")
result <- split(x, x$chrom)
for (i in c(6,10,17)){                  
   write.table(result[[i]], file = paste0("lung_chr",i,sep=""), quote = FALSE, row.names = FALSE)
}
newdata <- rbind(result[[6]],result[[10]],result[[17]])
write.table(nenwdata, file = paste0("finngen_lung_filter.txt"), quote = FALSE, row.names = FALSE)

结果意义：利用LDlinkR-LDproxy找到Exposore GrimAge4个SNP在所有Finngen5个癌种的GWASdataset都存在的proxy。(5个癌种breast,prostate, colorectal, ovarian and lung cancers)

结果如下: GrimAge的4个snp在finngen dataset中都能找到不需要proxy。
  chr.outcome pos.outcome other_allele.outcome effect_allele.outcome        SNP pval.outcome beta.outcome se.outcome eaf.outcome
1           6    31175734                    G                     A   rs887466    0.0482259   -0.0522429  0.0264471    0.472498
2           6   109297501                    C                     T  rs9386796    0.2056880    0.0334395  0.0264237    0.505749
3          10    99520522                    A                     G rs17094148    0.8839680    0.0039362  0.0269711    0.392376
4          17    39987295                    C                     T  rs4065321    0.1607330    0.0370189  0.0263929    0.480417
                                     outcome mr_keep.outcome pval_origin.outcome id.outcome
1 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
2 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
3 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
4 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA

尽管如此，还是运行一下该步骤的代码。结果如下：

