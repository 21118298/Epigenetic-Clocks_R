Finngen并没有在MR-Base platform，因此需要用‘LDlinkR’ R package to find proxies。

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
  chr.outcome pos.outcome other_allele.outcome effect_allele.outcome        SNP pval.outcome beta.outcome se.outcome eaf.outcome outcome mr_keep.outcome pval_origin.outcome id.outcome
1           6    31175734                    G                     A   rs887466    0.0482259   -0.0522429  0.0264471    0.472498 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
2           6   109297501                    C                     T  rs9386796    0.2056880    0.0334395  0.0264237    0.505749 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
3          10    99520522                    A                     G rs17094148    0.8839680    0.0039362  0.0269711    0.392376 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA
4          17    39987295                    C                     T  rs4065321    0.1607330    0.0370189  0.0263929    0.480417 lung cancer (excluding cancer in controls)            TRUE            reported     zIV0uA


尽管如此，还是运行一下该步骤的代码。1-4是exp_dat对应的snp，5-8是proxy的format_data后的结果。结果如下：

  chr.outcome pos.outcome other_allele.outcome effect_allele.outcome        SNP pval.outcome beta.outcome se.outcome eaf.outcome outcome mr_keep.outcome pval_origin.outcome id.outcome
1           6    31175734                    G                     A   rs887466    0.0482259  -0.05224290  0.0264471    0.472498 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
2           6   109297501                    C                     T  rs9386796    0.2056880   0.03343950  0.0264237    0.505749 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
3          10    99520522                    A                     G rs17094148    0.8839680   0.00393620  0.0269711    0.392376 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
4          17    39987295                    C                     T  rs4065321    0.1607330   0.03701890  0.0263929    0.480417 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
5          17    39987515                    C                     T rs12940405    0.1598920   0.03709510  0.0263940    0.480414 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
6           6   109296470                    C                     G  rs5005289    0.2052730   0.03347070  0.0264242    0.505742 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
7           6    31175917                    G                     A  rs1265155    0.0441448  -0.05326680  0.0264652    0.468238 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
8          10    99518968                    C                     T rs11190133    0.8673190   0.00450544  0.0269683    0.392278 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R


以下是proxy替换后的结果：
  chr.outcome pos.outcome other_allele.outcome effect_allele.outcome        SNP pval.outcome beta.outcome se.outcome eaf.outcome  outcome mr_keep.outcome pval_origin.outcome id.outcome
1          17    38143548                    C                     T  rs4065321    0.1598920   0.03709510  0.0263940    0.480414 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
2           6   109618704                    C                     T  rs9386796    0.2052730   0.03347070  0.0264242    0.505742 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
3           6    31143511                    G                     A   rs887466    0.0441448  -0.05326680  0.0264652    0.468238 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R
4          10   101280279                    A                     G rs17094148    0.8673190   0.00450544  0.0269683    0.392278 lung cancer (excluding cancer in controls)            TRUE            reported     1stH7R

step2复现成功。
