result <- split(x, x$chr)          #由于运行clump_data函数遇到timeout 300 seconds issue, 解决办法按染色体分割
for (i in 1:22){                   #先查看这个文件染色体列unique值，发现没有X，Y。因此循环从1到22
  write.table(result[[i]], file = paste0("abc_chr",i,sep=""), quote = FALSE, row.names = FALSE) #按染色体出文件，备用
  newdata <-subset(result[[i]], pval.exposure < 5e-08)                                          #直接clump_data()染色体文件报错 trafficing的问题，因此先筛选 pvalue<5e-08,再clump
  y <- clump_data(newdata, clump_p1 = 5e-08, clump_p2 = 5e-08)                                  #clump
  write.table(y, file = paste0(file,"_chr",i,sep=""), quote = FALSE, row.names = FALSE)         #出文件，备用
  #x <- rbind(assign(paste0(file,"_chr",i,sep=""),subset(result[[i]], pval.exposure < 5e-08)))
  }
  z <- rbind(chr6,chr17,chr10) 

一、最后的结果是4条 GrimAge，与文献所述相同，结果如下：

          SNP chr.exposure pos.exposure effect_allele.exposure other_allele.exposure eaf.exposure beta.exposure se.exposure   pval.exposure samplesize.exposure exposure mr_keep.exposure pval_origin.exposure id.exposure           r2            F
1   rs4065321           17     38143548                      T                     C       0.5333       -0.1703      0.0296 8.776e-09               32420 DNA_methylation_ageing             TRUE             reported     GrimAge 0.0010199753 33.09932
11  rs9386796            6    109618704                      T                     C       0.4597        0.1983      0.0294 1.640e-11               32418 DNA_methylation_ageing             TRUE             reported     GrimAge 0.0014013787 45.49084
23   rs887466            6     31143511                      A                     G       0.3815       -0.1928      0.0310 5.092e-10               30674 DNA_methylation_ageing             TRUE             reported     GrimAge 0.0012594269 38.67785
12 rs17094148           10    101280279                      A                     G       0.7067       -0.1800      0.0323 2.547e-08               31838 DNA_methylation_ageing             TRUE             reported     GrimAge 0.0009744752 31.05365



二、进行r2和F检验结果，均与文献所述相同。
variance_GrimAge = 0.004655256， #文献为0.47%
Fmin_GrimAge = 31.05365，        #文献为31
Fmax_GrimAge = 45.49084          #文献为45
Fmedian_GrimAge= 35.88859        #文献为36

Step1复现成功。

结果意义：在GrimAge GWAS中筛选pvalue<5*10-8,连锁不平衡的位点，用于TwoSampleMR分析。r2和F检验验证所得位点的显著性和精确性。
