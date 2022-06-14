结果如下：
id.exposure id.outcome                                    outcome               exposure                    method nsnp           b         se      pval
1     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing Inverse variance weighted    4  0.08087056 0.10288551 0.4318530
2     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing                  MR Egger    4  2.84249600 1.20909319 0.1430953
3     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing           Weighted median    4  0.11590759 0.08704695 0.1830073
4     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing             Weighted mode    4  0.16921787 0.12305455 0.2627859

结果意义：
MR分析输入exposure和outcome文件，Exposure文件为遗传变异与暴露因素的gwas结果，Outcome文件为遗传变异与结局变量的gwas结果。这两个文件我们分别在step1和2format好。
Harmonise data主要目的用于将SNP位点统一调整成正链，并且根据allele和频率判断两个gwas结果中的SNP位点是否一致，不一致的进行去除。
通过4种算法计算出MR结果不同，我们主要看你方差加权法IVW的结果，其他作为稳定性和可靠性的参考。b相当于OR值，是正向或者反向量化的一个指标。
IVW MR法是假设在所有遗传变异都是有效工具变量的基本前提下进行的有效分析，具有较强的因果关系检测额能力。
