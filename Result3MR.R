结果意义：
MR分析输入exposure和outcome文件，Exposure文件为遗传变异与暴露因素的gwas结果，Outcome文件为遗传变异与结局变量的gwas结果。这两个文件我们分别在step1和2format好。
Harmonise data主要目的用于将SNP位点统一调整成正链，并且根据allele和频率判断两个gwas结果中的SNP位点是否一致，不一致的进行去除。
通过4种算法计算出MR结果不同，我们主要看逆方差加权法IVW的结果，其他作为稳定性和可靠性的参考。b相当于OR值，是正向或者反向量化的一个指标。
IVW MR法是假设在所有遗传变异都是有效工具变量的基本前提下进行的有效分析，具有较强的因果关系检测额能力。

Exposure=GrimAge, Outcome=Finngen- Lung cancer & Colon cancer结果如下：
id.exposure id.outcome                                    outcome               exposure                    method nsnp           b         se      pval
1     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing Inverse variance weighted    4  0.08087056 0.10288551 0.4318530
2     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing                  MR Egger    4  2.84249600 1.20909319 0.1430953
3     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing           Weighted median    4  0.11590759 0.08704695 0.1830073
4     GrimAge     1stH7R lung cancer (excluding cancer in controls) DNA_methylation_ageing             Weighted mode    4  0.16921787 0.12305455 0.2627859

1     GrimAge     VTERLA colorectal cancer (excluding cancer in controls) DNA_methylation_ageing Inverse variance weighted    4 0.06011094 0.07159065 0.4011058
2     GrimAge     VTERLA colorectal cancer (excluding cancer in controls) DNA_methylation_ageing                  MR Egger    4 1.46685020 1.21036524 0.3492941
3     GrimAge     VTERLA colorectal cancer (excluding cancer in controls) DNA_methylation_ageing           Weighted median    4 0.09028742 0.08262885 0.2745315
4     GrimAge     VTERLA colorectal cancer (excluding cancer in controls) DNA_methylation_ageing             Weighted mode    4 0.11307694 0.11798135 0.4085507

Exposure=GrimAge, Outcome=MRBase- Lung cancer结果如下：
id.exposure id.outcome                                    outcome               exposure                    method nsnp           b         se      pval
"GrimAge" "ieu-a-965" "Lung adenocarcinoma || id:ieu-a-965" "DNA_methylation_ageing" "Inverse variance weighted" 3 0.0325053968745173 0.0891784022870953 0.715485768514375
"GrimAge" "ieu-a-965" "Lung adenocarcinoma || id:ieu-a-965" "DNA_methylation_ageing" "MR Egger" 3 -0.365195082996395 1.37408918199732 0.83462694970153
"GrimAge" "ieu-a-965" "Lung adenocarcinoma || id:ieu-a-965" "DNA_methylation_ageing" "Weighted median" 3 0.0373098350486472 0.104599227712093 0.721321483870316
"GrimAge" "ieu-a-965" "Lung adenocarcinoma || id:ieu-a-965" "DNA_methylation_ageing" "Weighted mode" 3 0.0517349504203492 0.123990551372079 0.717019829828505
"GrimAge" "ieu-a-966" "Lung cancer || id:ieu-a-966" "DNA_methylation_ageing" "Inverse variance weighted" 3 -0.0385126915922153 0.0572361481040768 0.501027724858474
"GrimAge" "ieu-a-966" "Lung cancer || id:ieu-a-966" "DNA_methylation_ageing" "MR Egger" 3 0.829277415004235 0.884466583404515 0.520494529837516
"GrimAge" "ieu-a-966" "Lung cancer || id:ieu-a-966" "DNA_methylation_ageing" "Weighted median" 3 -0.0366078923010231 0.0689839804889188 0.595645837102256
"GrimAge" "ieu-a-966" "Lung cancer || id:ieu-a-966" "DNA_methylation_ageing" "Weighted mode" 3 -0.0175680318253957 0.0839428862009827 0.853607115394654
"GrimAge" "ieu-a-967" "Squamous cell lung cancer || id:ieu-a-967" "DNA_methylation_ageing" "Inverse variance weighted" 3 -0.0724792257783886 0.154905194709269 0.639860326997747
"GrimAge" "ieu-a-967" "Squamous cell lung cancer || id:ieu-a-967" "DNA_methylation_ageing" "MR Egger" 3 3.11027302082866 1.37152067786909 0.264397634248977
"GrimAge" "ieu-a-967" "Squamous cell lung cancer || id:ieu-a-967" "DNA_methylation_ageing" "Weighted median" 3 -0.189818411460638 0.11462730810857 0.0977296115441971
"GrimAge" "ieu-a-967" "Squamous cell lung cancer || id:ieu-a-967" "DNA_methylation_ageing" "Weighted mode" 3 -0.241428428018679 0.182168817088319 0.316202968011809
