Meta-analysis 从相关研究中整合定量数据，并产生能够总结研究整体的结果。
fixed effect假定所有研究都估计出相同的治疗效果，只是因为取样的偶然（chance）差异导致观察到的treatment effect不同。
合并之前，首先应该进行异质性检验，保证现有的各独立研究间的结果的不同仅仅是由于抽样误差造成的。

Number of studies combined: k = 6  #本次纳入meta分析的研究（或效应量）有6个（lung cancer两个），因此k=6。


                        OR           95%-CI     z p-value
Common effect model 0.9613 [0.9357; 0.9875] -2.88  0.0040      #综合效应量的结果
#Common effect model计算的综合效应量为0.9613，95%CI为[0.9357; 0.9875]，z检验值为-2.88，p<0.0040。【结果显著】

#以下是异质性检验的结果。
Quantifying heterogeneity:
 tau^2 = 0 [0.0000; 0.0136]; tau = 0 [0.0000; 0.1166]
 I^2 = 0.0% [0.0%; 74.6%]; H = 1.00 [1.00; 1.99]

Test of heterogeneity:
    Q d.f. p-value
 4.07    5  0.5390
#主要看Q检验（检验水准为α=0.1）的结果，以及I^2值（0～100%），一般来讲，I^2值0%-25%为低异质性, 50%-75%为中等异质性，>75%为高异质性。
#异质性越高，代表纳入的研究（或效应量）间差异越大。

Results for subgroups (common effect model):
                               k     OR           95%-CI    Q  I^2  tau^2    tau
subgroup = Ovarian cancer      1 0.9651 [0.8985; 1.0367] 0.00   --     --     --
subgroup = Breast cancer       1 0.9573 [0.9256; 0.9901] 0.00   --     --     --
subgroup = Lung cancer         2 0.9898 [0.8973; 1.0917] 1.03 2.7% 0.0002 0.0140
subgroup = Prostate cancer     1 0.9284 [0.8555; 1.0075] 0.00   --     --     --
subgroup = Colorectal cancer   1 1.0620 [0.9229; 1.2219] 0.00   --     --     --

Test for subgroup differences (common effect model):
                  Q d.f. p-value
Between groups 3.04    4  0.5504
Within groups  1.03    1  0.3106

Details on meta-analytical method:
- Inverse variance method
- DerSimonian-Laird estimator for tau^2
- Jackson method for confidence interval of tau^2 and tau
#五个亚组各自的综合效应量，可以看到：固定效应模型的亚组分析没有发现组间（between groups）显著的差异，表明数据采集的地区并不调节综合效应的结果。


IVW_GrimAge矫正后结果：
            b         se       pval           outcome exposure                    method      p.fdr      p.bon     p.hoch p.sig
1 -0.03548549 0.03648373 0.33073289    Ovarian cancer  GrimAge Inverse variance weighted 0.50138228 1.00000000 0.80221165      
2 -0.04368244 0.01718616 0.01103078     Breast cancer  GrimAge Inverse variance weighted 0.05515389 0.05515389 0.05515389     *
3 -0.01029789 0.05001739 0.83687975       Lung cancer  GrimAge Inverse variance weighted 0.83687975 1.00000000 0.83687975      
4 -0.07429318 0.04169544 0.07478090   Prostate cancer  GrimAge Inverse variance weighted 0.18695225 0.37390450 0.29912360      
5  0.06011094 0.07159065 0.40110583 Colorectal cancer  GrimAge Inverse variance weighted 0.50138228 1.00000000 0.80221165      
