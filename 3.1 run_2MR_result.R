data_FINNGEN <- harmonise_data(
   exposure_dat = exp_dat, 
   outcome_dat = outcome
 )
Harmonising DNA_methylation_ageing (GrimAge) and lung cancer (excluding cancer in controls) (zIV0uA)
Removing the following SNPs for incompatible alleles:
rs4065321
> results_FINNGEN <- mr(data_FINNGEN, method_list = c("mr_ivw", "mr_egger_regression", 
                                         "mr_weighted_median", "mr_weighted_mode"))
Analysing 'GrimAge' on 'zIV0uA'

> results_FINNGEN
  id.exposure id.outcome                                    outcome               exposure                    method nsnp
1     GrimAge     zIV0uA lung cancer (excluding cancer in controls) DNA_methylation_ageing Inverse variance weighted    3
2     GrimAge     zIV0uA lung cancer (excluding cancer in controls) DNA_methylation_ageing                  MR Egger    3
3     GrimAge     zIV0uA lung cancer (excluding cancer in controls) DNA_methylation_ageing           Weighted median    3
4     GrimAge     zIV0uA lung cancer (excluding cancer in controls) DNA_methylation_ageing             Weighted mode    3
          b         se       pval
1 0.1615007 0.08058135 0.04504917
2 1.9807200 2.01513549 0.50548294
3 0.1723837 0.09897681 0.08156846
4 0.1952318 0.11767002 0.23895296
