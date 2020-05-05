%% Results shown in Supplementary Table 4
clear
load('Result2.mat')
%% scenario I

% all enzymes
[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_pFBA)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

% enzymes with data on Vmax
[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_pFBA_temp)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_pFBA_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_pFBA_temp)),find(~isnan(kapp_Vmax_combined_all_max))))

%% scenario II
kmax_vivo_ratios=kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios};
kmax_vivo_ratios_temp=kmax_vivo_ratios_temp{best_row_ss_ratios,best_col_co_ratios};
co_ratio(best_col_co_ratios)    
ss_ratio(best_row_ss_ratios)

% all enzymes
[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_ratios)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

% enzymes with data on Vmax
[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_ratios_temp)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_temp),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_ratios_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_ratios_temp)),find(~isnan(kapp_Vmax_combined_all_max))))

%% scenario III
kmax_vivo_ratios_Vmax=kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax};
kmax_vivo_ratios_Vmax_temp=kmax_vivo_ratios_Vmax_temp{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax};
co_ratio(best_col_co_ratio_Vmax)    
ss_ratio(best_row_ss_ratios_Vmax)

% enzymes without data on Vmax
[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_ratios_Vmax_temp)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_ratios_Vmax_temp),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

%% scenario IV
kmax_vivo_light=kmax_vivo_light{best_row_ss_light,best_col_co_light};
kmax_vivo_light_temp=kmax_vivo_light_temp{best_row_ss_light,best_col_co_light};
co_ratio(best_col_co_light)    
ss_ratio(best_row_ss_light)

% all enzymes
[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_light)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

% enzymes with data on Vmax
[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_mean),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_mean),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_light_temp)),find(~isnan(kcat_BRENDA_mean))))

[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_max),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_median),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light_temp),log(kcat_BRENDA_median),'rows','pairwise','type','Spearman')

[C,P]=corr(log(kmax_vivo_light_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise')
[C,P]=corr(log(kmax_vivo_light_temp),log(kapp_Vmax_combined_all_max),'rows','pairwise','type','Spearman')
length(intersect(find(~isnan(kmax_vivo_light_temp)),find(~isnan(kapp_Vmax_combined_all_max))))


