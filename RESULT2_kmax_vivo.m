%% $k_{\mathrm{max}}^{\mathrm{vivo}}$ and _in vitro_ $k_{\mathrm{cat}}$

clear
%% 
% $k_{\mathrm{max}}^{\mathrm{vivo}}$=$\frac{v\left(C\right)}{E\left(C\right)}$

abundance_data
kcat_data_plant

abundance_per_EC_mg_gDW_combined_all = [abundance_per_EC_mg_gDW_Pyl abundance_per_EC_mg_gDW_Piques abundance_per_EC_mg_gDW_Seaton];
%% Flux estimation $v\left(C\right)$
%% (1) pFBA and measured growth rate

flux_estimation_pFBA
Diff_flux_pFBA = Diff_flux;

v_umol_gDW_sec(v_umol_gDW_sec==0) = nan;

kapp_pFBA = v_umol_gDW_sec./abundance_per_EC_mg_gDW_combined_all;

kmax_vivo_pFBA = max(kapp_pFBA,[],2,'omitnan');
clear v_umol_gDW_sec Diff_flux
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and _in vitro _$k_{\mathrm{cat}}$ - all enzymes

disp('Number of enzymes:')
length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(kmax_vivo_pFBA))))

disp('Number of enzymes with etha > 1:')
sum(kcat_BRENDA_max<kmax_vivo_pFBA)

disp('Correlation kmax_vivo to mean kcat:')
C=corr(log(kcat_BRENDA_mean),log(kmax_vivo_pFBA),'rows','pairwise')
Csp=corr(log(kcat_BRENDA_mean),log(kmax_vivo_pFBA),'type','Spearman','rows','pairwise')

disp('Correlation kmax_vivo to max kcat:')
C=corr(log(kcat_BRENDA_max),log(kmax_vivo_pFBA),'rows','pairwise')
Csp=corr(log(kcat_BRENDA_mean),log(kmax_vivo_pFBA),'type','Spearman','rows','pairwise')
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and _in vitro_ $k_{\mathrm{cat}}$ - enzymes with Vmax, kcat and v estimate

disp('Number of enzymes:')
load('Result1.mat',"kapp_Vmax_combined_all_max")
kmax_vivo_pFBA_temp=kmax_vivo_pFBA;
kmax_vivo_pFBA_temp(isnan(kapp_Vmax_combined_all_max)) = nan;
kmax_vivo_pFBA_temp(isnan(kcat_BRENDA_max)) = nan;

length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(kmax_vivo_pFBA_temp))))

disp('Correlation kmax_vivo to mean kcat:')
C=corr(log(kcat_BRENDA_mean),log(kmax_vivo_pFBA_temp),'rows','pairwise')
Csp=corr(log(kcat_BRENDA_mean),log(kmax_vivo_pFBA_temp),'type','Spearman','rows','pairwise')

disp('Correlation kmax_vivo to max kcat:')
C=corr(log(kcat_BRENDA_max),log(kmax_vivo_pFBA_temp),'rows','pairwise')
Csp=corr(log(kcat_BRENDA_max),log(kmax_vivo_pFBA_temp),'type','Spearman','rows','pairwise')
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}^{\mathrm{Vmax}}$ - 22 enzymes with Vmax

disp('Number of enzymes:')
length(intersect(find(~isnan(kapp_Vmax_combined_all_max)),find(~isnan(kmax_vivo_pFBA_temp))))

disp('Correlation kmax_vivo to kcat_Vmax:')
C=corr(log(kapp_Vmax_combined_all_max),log(kmax_vivo_pFBA_temp),'rows','pairwise')
Csp=corr(log(kapp_Vmax_combined_all_max),log(kmax_vivo_pFBA_temp),'type','Spearman','rows','pairwise')
%% (2) pFBA and constraint c/o ratio, starch/sucrose ratio, and growth

flux_estimation_co_ratios
Diff_flux_ratios = Diff_flux;

[numRow,numCol]=size(v_umol_gDW_sec)
for i=1:numRow
    for j=1:numCol
        v_umol_gDW_sec{i,j}(v_umol_gDW_sec{i,j}==0)=nan;
    end
end

kapp_ratios = cellfun(@(x) x./abundance_per_EC_mg_gDW_combined_all,v_umol_gDW_sec,'un',0);
kmax_vivo_ratios = cellfun(@(x) max(x,[],2,'omitnan'),kapp_ratios,'UniformOutput',false);
clear v_umol_gDW_sec Diff_flux
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}$ - all enzymes

disp('Number of enzymes:')
cellfun(@(x) length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(x)))), kmax_vivo_ratios,'UniformOutput',false)

disp('Number of enzymes with etha > 1:')
cellfun(@(x) sum(kcat_BRENDA_max<x), kmax_vivo_ratios)

disp('Correlation kmax_vivo to mean kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'rows','pairwise'), kmax_vivo_ratios,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")

disp('Correlation kmax_vivo to max kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'rows','pairwise'), kmax_vivo_ratios,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}$ - enzymes with Vmax, kcat and v estimate

disp('Number of enzymes:')
load('Result1.mat',"kapp_Vmax_combined_all_max")
kmax_vivo_ratios_temp=kmax_vivo_ratios;
for i=1:numRow
    for j=1:numCol
        kmax_vivo_ratios_temp{i,j}(isnan(kapp_Vmax_combined_all_max))=nan;
        kmax_vivo_ratios_temp{i,j}(isnan(kcat_BRENDA_max))=nan;
    end
end

cellfun(@(x) length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(x)))), kmax_vivo_ratios_temp,'UniformOutput',false)

disp('Correlation kmax_vivo to mean kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")

disp('Correlation kmax_vivo to max kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}^{\mathrm{Vmax}}$ - enzymes with Vmax, kcat and v estimate

disp('Number of enzymes:')
cellfun(@(x) length(intersect(find(~isnan(kapp_Vmax_combined_all_max)),find(~isnan(x)))), kmax_vivo_ratios_temp,'UniformOutput',false)

disp('Correlation kmax_vivo to kcat_Vmax:')
C=cell2mat(cellfun(@(x) corr(log(kapp_Vmax_combined_all_max),log(x),'rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kapp_Vmax_combined_all_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios_temp,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")
[best_row_ss_ratios,best_col_co_ratios]=find(C==max(C,[],"all"))
% 
%% (3) pFBA and constraint c/o ratio, starch/sucrose ratio, and growth, Vmax as upper bounds

flux_estimation_co_ratios_Vmax
Diff_flux_ratios_Vmax = Diff_flux;

[numRow,numCol]=size(v_umol_gDW_sec);
for i=1:numRow
    for j=1:numCol
        v_umol_gDW_sec{i,j}(v_umol_gDW_sec{i,j}==0)=nan;
    end
end

kapp_ratios_Vmax = cellfun(@(x) x./abundance_per_EC_mg_gDW_combined_all,v_umol_gDW_sec,'un',0);
kmax_vivo_ratios_Vmax = cellfun(@(x) max(x,[],2,'omitnan'),kapp_ratios_Vmax,'UniformOutput',false);
clear v_umol_gDW_sec Diff_flux
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}$ -  correlation over enzymes without known Vmax

load('Result1.mat',"kapp_Vmax_combined_all_max")
kmax_vivo_ratios_Vmax_temp=kmax_vivo_ratios_Vmax;
for i=1:numRow
    for j=1:numCol
        kmax_vivo_ratios_Vmax_temp{i,j}(~isnan(kapp_Vmax_combined_all_max))=nan;
    end
end

disp('Number of enzymes:')
cellfun(@(x) length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(x)))), kmax_vivo_ratios_Vmax_temp,'UniformOutput',false)

disp('Number of enzymes with etha > 1:')
cellfun(@(x) sum(kcat_BRENDA_max<x), kmax_vivo_ratios_Vmax_temp)

disp('Correlation kmax_vivo to mean kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'rows','pairwise'), kmax_vivo_ratios_Vmax_temp,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios_Vmax_temp,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")

disp('Correlation kmax_vivo to max kcat:')
C=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'rows','pairwise'), kmax_vivo_ratios_Vmax_temp,'UniformOutput',false))
Csp=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_ratios_Vmax_temp,'UniformOutput',false))
max(C,[],"all")
max(Csp,[],"all")
[best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax]=find(C==max(C,[],"all"))
%% (4) pFBA and constraint c/o ratio, starch/sucrose ratio, and growth, light and CO2 uptake

flux_estimation_light_CO2
Diff_flux_light = Diff_flux;

[numRow,numCol]=size(v_umol_gDW_sec);
for i=1:numRow
    for j=1:numCol
        v_umol_gDW_sec{i,j}(v_umol_gDW_sec{i,j}==0)=nan;
    end
end

kapp_light = cellfun(@(x) x./abundance_per_EC_mg_gDW_combined_all,v_umol_gDW_sec,'un',0);
kmax_vivo_light = cellfun(@(x) max(x,[],2,'omitnan'),kapp_light,'UniformOutput',false);
clear v_umol_gDW_sec Diff_flux
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}$ - all enzymes

disp('Number of enzymes:')
X=cellfun(@(x) length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(x)))), kmax_vivo_light,'UniformOutput',false)

disp('Number of enzymes with etha > 1:')
cellfun(@(x) sum(kcat_BRENDA_max<x), kmax_vivo_light)

disp('Correlation kmax_vivo to mean kcat:')
C_mean=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'rows','pairwise'), kmax_vivo_light,'UniformOutput',false))
Csp_mean=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_light,'UniformOutput',false))
max(C_mean,[],'all')
max(Csp_mean,[],'all')

disp('Correlation kmax_vivo to max kcat:')
C_max=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'rows','pairwise'), kmax_vivo_light,'UniformOutput',false))
Csp_max=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_light,'UniformOutput',false))
max(C_max,[],'all')
max(Csp_max,[],'all')
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}$ - enzymes with Vmax, kcat and v estimate

disp('Number of enzymes:')
load('Result1.mat',"kapp_Vmax_combined_all_max")
kmax_vivo_light_temp=kmax_vivo_light;
for i=1:numRow
    for j=1:numCol
        kmax_vivo_light_temp{i,j}(isnan(kapp_Vmax_combined_all_max))=nan;
        kmax_vivo_light_temp{i,j}(isnan(kcat_BRENDA_max))=nan;
    end
end

X=cellfun(@(x) length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(x)))), kmax_vivo_light_temp,'UniformOutput',false)

disp('Correlation kmax_vivo to mean kcat:')
C_all_data_mean=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
Csp_all_data_mean=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_mean),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
max(C_all_data_mean,[],'all')
max(Csp_all_data_mean,[],'all')

disp('Correlation kmax_vivo to max kcat:')
C_all_data_max=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
Csp_all_data_max=cell2mat(cellfun(@(x) corr(log(kcat_BRENDA_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
max(C_all_data_max,[],'all')
max(Csp_all_data_max,[],'all')
% Correlation $k_{\mathrm{max}}^{\mathrm{vivo}}$ and $k_{\mathrm{cat}}^{\mathrm{Vmax}}$ - enzymes with Vmax, kcat and v estimate

disp('Number of enzymes:')
cellfun(@(x) length(intersect(find(~isnan(kapp_Vmax_combined_all_max)),find(~isnan(x)))), kmax_vivo_light_temp,'UniformOutput',false)

disp('Correlation kmax_vivo to kcat_Vmax:')
C_Vmax=cell2mat(cellfun(@(x) corr(log(kapp_Vmax_combined_all_max),log(x),'rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
Csp_Vmax=cell2mat(cellfun(@(x) corr(log(kapp_Vmax_combined_all_max),log(x),'type','Spearman','rows','pairwise'), kmax_vivo_light_temp,'UniformOutput',false))
max(C_Vmax,[],"all")
max(Csp_Vmax,[],"all")
[best_row_ss_light,best_col_co_light]=find(C_Vmax==max(C_Vmax,[],"all"))

length(intersect(find(~isnan(kapp_Vmax_combined_all_max)),find(~isnan(kmax_vivo_light_temp{best_row_ss_light,best_col_co_light}))))

save('Result2.mat')
%% How similar are the estimated $k_{\mathrm{max}}^{\mathrm{vivo}}$ ?

[C,P]=corr([kmax_vivo_pFBA kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios} ...
    kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax} ...
    kmax_vivo_light{best_row_ss_light,best_col_co_light}],'rows','pairwise');
C_upper=triu(C);
C_upper(C_upper==0) = nan;
P_upper=triu(P);
P_upper(C_upper==0) = nan;

figure 
subplot(1,2,1)
heatmap(C_upper,{'growth', 'flux ratios', 'Vmax', 'uptake rates'},{'pFBA', 'ratios', 'Vmax', 'light'},'%0.2f','Colormap',jet,'ColorBar',0,'MinColorValue',-1,'MaxColorValue',1,'NaNColor',[1 1 1],'GridLines', ':')

[C,P]=corr([kmax_vivo_pFBA kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios} ...
    kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax} ...
    kmax_vivo_light{best_row_ss_light,best_col_co_light}],'type','Spearman','rows','pairwise');
C_upper=triu(C);
C_upper(C_upper==0) = nan;
P_upper=triu(P);
P_upper(C_upper==0) = nan;

subplot(1,2,2)
heatmap(C_upper,{'growth', 'flux ratios', 'Vmax', 'uptake rates'},{'pFBA', 'ratios', 'Vmax', 'light'},'%0.2f','Colormap',jet,'ColorBar',0,'MinColorValue',-1,'MaxColorValue',1,'NaNColor',[1 1 1],'GridLines', ':')

id=find(all(~isnan([kmax_vivo_pFBA kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios} ...
    kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax} ...
    kmax_vivo_light{best_row_ss_light,best_col_co_light}]')));

CV_kmax_vivo_v_estimates=std([kmax_vivo_pFBA kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios} ...
    kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax} ...
    kmax_vivo_light{best_row_ss_light,best_col_co_light}]')./mean([kmax_vivo_pFBA kmax_vivo_ratios{best_row_ss_ratios,best_col_co_ratios} ...
    kmax_vivo_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax} ...
    kmax_vivo_light{best_row_ss_light,best_col_co_light}]')

mean(CV_kmax_vivo_v_estimates,'omitnan')
median(CV_kmax_vivo_v_estimates,'omitnan')
max(CV_kmax_vivo_v_estimates,[],'omitnan')
sum(~isnan(CV_kmax_vivo_v_estimates))
%% Figures
%% Figure 4

% Figure 4
load("Names_of_example_enzymes.mat")
figure
subplot(1,3,1)
plot(kapp_Vmax_combined_all_max,kmax_vivo_light{best_row_ss_light,best_col_co_light},'.');
text(kapp_Vmax_combined_all_max,kmax_vivo_light{best_row_ss_light,best_col_co_light},NAMES,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8)

xlabel('k_{cat}^{Vmax} [s^{-1}]')
ylabel('k_{max}^{vivo} [s^{-1}]')
set(gca,'XScale','Log','YScale','Log')
hold on
plot([1e-2 1e6],[1e-2 1e6],'k')
hold on
[~,m,b]=regression(log(kapp_Vmax_combined_all_max)',log(kmax_vivo_light{best_row_ss_light,best_col_co_light})');
x=log([1e-4; 1e-3; 1e-2; 1e-1;1e0; 1e1; 1e2;1e3;1e4;1e5;1e6; 1e7]);
y=(m.*x)+b;
plot(exp(x),exp(y),'-b')
xlim([1e-2 1e6])
ylim([1e-2 1e6])

subplot(1,3,2)
plot(kcat_BRENDA_mean,kmax_vivo_light{best_row_ss_light,best_col_co_light},'.');
text(kcat_BRENDA_mean,kmax_vivo_light{best_row_ss_light,best_col_co_light},NAMES,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8)

xlabel('in vitro k_{cat} [s^{-1}]')
% ylabel('k_{app}')
set(gca,'XScale','Log','YScale','Log')
hold on
xlim([1e-2 1e6])
ylim([1e-2 1e6])
plot([1e-4 1e7],[1e-4 1e7],'k')
hold on
[~,m,b]=regression(log(kcat_BRENDA_mean)',log(kmax_vivo_light{best_row_ss_light,best_col_co_light})');
x=log([1e-4; 1e-3; 1e-2; 1e-1;1e0; 1e1; 1e2;1e3;1e4;1e5;1e6; 1e7]);
y=(m.*x)+b;
plot(exp(x),exp(y),'-b')

x=log(kmax_vivo_light{best_row_ss_light,best_col_co_light}(kcat_BRENDA_mean>0));
y=log(kcat_BRENDA_mean(kcat_BRENDA_mean>0));
NAMES_1=NAMES((kcat_BRENDA_mean>0));

subplot(1,3,3)
plot(exp(y(Comp.number(kcat_BRENDA_mean>0)==1)),exp(x(Comp.number(kcat_BRENDA_mean>0)==1)),'.');
text(exp(y(Comp.number(kcat_BRENDA_mean>0)==1)),exp(x(Comp.number(kcat_BRENDA_mean>0)==1)),NAMES_1(Comp.number(kcat_BRENDA_mean>0)==1),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',8)

xlabel('in vitro k_{cat} [s^{-1}]')
% ylabel('k_{app}')
set(gca,'XScale','Log','YScale','Log')
hold on
xlim([1e-2 1e6])
ylim([1e-2 1e6])
plot([1e-4 1e7],[1e-4 1e7],'k')
hold on
[~,m,b]=regression(y(Comp.number(kcat_BRENDA_mean>0)==1)',x(Comp.number(kcat_BRENDA_mean>0)==1)');
xr=log([1e-4; 1e-3; 1e-2; 1e-1;1e0; 1e1; 1e2;1e3;1e4;1e5;1e6; 1e7]);
yr=(m.*xr)+b;
plot(exp(xr),exp(yr),'-b')

disp('Correlation for enzymes in single compartment:')
[C,P]=corr(x(Comp.number(kcat_BRENDA_mean>0)==1),y(Comp.number(kcat_BRENDA_mean>0)==1),'rows','pairwise')
disp('Number of enzymes:')
length(intersect(find(~isnan(x(Comp.number(kcat_BRENDA_mean>0)==1))),find(~isnan(y(Comp.number(kcat_BRENDA_mean>0)==1)))))
%% 
% %% Supplementary Figure 3

disp('Flux variability in flux estimates:')
figure
subplot(2,2,1)
x=max((Diff_flux_pFBA./max(model_irr.ub)),[],2,'omitnan');
[H,b]=hist(x(:),50);
bar(b,(H/sum(H)),1)
ylabel('Fraction')
title('A')

subplot(2,2,2)
x=max((Diff_flux_ratios{best_row_ss_ratios,best_col_co_ratios}./max(model_irr.ub)),[],2,'omitnan');
[H,b]=hist(x(:),50);
bar(b,(H/sum(H)),1)
title('B')

subplot(2,2,3)
x=max((Diff_flux_ratios_Vmax{best_row_ss_ratios_Vmax,best_col_co_ratio_Vmax}./max(model_irr.ub)),[],2,'omitnan');
[H,b]=hist(x(:),50);
bar(b,(H/sum(H)),1)
xlabel('Flux variability with respect to maximum flux [%]')
ylabel('Fraction')
title('C')

subplot(2,2,4)
x=max((Diff_flux_light{best_row_ss_light,best_col_co_light}./max(model_irr.ub)),[],2,'omitnan');
[H,b]=hist(x(:),50);
bar(b,(H/sum(H)),1)
%xlabel('Flux variability with respect to maximum flux [%]')
title('D')
%% Supplementary Figure 4

% SUPPLEMENTARY FIGURE 4 
% k_max_vivo ~ in vitro kcat (all enzymes)
figure
subplot(2,2,1)
heatmap(C_Vmax,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',0,'MinColorValue',-0.1,'MaxColorValue',0.7,'NaNColor',[1 1 1],'GridLines', ':')
% xlabel('carboxylation/oxygenation ratio')
ylabel('starch/sucrose synthases')
title('k_{max}^{vivo} and k_{cat}^{Vmax}')

subplot(2,2,2)
heatmap(C_all_data_mean,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',1,'MinColorValue',-0.1,'MaxColorValue',0.7,'NaNColor',[1 1 1],'GridLines', ':')
% xlabel('carboxylation/oxygenation ratio')
% ylabel('starch/sucrose synthase ratio')
title('k_{max}^{vivo} and mean k_{cat} for enzymes with available V_{max}')

subplot(2,2,3)
heatmap(C_mean,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',0,'MinColorValue',-0.1,'MaxColorValue',0.7,'NaNColor',[1 1 1],'GridLines', ':')
xlabel('carboxylation/oxygenation')
ylabel('starch/sucrose synthases')
title('k_{max}^{vivo} and mean k_{cat}')

subplot(2,2,4)
heatmap(C_max,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',1,'MinColorValue',-0.1,'MaxColorValue',0.7,'NaNColor',[1 1 1],'GridLines', ':')
xlabel('carboxylation/oxygenation')
% ylabel('starch/sucrose synthase ratio')
title('k_{max}^{vivo} and maximum k_{cat}')