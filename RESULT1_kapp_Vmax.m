%% $k_{\mathrm{app}}^{\mathrm{Vmax}}$ estimation

clear
%% 
% $k_{\mathrm{app}}^{\mathrm{Vmax}}$=$\frac{V_{\mathrm{max}} \left(C\right)}{E\left(C\right)}$

abundance_data
activity_data
%% Sulpice et al. data set

kapp_Vmax_Seaton = (vmax_Seaton./abundance_per_EC_mg_gFW_Seaton)/60;
%% Piques et al. data set

kapp_Vmax_Piques = (vmax_Piques./abundance_per_EC_mg_gFW_Piques)/60;
%% Pyl et al. data set

kapp_Vmax_Pyl = (vmax_Pyl./abundance_per_EC_mg_gFW_Pyl)/60;
%% 
% 
% Combine data set to one (Piques Pyl Seaton data)
vmax_combined_all = [vmax_Pyl vmax_Piques vmax_Seaton];
kapp_Vmax_combined_all = [kapp_Vmax_Pyl kapp_Vmax_Piques kapp_Vmax_Seaton];
abundance_per_EC_mg_gFW_combined_all = [abundance_per_EC_mg_gFW_Pyl abundance_per_EC_mg_gFW_Piques abundance_per_EC_mg_gFW_Seaton];
abundance_per_EC_mg_gDW_combined_all = [abundance_per_EC_mg_gDW_Pyl abundance_per_EC_mg_gDW_Piques abundance_per_EC_mg_gDW_Seaton];
%% RESULT 1
%% VARIANCE $k_{\mathrm{app}}^{\mathrm{Vmax}}$

C=triu(corr(kapp_Vmax_combined_all,'rows','pairwise'),1);

disp('Minimum correlation kapp_Vmax across conditions:')
min(C(C>0),[],'all')

disp('Average correlation kapp_Vmax across conditions:')
mean(C(C>0),'all')

save('Result1.mat')
%% Supplementary Figure 2

figure
label={'Pyl ED 12°C/12°C';'Pyl ED 16°C/16°C';'Pyl ED 24°C/12°C';'Pyl ED 24°C/16°C';'Pyl ED 24°C/24°C';'Piques ED';'Seaton ED 6-h';'Seaton ED 8-h';'Seaton ED 12-h';'Seaton ED 18-h'};
heatmap(corr(kapp_Vmax_combined_all,'rows','pairwise'),...
    label,label,[],'Colormap',jet,'ColorBar',0,'MinColorValue',0,'MaxColorValue',1,'NaNColor',[1 1 1],'GridLines', ':','TickAngle',45)

var_id = find(sum(kapp_Vmax_combined_all'>0));
kcat_estimated_vmax_combined = kapp_Vmax_combined_all(var_id,:);
vmax_combined = vmax_combined_all(var_id,:);
abundance_combined = abundance_per_EC_mg_gFW_combined_all(var_id,:);

disp('----------------------------------------------')
disp('CV of estimated kapp_Vmax across conditions:')
CV=std(kapp_Vmax_combined_all(~all(isnan(kapp_Vmax_combined_all')),:)',0,1,"omitnan")./mean(kapp_Vmax_combined_all(~all(isnan(kapp_Vmax_combined_all')),:)','omitnan')
disp('average CV of estimated kapp_Vmax across conditions:')
mean(CV)

%% Supplementary Figure 1

figure
boxplot(kapp_Vmax_combined_all(~all(isnan(kapp_Vmax_combined_all')),:)')
ylabel('k_{app}^{Vmax} [s^{-1}]')
xlabel('EC number')
xlim([0 length(var_id)+1])
set(gca,'XTick',1:length(var_id),'XTickLabel',EC(var_id),'XTickLabelRotation',45,'YScale','Log')