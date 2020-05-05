clear
load_model

kcat=readtable('Data/kcat_final_no_mutant.csv','Delimiter','\t'); % kcat BRENDA
for i=1:length(EC)
    if isempty(find(strcmp(kcat.EC,EC(i))))
        if ~isempty(strfind(EC{i},'/'))
            EC_temp=strsplit(EC{i},'/');
            kcat_temp=[]; t_temp=[];
            for j=1:length(EC_temp)
                kcat_temp=[kcat_temp; kcat.kcat_1_sec(find(strcmp(kcat.EC,EC_temp(j))))];
                t_temp=[t_temp; kcat.tax_dist(find(strcmp(kcat.EC,EC_temp(j))))];
            end
            if ~isempty(kcat_temp)
                kcat_BRENDA(i,1) = mean(kcat_temp)*60;
                kcat_BRENDA_taxa(i,1) = max(t_temp);
            else
                kcat_BRENDA(i,1) = nan;
            end
        elseif ~isempty(strfind(EC{i},'+'))
            EC_temp=strsplit(EC{i},' + ');
            kcat_temp=[];t_temp=[];
            for j=1:length(EC_temp)
                kcat_temp=[kcat_temp; kcat.kcat_1_sec(find(strcmp(kcat.EC,EC_temp(j))))];
                t_temp=[t_temp; kcat.tax_dist(find(strcmp(kcat.EC,EC_temp(j))))];
            end
            if ~isempty(kcat_temp)
                kcat_BRENDA(i,1) = mean(kcat_temp)*60;
                kcat_BRENDA_taxa(i,1) = max(t_temp);
            else
                kcat_BRENDA(i,1) = nan;
            end
        else
            kcat_BRENDA(i,1) = nan;
        end
    else
        kcat_BRENDA(i,1) = mean(kcat.kcat_1_sec(find(strcmp(kcat.EC,EC(i)))));
        kcat_BRENDA_taxa(i,1) = min(kcat.tax_dist(find(strcmp(kcat.EC,EC(i)))));
    end
end

kcat_p=kcat_BRENDA;
kcat_p(kcat_BRENDA_taxa>18)=nan;

kcat_np=kcat_BRENDA;
kcat_np(kcat_BRENDA_taxa<=18)=nan;

load('Result1.mat','vmax_combined_all')
load('Result2.mat','best_row_ss_light','best_col_co_light','v_umol_gDW_day','abundance_per_EC_mg_gDW_Pyl','abundance_per_EC_mg_gFW_Piques','abundance_per_EC_mg_gFW_Seaton','Z','Comp','MW_*','kmax_vivo_light')

abundance_combined_all=[abundance_per_EC_mg_gDW_Pyl,abundance_per_EC_mg_gFW_Piques,abundance_per_EC_mg_gFW_Seaton];
abundance_combined_all=max(abundance_combined_all,[],2,'omitnan');
abundance_combined_all(abundance_combined_all==0)=nan;

vmax_combined_all=max(vmax_combined_all,[],2,'omitnan');
vmax_combined_all(vmax_combined_all==0)=nan;

Flux=v_umol_gDW_day{best_row_ss_light,best_col_co_light}/86400;
Flux=max(Flux,[],2,'omitnan');
Flux(Flux==0)=nan;

%% 
% 
%%%%%%% RESULTS %%%%%%%%%%

% v - flux
v=find(~isnan(Flux));
nan_v=find(isnan(Flux));
% a - abundance
abundance_combined_all(abundance_combined_all==0)=nan;
a=find(~isnan(abundance_combined_all));
nan_a=find(isnan(abundance_combined_all));
% kp - kcat plant
kp=find(~isnan(kcat_p));
nan_kp=find(isnan(kcat_p));
% kn - kcat non-plant
kn=find(~isnan(kcat_np));
nan_kn=find(isnan(kcat_np));
% M - Vmax
M=find(~isnan(vmax_combined_all));
nan_M=find(isnan(vmax_combined_all));

disp('Data presented in Figure 1 main text')

% all kn
all_kp=length(intersect(intersect(intersect(intersect(v,a),kp),M),nan_kn)) %15
% all knp
all_kn=length(intersect(intersect(intersect(intersect(v,a),kn),M),nan_kp)) %26

% v,a,kp
vakp=length(intersect(intersect(intersect(intersect(v,a),kp),nan_M),nan_kn)) %14
% v,a,kn
vakn=length(intersect(intersect(intersect(intersect(v,a),kn),nan_M),nan_kp)) %25

%vaM
vaM=length(intersect(intersect(intersect(intersect(v,a),nan_kp),M),nan_kn)) %13

% v,kp,M
vkpM=length(intersect(intersect(intersect(intersect(v,nan_a),kp),M),nan_kn)) %12
% v,kn,M
vknM=length(intersect(intersect(intersect(intersect(v,nan_a),kn),M),nan_kp)) %23

% a,kp,M
akpM=length(intersect(intersect(intersect(intersect(nan_v,a),kp),M),nan_kn)) %11
% a,kn,M
aknM=length(intersect(intersect(intersect(intersect(nan_v,a),kn),M),nan_kp)) %22

% v,a
va=length(intersect(intersect(intersect(intersect(v,a),nan_kp),nan_M),nan_kn)) %10

% v,kp
vkp=length(intersect(intersect(intersect(intersect(nan_a,v),kp),nan_M),nan_kn)) %9
% v,kn
vkn=length(intersect(intersect(intersect(intersect(nan_a,v),kn),nan_M),nan_kp)) %20

% v,M
vM=length(intersect(intersect(intersect(intersect(v,nan_a),nan_kp),M),nan_kn)) %8

% a,kp
akp=length(intersect(intersect(intersect(intersect(nan_v,a),kp),nan_M),nan_kn)) %7
% a,kn
akn=length(intersect(intersect(intersect(intersect(nan_v,a),kn),nan_M),nan_kp)) %18

% a,M
aM=length(intersect(intersect(intersect(intersect(nan_v,a),nan_kp),M),nan_kn)) %6

% kp,M
kpM=length(intersect(intersect(intersect(intersect(nan_v,nan_a),kp),M),nan_kn)) %7
% kn,M
knM=length(intersect(intersect(intersect(intersect(nan_v,nan_a),kn),M),nan_kp)) %18

v_ref=length(intersect(intersect(intersect(intersect(v,nan_a),nan_kn),nan_M),nan_kp))
a_ref=length(intersect(intersect(intersect(intersect(nan_v,a),nan_kn),nan_M),nan_kp))
M_ref=length(intersect(intersect(intersect(intersect(nan_v,nan_a),nan_kn),M),nan_kp))
kp_ref=length(intersect(intersect(intersect(intersect(nan_v,nan_a),nan_kn),nan_M),kp))
kn_ref=length(intersect(intersect(intersect(intersect(nan_v,nan_a),kn),nan_M),nan_kp))