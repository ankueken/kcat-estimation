%% Analysis kapp and plant kcat for varying carboxylation to oxygenation and starch to sucrose ratios 

clear
%% Load model and set physiological constraints
%% 
% * AraCORE model 
% * photoautotrophic growth
% * unit: umol/gDW/d

model = readCbModel('Data/AraCORE_Arnold/ArabidopsisCoreModel.xml');
model.rxnGeneMat = zeros(length(model.rxns),length(model.genes));
for i=1:length(model.genes)
    for j=1:length(model.rxns)
        if ~isempty(strfind(model.rules{j},strcat('x(',num2str(i),')')))
            model.rxnGeneMat(j,i)=1;
        end
    end
end

%% bounds for photoautotrophic growth
bounds = readtable('Data/AraCORE_Arnold/AraCORE_photoautotrophic_bounds.csv');
model.lb = bounds{:,2};model.ub = bounds{:,3};
model.rev = zeros(size(model.rxns));
model.rev(model.lb<0) = 1;

model_irr = split_rxns(model);
model_irr.c(:) = 0;

model_irr.csense=repmat('E',size(model_irr.b));
lb_orig = model_irr.lb;
ub_orig = model_irr.ub;

%%
DIRECTION=ones(size(model_irr.rxns));
for i=1:length(model_irr.rxns)
    REV=strcmp(model_irr.rxns{i}(end-3:end),'_rev');
    if REV==1
        DIRECTION(i)=-1;
    end
end

[EC,~,IC] = unique(model_irr.rxnECNumbers);

Z=zeros(0,size(model_irr.S,2));
for i=1:length(EC)
    Z(end+1,IC==i) = DIRECTION(IC==i);
    T=unique(cell2table(model_irr.subSystems(IC==i)));
    Pathways{i,1} = T{:,1};
end

for i=1:length(EC)
    enzyme_gene_matrix(i,:)=max(model_irr.rxnGeneMat(IC==i,:),[],1);
end

%% for those with available abundance which EC have multiple rxns assigned
multiple_rxns=sum(Z')'>1;

for i=1:length(EC)
    id=find(Z(i,1:379)~=0); % exclude transport rxns, biomass
    
    comp_temp=[];
    for j=1:length(id)
        comp_temp{j}=model_irr.rxns{id(j)}(end);
    end
    Comp.names{i} = strjoin(unique(comp_temp),',');
    Comp.number(i) = length(unique(comp_temp));
end
%% *Functions and Packages needed*
% Matlab Optimization Toolbox
% 
% Cobra Toolbox 3.0
%% Set solver preference

changeCobraSolver('glpk');
options = optimset('linprog');
options.Display = 'off';
%% Choose c/o and s/s range

co_ratio=1:0.2:4;
ss_ratio=1:0.2:3;
%% Load Results

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
                kcat_BRENDA_mean(i,1) = mean(kcat_temp)*60;
                kcat_BRENDA_min(i,1) = min(kcat_temp)*60;
                kcat_BRENDA_max(i,1) = max(kcat_temp)*60;
                kcat_BRENDA_taxa(i,1) = max(t_temp);
            else
                kcat_BRENDA_mean(i,1) = nan;
                kcat_BRENDA_min(i,1) = nan;
                kcat_BRENDA_max(i,1) = nan;
            end
        elseif ~isempty(strfind(EC{i},'+'))
            EC_temp=strsplit(EC{i},' + ');
            kcat_temp=[];t_temp=[];
            for j=1:length(EC_temp)
                kcat_temp=[kcat_temp; kcat.kcat_1_sec(find(strcmp(kcat.EC,EC_temp(j))))];
                t_temp=[t_temp; kcat.tax_dist(find(strcmp(kcat.EC,EC_temp(j))))];
            end
            if ~isempty(kcat_temp)
                kcat_BRENDA_mean(i,1) = mean(kcat_temp)*60;
                kcat_BRENDA_min(i,1) = min(kcat_temp)*60;
                kcat_BRENDA_max(i,1) = max(kcat_temp)*60;
                kcat_BRENDA_taxa(i,1) = max(t_temp);
            else
                kcat_BRENDA_mean(i,1) = nan;
                kcat_BRENDA_min(i,1) = nan;
                kcat_BRENDA_max(i,1) = nan;
            end
        else
            kcat_BRENDA_mean(i,1) = nan;
            kcat_BRENDA_max(i,1) = nan;
            kcat_BRENDA_min(i,1) = nan;
        end
    else
        kcat_BRENDA_mean(i,1) = mean(kcat.kcat_1_sec(find(strcmp(kcat.EC,EC(i)))));
        kcat_BRENDA_min(i,1) = min(kcat.kcat_1_sec(find(strcmp(kcat.EC,EC(i)))));
        kcat_BRENDA_max(i,1) = max(kcat.kcat_1_sec(find(strcmp(kcat.EC,EC(i)))));
        kcat_BRENDA_taxa(i,1) = min(kcat.tax_dist(find(strcmp(kcat.EC,EC(i)))));
    end
end

kcat_BRENDA_max(kcat_BRENDA_taxa>18)=nan;
kcat_BRENDA_min(kcat_BRENDA_taxa>18)=nan;
kcat_BRENDA_mean(kcat_BRENDA_taxa>18)=nan;
kcat_BRENDA_median(kcat_BRENDA_taxa>18)=nan;

load('Result2_kcat_from_vmax_from_BRENDA.mat','gDW_per_gFW','EA','EA2','EA3','iu*','kcat_estimated_vmax_combined_all_mean','abundance_per_EC_mg_gDW_Pyl','abundance_per_EC_mg_gFW_Piques','abundance_per_EC_mg_gFW_Seaton')
taxa=readtable('Data/taxaWB.csv','Delimiter','\t');

files=strcat('MW_',strrep(strrep(strrep(EC,'.','_'),'/','-'),' + ','-'),'.csv');

for f=2:length(files)
if exist(strcat('Data/MW/',files{f}))
    temp = readtable(strcat('Data/MW/',files{f}),'Delimiter','\t');
    temp.Var9=nan(size(temp.Var8));
    Org = unique(temp{:,6});
    for o=1:length(Org)
        if ~isempty(find(strcmp(taxa.Var3,Org(o)))) && ~strcmp(Org(o),'')
            temp.Var9(strcmp(temp.Var6,Org(o)))=taxa.Var1(find(strcmp(taxa.Var3,Org(o))));
        end
    end
    temp.Var3(temp.Var3==-999)=nan;
    MW_mg_per_nmol(f,1) = (median(temp.Var3(temp.Var9==min(temp.Var9,[],1,'omitnan')))*1000)/1e+9; % from Da to mg/nmol
else
    MW_mg_per_nmol(f,1) = nan;
end
end
MW_mg_per_nmol(isnan(MW_mg_per_nmol))=median(MW_mg_per_nmol,"all","omitnan");

kcat_BRENDA_max_vmax=kcat_BRENDA_max;
kcat_BRENDA_max_vmax(isnan(kcat_estimated_vmax_combined_all_mean))=nan;
kcat_BRENDA_mean_vmax=kcat_BRENDA_mean;
kcat_BRENDA_mean_vmax(isnan(kcat_estimated_vmax_combined_all_mean))=nan;

kcat_estimated_vmax_combined_all_mean_vmax=kcat_estimated_vmax_combined_all_mean;
kcat_estimated_vmax_combined_all_mean_vmax(isnan(kcat_BRENDA_max_vmax))=nan;
%% Estimate k_vivo_max
% Fixed light, CO2
%% 
% * set range for light and CO2 uptake according to measurements
% * we use measured enzyme activity to bound flux (v<=v_max)
% * condition specific v_max sampled from corresponding measured range (Piques) 
% or sampled from range v_max+-20% (Pyl) 
%% 
% Fixed biomass
%% 
% * sample biomass values around (Pyl, condition-specific measured biomass +- 
% 20%) or within (Piques, value from measured range) measurements
% * we use condition-specific measured enzyme activity to bound flux (v<=v_max_k)
% * SA_vivo_max for different c/o and s/s ratios

carb = find(strcmp(model_irr.rxns,'RBC_h'));
oxy = find(strcmp(model_irr.rxns,'RBO_h'));
starch1 = find(strcmp(model_irr.rxns,'StS_h1'));
starch2 = find(strcmp(model_irr.rxns,'StS_h2'));
starch3 = find(strcmp(model_irr.rxns,'StS_h3'));
suc = find(strcmp(model_irr.rxns,'SucS_c'));
suc_syn = find(strcmp(model_irr.rxns,'SucS_c_rev'));

%% max enzyme activity
% EA - Pyl, EA2 - Piques
EA.data_umol_gFW_d = (EA.data/1000)*1440; % nmol/gFW/min -> umol/gFW/d
EA.data_umol_gDW_d = EA.data_umol_gFW_d./repmat(gDW_per_gFW',1,size(EA.data,2)); % umol/gFW/d -> umol/gDW/d

EA2.data_umol_gFW_d = (EA2.data/1000)*1440; % nmol/gFW/min -> umol/gFW/d
EA2.data_umol_gDW_d = EA2.data_umol_gFW_d/0.088; % umol/gFW/d -> umol/gDW/d
EA2.data_umol_gFW_d_lb = (EA2.data_lb/1000)*1440; % nmol/gFW/min -> umol/gFW/d
EA2.data_umol_gDW_d_lb = EA2.data_umol_gFW_d_lb/0.088; % umol/gFW/d -> umol/gDW/d
EA2.data_umol_gFW_d_ub = (EA2.data_ub/1000)*1440; % nmol/gFW/min -> umol/gFW/d
EA2.data_umol_gDW_d_ub = EA2.data_umol_gFW_d_ub/0.088; % umol/gFW/d -> umol/gDW/d

EA3.data_umol_gFW_d = (EA3.data/1000)*1440; % nmol/gFW/min -> umol/gFW/d
EA3.data_umol_gDW_d = EA3.data_umol_gFW_d./0.088; % umol/gFW/d -> umol/gDW/d

abundance_per_EC_mg_gDW_Seaton=abundance_per_EC_mg_gFW_Seaton/0.088;

model_orig=model_irr;

%% max enzyme activity
v_max_max=nan(length(EC),10);
v_max_max(iu_EC2,1:5)=EA.data_umol_gDW_d(:,iu_Pyl)';
v_max_max(iu_EC1,6)=EA2.data_umol_gDW_d(iu_Piques);
v_max_max(iu_EC3,7:10)=EA3.data_umol_gDW_d(iu_Seaton,:);
v_max_max=max(v_max_max,[],2,"omitnan");

col=0;
for co_ratio_max=co_ratio
    col=col+1;row=0;
    for ss_ratio_max=ss_ratio
        row=row+1;v=[];k_vivo=[];
        
        model_irr=model_orig;
        
        model_irr.S(end+1,[carb oxy]) = [1 -co_ratio_max];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'c/o ratio ub'};
        model_irr.metNames(end+1) = {'c/o ratio ub'};
        
        model_irr.S(end+1,[starch1 starch2 starch3 suc suc_syn]) = [1 1 1 ss_ratio_max -ss_ratio_max];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'starch/suc ratio ub'};
        model_irr.metNames(end+1) = {'starch/suc ratio ub'};
        
        %% program 1 - Pyl
        % carbon limiting biomass
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;
        model_irr.c(678)=1;
        model_irr.ub(678)=1000;
        
        id=find(~isnan(v_max_max));
        
        model_c = model_irr;
        model_c.S = [model_irr.S;
            Z(id,:)];
        model_c.b = [model_irr.b; v_max_max(id)];
        model_c.csense = [model_irr.csense; repmat('L',length(id),1);];
        for i=1:size(model_c.S,1)
            model_c.mets{i,1} = strcat('M',num2str(i));
        end
        for i=1:size(model_c.S,2)
            model_c.rxns{i,1} = strcat('R',num2str(i));
        end
        [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model_c.c,model_c.S(model_c.csense=='L',:),model_c.b(model_c.csense=='L'),model_c.S(model_c.csense=='E',:),model_c.b(model_c.csense=='E'),model_c.lb,model_c.ub,options);
        bio(row,col)=abs(Sol.f);
    end
end

ref_bio=max(max(bio,[],1,'omitnan'));
ref_bio=ref_bio+(0.1*ref_bio);
k_max_vivo=[];Diff_flux=[];col=0;v_min=[];v_max=[];
for co_ratio_max=co_ratio
    col=col+1;row=0;
    for ss_ratio_max=ss_ratio
        row=row+1;v=[];k_vivo=[];
        
        model_irr=model_orig;
        
        model_irr.S(end+1,[carb oxy]) = [1 -co_ratio_max];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'c/o ratio ub'};
        model_irr.metNames(end+1) = {'c/o ratio ub'};
        
        model_irr.S(end+1,[starch1 starch2 starch3 suc suc_syn]) = [1 1 1 ss_ratio_max -ss_ratio_max];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'starch/suc ratio ub'};
        model_irr.metNames(end+1) = {'starch/suc ratio ub'};
        
        %% program 1 - Pyl
        % carbon limiting biomass
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;
        
        bio = ([0.15 0.2 0.23 0.25 0.24]/0.3065)*ref_bio;
        
        [~,id,id2]=intersect(EC,EA.EC);
        
        for k=1:size(abundance_per_EC_mg_gDW_Pyl,2)
            
            v_max_sampled = EA.data_umol_gDW_d(k,id2)';
            
            model_irr.lb(678)=bio(k)-(bio(k)*0.2);
            model_irr.ub(678)=bio(k)+(bio(k)*0.2);
%             
            model_c = model_irr;
            model_c.S = [model_irr.S;
                Z(id,:)];
            model_c.b = [model_irr.b; v_max_sampled];
            model_c.csense = [model_irr.csense; repmat('L',length(id2),1);];
            for i=1:size(model_c.S,1)
                model_c.mets{i,1} = strcat('M',num2str(i));
            end
            for i=1:size(model_c.S,2)
                model_c.rxns{i,1} = strcat('R',num2str(i));
            end
            gene_associated_rxns=find(sum(model_irr.rxnGeneMat(:,find(cellfun(@isempty,strfind(model_irr.genes,'AT'))==0))')>0);
            model_c.c(gene_associated_rxns)=-1;
            [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model_c.c,model_c.S(model_c.csense=='L',:),model_c.b(model_c.csense=='L'),model_c.S(model_c.csense=='E',:),model_c.b(model_c.csense=='E'),model_c.lb,model_c.ub,options);
             
            if ~isempty(Sol.x) && Sol.stat==1
                v(:,end+1) = Sol.x;
                v(v<1e-8)=0;
                [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_c,100,'max');
                Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);
                k_vivo(:,end+1) = abs(Z*v(1:680,end))./abundance_per_EC_mg_gDW_Pyl(:,k);
            end
        end
        
        %% program 2 - Piques
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;
        
        bio2=(0.175/0.3065)*ref_bio;
        
        [~,id,id2]=intersect(EC,EA2.EC);
        
        v_max_sampled = EA2.data_umol_gDW_d(id2)';
        model_irr.lb(680) = bio2-(bio2*0.2);
        model_irr.ub(680) = bio2+(bio2*0.2);
%         
        model_c = model_irr;
        model_c.S = [model_irr.S;
            Z(id,:)];
        model_c.b = [model_irr.b; v_max_sampled];
        model_c.csense = [model_irr.csense; repmat('L',length(id2),1);];
        for i=1:size(model_c.S,1)
            model_c.mets{i,1} = strcat('M',num2str(i));
        end
        for i=1:size(model_c.S,2)
            model_c.rxns{i,1} = strcat('R',num2str(i));
        end
        gene_associated_rxns=find(sum(model_irr.rxnGeneMat(:,find(cellfun(@isempty,strfind(model_irr.genes,'AT'))==0))')>0);
        model_c.c(gene_associated_rxns)=-1;
        [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model_c.c,model_c.S(model_c.csense=='L',:),model_c.b(model_c.csense=='L'),model_c.S(model_c.csense=='E',:),model_c.b(model_c.csense=='E'),model_c.lb,model_c.ub,options);
            
        if ~isempty(Sol.x) && Sol.stat==1
            v(:,end+1) = Sol.x;
            v(v<1e-8)=0;
            [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_c,100,'max');
            Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);            
            k_vivo(:,end+1) = abs(Z*v(1:680,end))./(abundance_per_EC_mg_gFW_Piques/0.088);
        end
        
        %% program 3 - Seaton
        % carbon limiting biomass
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;
        
        bio = ([0.1135 0.1708 0.2600 0.3065]/0.3065)*ref_bio;
        
        [~,id,id2]=intersect(EC,EA3.EC);
        
        for k=1:size(abundance_per_EC_mg_gDW_Seaton,2)
            
            v_max_sampled = EA3.data_umol_gDW_d(id2,k);
            
            if k==4
                model_irr.lb(680)=bio(k)-(bio(k)*0.2);
                model_irr.ub(680)=bio(k)+(bio(k)*0.2);
            elseif k==3 
                model_irr.lb(680)=bio(k)-(bio(k)*0.2);
                model_irr.ub(680)=bio(k)+(bio(k)*0.2);
            else
                model_irr.lb(678)=bio(k)-(bio(k)*0.2);
                model_irr.ub(678)=bio(k)+(bio(k)*0.2);
            end
            
            model_c = model_irr;
            model_c.S = [model_irr.S;
                Z(id,:)];
            model_c.b = [model_irr.b; v_max_sampled];
            model_c.csense = [model_irr.csense; repmat('L',length(id2),1);];
            for i=1:size(model_c.S,1)
                model_c.mets{i,1} = strcat('M',num2str(i));
            end
            for i=1:size(model_c.S,2)
                model_c.rxns{i,1} = strcat('R',num2str(i));
            end
            gene_associated_rxns=find(sum(model_irr.rxnGeneMat(:,find(cellfun(@isempty,strfind(model_irr.genes,'AT'))==0))')>0);
            model_c.c(gene_associated_rxns)=-1;
            [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(-model_c.c,model_c.S(model_c.csense=='L',:),model_c.b(model_c.csense=='L'),model_c.S(model_c.csense=='E',:),model_c.b(model_c.csense=='E'),model_c.lb,model_c.ub,options);
             
            if ~isempty(Sol.x) && Sol.stat==1
                v(:,end+1) = Sol.x;
                v(v<1e-8)=0;
                [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_c,100,'max');
                Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);                
                k_vivo(:,end+1) = abs(Z*v(1:680,end))./abundance_per_EC_mg_gDW_Seaton(:,k);
            end
        end        
        
        if ~isempty(v)
            [k_max_vivo(:,end+1),I] = max(k_vivo,[],2,'omitnan');
            for i=1:length(I)
                vDirection{row,col}(i,1)=Z(i,:)*v(:,I(i));
            end
            k_max_vivo(k_max_vivo==0)=nan;
            
            MW_mg_per_umol = MW_mg_per_nmol*1000;
            k_cat_vivo = k_max_vivo(:,end).*MW_mg_per_umol;
            
            var1 = k_cat_vivo/86400;
            
            [C_mean(row,col),P_mean(row,col)] = corr(var1,kcat_estimated_vmax_combined_all_mean,'rows','pairwise');
            [C_mean_log(row,col),P_mean_log(row,col)] = corr(log(var1),log(kcat_estimated_vmax_combined_all_mean),'rows','pairwise');
           
            [C_mean_log_vmax(row,col),P_mean_log_vmax(row,col)] = corr(log(var1),log(kcat_estimated_vmax_combined_all_mean_vmax),'rows','pairwise');
            
            [C_mean_log_B_max_vmax(row,col),P_mean_log_B_max_vmax(row,col)] = corr(log(var1),log(kcat_BRENDA_max_vmax),'rows','pairwise');
            [C_mean_log_B_mean_vmax(row,col),P_mean_log_B_mean_vmax(row,col)] = corr(log(var1),log(kcat_BRENDA_mean_vmax),'rows','pairwise');

            [C_mean_B_max(row,col),P_mean_B_max(row,col)] = corr(var1,kcat_BRENDA_max,'rows','pairwise');
            [C_mean_log_B_max(row,col),P_mean_log_B_max(row,col)] = corr(log(var1),log(kcat_BRENDA_max),'rows','pairwise');
            [C_mean_B_mean(row,col),P_mean_B_mean(row,col)] = corr(var1,kcat_BRENDA_mean,'rows','pairwise');
            [C_mean_log_B_mean(row,col),P_mean_log_B_mean(row,col)] = corr(log(var1),log(kcat_BRENDA_mean),'rows','pairwise');
            [C_mean_B_min(row,col),P_mean_B_min(row,col)] = corr(var1,kcat_BRENDA_min,'rows','pairwise');
            [C_mean_log_B_min(row,col),P_mean_log_B_min(row,col)] = corr(log(var1(kcat_BRENDA_min>0)),log(kcat_BRENDA_min(kcat_BRENDA_min>0)),'rows','pairwise');
            
            Var1{row,col}=var1;
        end
    end
end
%% RESULT 3 - Correlation kapp and kcat from vmax

disp('Correlation|P-Value|Number')
[r_mean_log,c_mean_log]=find(C_mean_log_vmax==max(C_mean_log_vmax(:),[],1,'omitnan'));
[C_mean_log_vmax(r_mean_log,c_mean_log) P_mean_log_vmax(r_mean_log,c_mean_log)...
    length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_estimated_vmax_combined_all_mean_vmax))))]

C_mean_log_B_max(:,1)=nan;C_mean_log_B_min(:,1)=nan;C_mean_log_B_mean(:,1)=nan;C_mean_log(:,1)=nan;C_mean_log_B_mean_vmax(:,1)=nan;

[C_mean_log_B_max_vmax(r_mean_log,c_mean_log) P_mean_log_B_max_vmax(r_mean_log,c_mean_log) length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_max_vmax))));...
    C_mean_log_B_mean_vmax(r_mean_log,c_mean_log) P_mean_log_B_mean_vmax(r_mean_log,c_mean_log) length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_mean_vmax))));...
    C_mean_log_B_max(r_mean_log,c_mean_log) P_mean_log_B_max(r_mean_log,c_mean_log) length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_max))));...
    C_mean_log_B_mean(r_mean_log,c_mean_log) P_mean_log_B_mean(r_mean_log,c_mean_log) length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_mean))));...
    C_mean_log_B_min(r_mean_log,c_mean_log) P_mean_log_B_min(r_mean_log,c_mean_log) length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_min))))]

disp('kapp<=kcat')
sum(Var1{r_mean_log,c_mean_log}<kcat_BRENDA_mean)/length(intersect(find(~isnan(Var1{r_mean_log,c_mean_log})),find(~isnan(kcat_BRENDA_max))))
%% 
% 
figure
subplot(2,2,1)
heatmap(C_mean_log_vmax,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',0,'MinColorValue',0,'MaxColorValue',0.75,'NaNColor',[1 1 1],'GridLines', ':')
% xlabel('carboxylation/oxygenation ratio')
ylabel('starch/sucrose synthases')
title('k_{app} and k_{cat}')

subplot(2,2,2)
heatmap(C_mean_log_B_mean_vmax,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',1,'MinColorValue',0,'MaxColorValue',0.75,'NaNColor',[1 1 1],'GridLines', ':')
% xlabel('carboxylation/oxygenation ratio')
% ylabel('starch/sucrose synthase ratio')
title('k_{app} and mean k_{cat} for enzymes with available V_{max}')

subplot(2,2,3)
heatmap(C_mean_log_B_mean,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',0,'MinColorValue',0,'MaxColorValue',0.75,'NaNColor',[1 1 1],'GridLines', ':')
xlabel('carboxylation/oxygenation')
ylabel('starch/sucrose synthases')
title('k_{app} and mean k_{cat}')

subplot(2,2,4)
heatmap(C_mean_log_B_max,co_ratio,ss_ratio,[],'Colormap',jet,'ColorBar',1,'MinColorValue',0,'MaxColorValue',0.75,'NaNColor',[1 1 1],'GridLines', ':')
xlabel('carboxylation/oxygenation')
% ylabel('starch/sucrose synthase ratio')
title('k_{app} and maximum k_{cat}')

%% 
% 
figure
subplot(1,3,1)
plot(kcat_estimated_vmax_combined_all_mean_vmax,Var1{r_mean_log,c_mean_log},'.');
xlabel('k_{cat} [s^{-1}]')
ylabel('k_{app} [s^{-1}]')
set(gca,'XScale','Log','YScale','Log')
hold on
plot([1e-4 1e7],[1e-4 1e7],'k')
hold on
[~,m,b]=regression(log(kcat_estimated_vmax_combined_all_mean_vmax)',log(Var1{r_mean_log,c_mean_log})');
x=log([1e-4; 1e-3; 1e-2; 1e-1;1e0; 1e1; 1e2;1e3;1e4;1e5;1e6; 1e7]);
y=(m.*x)+b;
plot(exp(x),exp(y),'-b')
xlim([1e-4 1e7])
ylim([1e-4 1e7])

subplot(1,3,2)
plot(kcat_BRENDA_mean,Var1{r_mean_log,c_mean_log},'.');
xlabel('k_{cat} [s^{-1}]')
% ylabel('k_{app}')
set(gca,'XScale','Log','YScale','Log')
hold on
xlim([1e-4 1e7])
ylim([1e-4 1e7])
plot([1e-4 1e7],[1e-4 1e7],'k')
hold on
[~,m,b]=regression(log(kcat_BRENDA_mean)',log(Var1{r_mean_log,c_mean_log})');
x=log([1e-4; 1e-3; 1e-2; 1e-1;1e0; 1e1; 1e2;1e3;1e4;1e5;1e6; 1e7]);
y=(m.*x)+b;
plot(exp(x),exp(y),'-b')

x=log(Var1{r_mean_log,c_mean_log}(kcat_BRENDA_mean>0));
y=log(kcat_BRENDA_mean(kcat_BRENDA_mean>0));

subplot(1,3,3)
plot(exp(y(Comp.number(kcat_BRENDA_mean>0)==1)),exp(x(Comp.number(kcat_BRENDA_mean>0)==1)),'.');
xlabel('k_{cat} [s^{-1}]')
% ylabel('k_{app}')
set(gca,'XScale','Log','YScale','Log')
hold on
xlim([1e-4 1e7])
ylim([1e-4 1e7])
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
% 
disp('Flux variability in flux estimates:')
figure
x=max((Diff_flux./1000)*100,[],2,'omitnan');
[H,b]=hist(x(:),50);
bar(b,(H/sum(H)),1)
xlabel('Flux variability with respect to maximum flux [%]')
ylabel('Fraction')
%% 
% 
save('Result3_kcat_from_flux_kcat_from_BRENDA.mat')