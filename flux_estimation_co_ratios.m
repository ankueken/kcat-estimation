%% flux distributions pFBA with fixed biomass, c/o and starch/sucrose ratio
load_model
% unit: mol/gDW/day

% Set solver preference
options = optimset('linprog');
options.Display = 'off';

%% MW 
%% unit
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

%% Choose c/o and s/s range

co_ratio=1:0.2:4;
ss_ratio=1:0.2:3;

v=cell(length(ss_ratio),length(co_ratio));
k_vivo=cell(length(ss_ratio),length(co_ratio));
kmax_vivo=cell(length(ss_ratio),length(co_ratio));
v_umol_gDW_day=cell(length(ss_ratio),length(co_ratio));
v_umol_gDW_sec=cell(length(ss_ratio),length(co_ratio));
v_min=cell(length(ss_ratio),length(co_ratio));
v_max=cell(length(ss_ratio),length(co_ratio));
Diff_flux=cell(length(ss_ratio),length(co_ratio));

% Fixed biomass (absolute)

carb = find(strcmp(model_irr.rxns,'RBC_h'));
oxy = find(strcmp(model_irr.rxns,'RBO_h'));
starch1 = find(strcmp(model_irr.rxns,'StS_h1'));
starch2 = find(strcmp(model_irr.rxns,'StS_h2'));
starch3 = find(strcmp(model_irr.rxns,'StS_h3'));
suc = find(strcmp(model_irr.rxns,'SucS_c'));
suc_syn = find(strcmp(model_irr.rxns,'SucS_c_rev'));

col_co=0;model_orig=model_irr;
for co = co_ratio
    col_co=col_co+1;row_ss=0;
    for ss=ss_ratio
        row_ss=row_ss+1;
        model_irr=model_orig;
        
        % set c/o and starch/surcose ratio
        model_irr.S(end+1,[carb oxy]) = [1 -co];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'c/o ratio ub'};
        model_irr.metNames(end+1) = {'c/o ratio ub'};
        
        model_irr.S(end+1,[starch1 starch2 starch3 suc suc_syn]) = [1 1 1 ss -ss];
        model_irr.b(end+1) = 0;
        model_irr.csense(end+1)='E';
        model_irr.mets(end+1) = {'starch/suc ratio ub'};
        model_irr.metNames(end+1) = {'starch/suc ratio ub'};
        
        %% program 1 - Pyl

        % carbon limiting biomass
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;

        % biomass measurements Pyl et al. [1/day]
        bio = [0.15 0.2 0.23 0.25 0.24];       

        model_irr.c(:)=0;
        gene_associated_rxns=find(sum(model_irr.rxnGeneMat(:,find(cellfun(@isempty,strfind(model_irr.genes,'AT'))==0))')>0);
        model_irr.c(gene_associated_rxns)=1;

        for k=1:size(abundance_per_EC_mg_gDW_Pyl,2)

            model_irr.lb(678)=bio(k);
            model_irr.ub(678)=bio(k);

            [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(model_irr.c,model_irr.S(model_irr.csense~='E',:),model_irr.b(model_irr.csense~='E'),model_irr.S(model_irr.csense=='E',:),model_irr.b(model_irr.csense=='E'),model_irr.lb,model_irr.ub,options);

            if ~isempty(Sol.x) && Sol.stat==1
                v{row_ss,col_co}(:,end+1) = Sol.x;
                v{row_ss,col_co}(v{row_ss,col_co}<1e-8)=0;
                
                [v_min{row_ss,col_co}(:,end+1),v_max{row_ss,col_co}(:,end+1)] = fluxVariability(model_irr,100,'min');
                Diff_flux{row_ss,col_co}(:,end+1) = v_max{row_ss,col_co}(:,end) - v_min{row_ss,col_co}(:,end);
                
                k_vivo{row_ss,col_co}(:,end+1) = abs(Z*v{row_ss,col_co}(1:680,end))./abundance_per_EC_mg_gDW_Pyl(:,k);
            else
                v{row_ss,col_co}(1:length(model_irr.rxns),end+1) = nan;
                k_vivo{row_ss,col_co}(1:223,end+1) = nan;
            end
        end

        %% program 2 - Piques
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;

        model_irr.lb(680) = 0.175; % growth [1/day]
        model_irr.ub(680) = 0.175;

        [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(model_irr.c,model_irr.S(model_irr.csense~='E',:),model_irr.b(model_irr.csense~='E'),model_irr.S(model_irr.csense=='E',:),model_irr.b(model_irr.csense=='E'),model_irr.lb,model_irr.ub,options);

        if ~isempty(Sol.x) && Sol.stat==1
            v{row_ss,col_co}(:,end+1) = Sol.x;
            v{row_ss,col_co}(v{row_ss,col_co}<1e-8)=0;
            
            [v_min{row_ss,col_co}(:,end+1),v_max{row_ss,col_co}(:,end+1)] = fluxVariability(model_irr,100,'min');
            Diff_flux{row_ss,col_co}(:,end+1) = v_max{row_ss,col_co}(:,end) - v_min{row_ss,col_co}(:,end);
                    
            k_vivo{row_ss,col_co}(:,end+1) = abs(Z*v{row_ss,col_co}(1:680,end))./(abundance_per_EC_mg_gDW_Piques);
        else
            v{row_ss,col_co}(1:length(model_irr.rxns),end+1) = nan;
            k_vivo{row_ss,col_co}(1:223,end+1) = nan;
        end

        %% program 3 - Seaton
        % carbon limiting biomass
        model_irr.lb(677:680)=0;
        model_irr.ub(677:680)=0;
        model_irr.c(:)=0;

        bio = [0.1135 0.1708 0.2600 0.3065]; % growth rate [1/day]

        for k=1:size(abundance_per_EC_mg_gDW_Seaton,2)

            if k==3 || k==4
                model_irr.lb(677:680)=0;
                model_irr.ub(677:680)=0;
                model_irr.lb(680)=bio(k);
                model_irr.ub(680)=bio(k);
            else
                model_irr.lb(677:680)=0;
                model_irr.ub(677:680)=0;
                model_irr.lb(678)=bio(k);
                model_irr.ub(678)=bio(k);
            end

            [Sol.x,Sol.f,Sol.stat,Sol.output]=linprog(model_irr.c,model_irr.S(model_irr.csense~='E',:),model_irr.b(model_irr.csense~='E'),model_irr.S(model_irr.csense=='E',:),model_irr.b(model_irr.csense=='E'),model_irr.lb,model_irr.ub,options);

            if ~isempty(Sol.x) && Sol.stat==1
                v{row_ss,col_co}(:,end+1) = Sol.x;
                v{row_ss,col_co}(v{row_ss,col_co}<1e-8)=0;
                
                [v_min{row_ss,col_co}(:,end+1),v_max{row_ss,col_co}(:,end+1)] = fluxVariability(model_irr,100,'min');
                Diff_flux{row_ss,col_co}(:,end+1) = v_max{row_ss,col_co}(:,end) - v_min{row_ss,col_co}(:,end);
         
                k_vivo{row_ss,col_co}(:,end+1) = abs(Z*v{row_ss,col_co}(1:680,end))./abundance_per_EC_mg_gDW_Seaton(:,k);
            else
                v{row_ss,col_co}(1:length(model_irr.rxns),end+1) = nan;
                k_vivo{row_ss,col_co}(1:223,end+1) = nan;
            end
        end

        if ~isempty(v{row_ss,col_co})
            MW_mg_per_umol = MW_mg_per_nmol*1000;
            v_umol_gDW_day{row_ss,col_co} = abs(Z*v{row_ss,col_co}(1:680,:)).*MW_mg_per_umol;

            v_umol_gDW_sec{row_ss,col_co} = v_umol_gDW_day{row_ss,col_co}/86400;
        end
        
    end % end ss ratio
end % end co ratio
