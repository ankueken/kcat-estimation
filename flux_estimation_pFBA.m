%% flux distributions fixed biomass and pFBA
load_model
% unit: mol/gDW/day

% Set solver preference
changeCobraSolver('glpk');
options = optimset('linprog');
options.Display = 'off';

v=[];k_vivo=[];kmax_vivo=[];v_min=[];v_max=[];Diff_flux=[];
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
        v(:,end+1) = Sol.x;
        v(v<1e-8)=0;
        
        [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_irr,100,'min');
        Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);         
        
        k_vivo(:,end+1) = abs(Z*v(1:680,end))./abundance_per_EC_mg_gDW_Pyl(:,k);
    else
        v(1:length(model_irr.rxns),end+1) = nan;
        k_vivo(1:223,end+1) = nan;
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
    v(:,end+1) = Sol.x;
    v(v<1e-8)=0;
    
    [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_irr,100,'min');
    Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);
        
    k_vivo(:,end+1) = abs(Z*v(1:680,end))./(abundance_per_EC_mg_gDW_Piques);
else
    v(1:length(model_irr.rxns),end+1) = nan;
    k_vivo(1:223,end+1) = nan;
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
        v(:,end+1) = Sol.x;
        v(v<1e-8)=0;
        
        [v_min(:,end+1),v_max(:,end+1)] = fluxVariability(model_irr,100,'min');
        Diff_flux(:,end+1) = v_max(:,end) - v_min(:,end);
        
        k_vivo(:,end+1) = abs(Z*v(1:680,end))./abundance_per_EC_mg_gDW_Seaton(:,k);
    else
        v(1:length(model_irr.rxns),end+1) = nan;
        k_vivo(1:223,end+1) = nan;
    end
end
        
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

if ~isempty(v)
    MW_mg_per_umol = MW_mg_per_nmol*1000;
    v_umol_gDW_day = abs(Z*v(1:680,:)).*MW_mg_per_umol;
    
    v_umol_gDW_sec = v_umol_gDW_day/86400;    
end
