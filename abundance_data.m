%% enzyme abundance

% unit: mg/gFW
% transformation mg/gFW into mg/gDW 
clear
load_model

%% Seaton et al. 

A = readtable('Data/Seaton_abundance/msb177962-sup-0003-tableev1_abundance.csv','Delimiter','\t');

protein_level=[19.0863 18.9733 18.2321 17.6531]; % mg protein/gFW
gDW_per_gFW = 0.088; % no measured gDW/gFW we assume 0.088 gDW/gFW (Tschoep et al., 2009)
val = protein_level./gDW_per_gFW;

for s=1:length(protein_level)
    % relative abundance per sample (mol %)
    tot = sum(A{:,s+1},'omitnan');
    A{:,s+1} = A{:,s+1}/tot;
    
    mapECvalues=zeros(size(enzyme_gene_matrix));
    for i=1:length(model_irr.genes)
        id = find(strcmp(A.locus_id,model_irr.genes(i)));
        if ~isempty(id)
            mapECvalues(enzyme_gene_matrix(:,i)~=0,i) = sum(A{id,s+1},'omitnan');
        end
    end
    % we use total protein data to calculate mol % into mg protein/gFW
    abundance_per_EC_mg_gFW_Seaton(:,s)=sum(mapECvalues,2).*protein_level(s);
    abundance_per_EC_mg_gFW_Seaton(abundance_per_EC_mg_gFW_Seaton==0) = nan;
    % we use total protein data to calculate mol % into mg protein/gDW
    abundance_per_EC_mg_gDW_Seaton(:,s)=sum(mapECvalues,2).*val(s);
    abundance_per_EC_mg_gDW_Seaton(abundance_per_EC_mg_gDW_Seaton==0) = nan;
end

%% Piques et al.

A = readtable('Data/Piques_Data/abundance_msb200968-s7.csv','HeaderLines',3,'Delimiter','\t');

protein_level=15; % mg protein/gFW
gDW_per_gFW = 0.088; % no measured gDW/gFW we assume 0.088 gDW/gFW (Tschoep et al., 2009)
val = protein_level/gDW_per_gFW;

% data in mol %
tot = sum(A.Var7,'omitnan');
A.Var7 = A.Var7/tot;

mapECvalues=zeros(size(enzyme_gene_matrix));
for i=1:length(model_irr.genes)
    id = find(strcmp(A.Var3,model_irr.genes(i)));
    if ~isempty(id)
        mapECvalues(enzyme_gene_matrix(:,i)~=0,i) = sum(A.Var7(id),'omitnan');
    end
end
% transform mol % into mg protein/gFW
abundance_per_EC_mg_gFW_Piques=sum(mapECvalues,2).*protein_level;
abundance_per_EC_mg_gFW_Piques(abundance_per_EC_mg_gFW_Piques==0) = nan;
% transform mol % into mg protein/gDW
abundance_per_EC_mg_gDW_Piques=sum(mapECvalues,2).*val;
abundance_per_EC_mg_gDW_Piques(abundance_per_EC_mg_gDW_Piques==0) = nan;

%% Pyl et al. 

A = readtable('Data/Pyl_Data/tpc097188SupplementalDS4.csv','HeaderLines',3,'Delimiter','\t');

protein_level=[20.84 17.91 13.43 12.91 12.91]; % mg protein/gFW Suppl. DataSet2
gDW_per_gFW=[109.9 98 87.5 88.2 88.1]/1000; % gDW/gFW
val = protein_level./gDW_per_gFW; % mg protein/gDW

for f=1:size(A,2)-1
    % unit mol %
    tot = sum(A{:,f+1});
    A{:,f+1} = A{:,f+1}/tot;
    
    mapECvalues=zeros(size(enzyme_gene_matrix));
    for i=1:length(model_irr.genes)
        id = find(strcmp(A.Genes,model_irr.genes(i)));
        if ~isempty(id)
            mapECvalues(enzyme_gene_matrix(:,i)~=0,i) = sum(A{id,f+1},'omitnan');
        end
    end
    % tranform unit from mol % into mg protein/gDW
    abundance_per_EC_mg_gDW_Pyl(:,f)=sum(mapECvalues,2).*val(f);
    % tranform unit from mol % into mg protein/gFW
    abundance_per_EC_mg_gFW_Pyl(:,f)=sum(mapECvalues,2).*protein_level(f);
end
abundance_per_EC_mg_gDW_Pyl(abundance_per_EC_mg_gDW_Pyl==0) = nan;
abundance_per_EC_mg_gFW_Pyl(abundance_per_EC_mg_gFW_Pyl==0) = nan;

