%% activity data

% unit: nmol/gFW/min
% molecular weight (mg/mol -> mg/nmol) from BRENDA to transform nmol/gFW/min into mg/gFW/min

%% Sulpice et al. (same conditions as Seaton et al. abundance data)

temp_EA = readtable('Data/Sulpice_vmax/enzyme_activity.xlsx');

EA_Sulpice.data=temp_EA{:,4:7}; % nmol/gFW/min  6,8,12,18-h

EA_Sulpice.std=temp_EA{:,9:end};

EA_Sulpice.data_lb = EA_Sulpice.data - EA_Sulpice.std;
EA_Sulpice.data_ub = EA_Sulpice.data + EA_Sulpice.std;

EA_Sulpice.EC=temp_EA.EC;
EA_Sulpice.names=temp_EA.PhotoperiodMean;

% find position of EC in EA_Sulpice.EC in entire EC list
[EC_fit_Seaton,id_EC,id_data_Seaton]=intersect(EC,EA_Sulpice.EC);

% change unit from nmol/gFW/min into mg/gFW/min 
clear files
for f=1:length(EC_fit_Seaton)
    files(f,1)=strcat('MW_',strrep(EC_fit_Seaton(f),'.','_'),'.csv');
end

taxa=readtable('Data/taxaWB.csv','Delimiter','\t');

for f=1:length(files)
    temp = readtable(strcat('Data/MW/',files{f}),'Delimiter','\t');
    temp.Var9=nan(size(temp.Var8));
    Org = unique(temp{:,6});
    for o=1:length(Org)
        if ~isempty(find(strcmp(taxa.Var3,Org(o)))) && ~strcmp(Org(o),'')
            temp.Var9(strcmp(temp.Var6,Org(o)))=taxa.Var1(find(strcmp(taxa.Var3,Org(o))));
        end
    end
    temp.Var3(temp.Var3==-999)=nan;
    MW_mg_per_nmol_Seaton(f,1) = (median(temp.Var3(temp.Var9==min(temp.Var9,[],1,'omitnan')))*1000)/1e+9; % from Da to mg/nmol
    % 1Da = 1 g/mol = 1e-6 mg/nmol
end

vmax_Seaton = nan(length(EC),size(EA_Sulpice.data,2));
vmax_Seaton(id_EC,:) = EA_Sulpice.data(id_data_Seaton,:).*repmat(MW_mg_per_nmol_Seaton,1,size(EA_Sulpice.data,2)); % nmol/gFW/min * mg/nmol = mg/gFW/min

vmax_Seaton_std = nan(length(EC),size(EA_Sulpice.data,2));
vmax_Seaton_std(id_EC,:) = EA_Sulpice.std(id_data_Seaton,:).*repmat(MW_mg_per_nmol_Seaton,1,size(EA_Sulpice.data,2)); % nmol/gFW/min * mg/nmol = mg/gFW/min

%% Piques et al.

temp_EA = readtable('Data/Piques_Data/msb200968-s4_enzyme_activity.csv','Delimiter',';','HeaderLines',2);

% data after 20h
EA_Piques.data = temp_EA.mean_nmol_min_1_g_1FW_TimeDuringTheDiurnalCycle_h_20(1:35);
EA_Piques.std = temp_EA.sd_20(1:35);

EA_Piques.data_lb = EA_Piques.data - EA_Piques.std;
EA_Piques.data_ub = EA_Piques.data + EA_Piques.std;

EA_Piques.EC = temp_EA.Var2(1:35);
EA_Piques.names = temp_EA.Enzyme_name(1:35);

% find position of EC in EA_Piques.EC in entire EC list
[EC_fit_Piques,id_EC,id_data_Piques]=intersect(EC,EA_Piques.EC);

% change unit from nmol/gFW/min into mg/gFW/min 
clear files
for f=1:length(EC_fit_Piques)
    files(f,1)=strcat('MW_',strrep(EC_fit_Piques(f),'.','_'),'.csv');
end

taxa=readtable('Data/taxaWB.csv','Delimiter','\t');

for f=1:length(files)
    temp = readtable(strcat('Data/MW/',files{f}),'Delimiter','\t');
    temp.Var9=nan(size(temp.Var8));
    Org = unique(temp{:,6});
    for o=1:length(Org)
        if ~isempty(find(strcmp(taxa.Var3,Org(o)))) && ~strcmp(Org(o),'')
            temp.Var9(strcmp(temp.Var6,Org(o)))=taxa.Var1(find(strcmp(taxa.Var3,Org(o))));
        end
    end
    temp.Var3(temp.Var3==-999)=nan;
    MW_mg_per_nmol_Piques(f,1) = (median(temp.Var3(temp.Var9==min(temp.Var9,[],1,'omitnan')))*1000)/1e+9; % from Da to mg/nmol
    % 1Da = 1 g/mol = 1e-6 mg/nmol
end

vmax_Piques = nan(length(EC),size(EA_Piques.data,2));
vmax_Piques(id_EC,:) = EA_Piques.data(id_data_Piques,:).*repmat(MW_mg_per_nmol_Piques,1,size(EA_Piques.data,2)); % nmol/gFW/min * mg/nmol = mg/gFW/min

vmax_Piques_std = nan(length(EC),size(EA_Piques.data,2));
vmax_Piques_std(id_EC,:) = EA_Piques.std(id_data_Piques,:).*repmat(MW_mg_per_nmol_Piques,1,size(EA_Piques.data,2)); % nmol/gFW/min * mg/nmol = mg/gFW/min

%% Pyl et al. 

temp_EA = readtable('Data/Pyl_Data/tpc097188SupplementalDS_enzyme_activity.csv','Delimiter',';');

EA_Pyl.data = str2double(table2array(temp_EA(2:end,3:7)));

EA_Pyl.EC = temp_EA.Var2(2:end);
EA_Pyl.names = temp_EA{2:end,1};

% find position of EC in EA_Pyl.EC in entire EC list
[EC_fit_Pyl,id_EC,id_data_Pyl]=intersect(EC,EA_Pyl.EC);

% change unit from nmol/gFW/min into mg/gFW/min 
clear files
for f=1:length(EC_fit_Pyl)
    files(f,1)=strcat('MW_',strrep(EC_fit_Pyl(f),'.','_'),'.csv');
end

taxa=readtable('Data/taxaWB.csv','Delimiter','\t');

for f=1:length(files)
    temp = readtable(strcat('Data/MW/',files{f}),'Delimiter','\t');
    temp.Var9=nan(size(temp.Var8));
    Org = unique(temp{:,6});
    for o=1:length(Org)
        if ~isempty(find(strcmp(taxa.Var3,Org(o)))) && ~strcmp(Org(o),'')
            temp.Var9(strcmp(temp.Var6,Org(o)))=taxa.Var1(find(strcmp(taxa.Var3,Org(o))));
        end
    end
    temp.Var3(temp.Var3==-999)=nan;
    MW_mg_per_nmol_Pyl(f,1) = (median(temp.Var3(temp.Var9==min(temp.Var9,[],1,'omitnan')))*1000)/1e+9; % from Da to mg/nmol
    % 1Da = 1 g/mol = 1e-6 mg/nmol
end

vmax_Pyl = nan(length(EC),size(EA_Pyl.data,2));
vmax_Pyl(id_EC,:) = EA_Pyl.data(id_data_Pyl,:).*repmat(MW_mg_per_nmol_Pyl,1,size(EA_Pyl.data,2)); % nmol/gFW/min * mg/nmol = mg/gFW/min

vmax_Pyl_std = nan(length(EC),size(EA_Pyl.data,2));
vmax_Pyl_std(id_EC,:) = (eps*ones(size(EA_Pyl.data(id_data_Pyl,:)))).*repmat(MW_mg_per_nmol_Pyl,1,size(EA_Pyl.data,2));

