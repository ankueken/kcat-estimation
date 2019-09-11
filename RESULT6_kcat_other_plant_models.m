clear
% Load models 

files=[dir('../Data/PlantModels/*.xml');dir('../Data/PlantModels/*.sbml')];

[EC,~,IC] = unique(model_irr.rxnECNumbers);

kcat=readtable('../Data/kcat_final_no_mutant.csv','Delimiter','\t'); % kcat BRENDA
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
