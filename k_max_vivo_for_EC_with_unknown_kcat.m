clear
load('Result2.mat')

% kcat unknown
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

kmax_vivo_light_best=kmax_vivo_light{best_row_ss_light,best_col_co_light};
%% enzymes with previously unknown kcat

unknown_kcat_id = intersect(find(isnan(kcat_BRENDA_mean)),find(~isnan(kmax_vivo_light_best)));

for i=1:length(unknown_kcat_id)
    RN_table{i,1}=strjoin(unique(model_irr.rxnNames(find(Z(unknown_kcat_id(i),:)==1)))',',');
end

for i=1:length(Pathways)
    Pathways{i,1}=strjoin(unique(Pathways{i})',',');
end

table(EC(unknown_kcat_id),RN_table,kmax_vivo_light_best(unknown_kcat_id),Pathways(unknown_kcat_id))
