clear
load('Result2.mat')
kmax_vivo_light_best=kmax_vivo_light{best_row_ss_light,best_col_co_light};

% number of enzymes for which k_max_vivo could be estimated and kcat known
% from plants
length(intersect(find(~isnan(kcat_BRENDA_max)),find(~isnan(kmax_vivo_light_best))))

for i=1:length(EC)
    RN_table{i,1}=strjoin(unique(model_irr.rxnNames(find(Z(i,:)==1)))',',');
end


for i=1:length(Pathways)
    Pathways{i,1}=strjoin(unique(Pathways{i})',',');
end

% largest k_max_vivo
[val,id] = sort(kmax_vivo_light_best);

table(EC(id),RN_table(id),val)

% CBC enzymes
id1=find(cellfun(@isempty,strfind(Pathways,'Calvin-Benson'))==0);
id2=find(cellfun(@isempty,strfind(Pathways,'starch synthesis'))==0);
id3=find(cellfun(@isempty,strfind(Pathways,'sucrose'))==0);

id=[id1;id2;id3];
[val,id_short]=sort(kmax_vivo_light_best(id));

table(EC(id(id_short)),RN_table(id(id_short)),kmax_vivo_light_best(id(id_short)))

% enzymes without estimated k_max_vivo, but known kcat and abundance
intersect(find(isnan(kmax_vivo_light_best)),...
    intersect(find(~all(isnan(abundance_per_EC_mg_gDW_combined_all'))),find(~isnan(kcat_BRENDA_max))))

% enzymes with largest discrapancy predicted and average measured
[discrap,id]=sort(abs((kmax_vivo_light_best./kcat_BRENDA_median)-1));
table(EC(id),RN_table(id),discrap,kmax_vivo_light_best(id),kcat_BRENDA_median(id),kcat_BRENDA_mean(id))
