%% Saturation

clear
load('Result2.mat','EC','model_irr','Z','kcat_BRENDA_max','kmax_vivo_light','best_row_ss_light','best_col_co_light')
for i=1:length(EC)
    RN_table{i,1}=strjoin(unique(model_irr.rxnNames(find(Z(i,:)==1)))',',');
end

kmax_vivo_light_best=kmax_vivo_light{best_row_ss_light,best_col_co_light};

load('Result1.mat','kapp_Vmax_combined_all_max')

sum((kmax_vivo_light_best>kcat_BRENDA_max))

sk = kmax_vivo_light_best./kapp_Vmax_combined_all_max;

id = find(~isnan(sk));

% Supplementary Table S7a
table(EC(id),RN_table(id),sk(id));

RXN_Formula=[];RXN_EC=[];RXN_table=[];
for i=1:length(id)
    rxns=find(Z(id(i),:)==1);
    for j=1:length(rxns)
        disp(EC{id(i)})
        disp(RN_table{id(i)})
        substrates=strjoin(model_irr.mets(model_irr.S(:,rxns(j))<0)','+');
        coeff_s=model_irr.S(model_irr.S(:,rxns(j))<0,rxns(j))
        products=strjoin(model_irr.mets(model_irr.S(:,rxns(j))>0)','+');
        coeff_p=model_irr.S(model_irr.S(:,rxns(j))>0,rxns(j))
        RXN_Formula{end+1,1}=strcat(substrates,'->',products);
        disp(RXN_Formula{end})
        RXN_EC{end+1,1}=EC{id(i)};
        RXN_table{end+1,1}=RN_table{id(i)};
    end
end

table(RXN_EC,RXN_table,sk)