%% LOAD MODEL
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
model.lb = bounds{:,2}*1000;model.ub = bounds{:,3}*1000;
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

%% for those enzymes with available abundance which have multiple rxns assigned
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

