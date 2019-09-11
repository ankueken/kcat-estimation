clear
taxa=readtable('Data/taxaWB.csv');

files=dir('Data/KM/*.csv');

N=cell(16,1);
for f=1:length(files)
    KM=readtable(strcat('Data/KM/', files(f).name));
    bestTaxa=25;
    KM_closest=[];
    for i=1:size(KM,1)
        if taxa.Var1(find(strcmp(taxa.Var3,KM.Var7(i))))> bestTaxa
        elseif taxa.Var1(find(strcmp(taxa.Var3,KM.Var7(i))))== bestTaxa
            KM_closest(end+1)=KM.Var3(i);
        elseif taxa.Var1(find(strcmp(taxa.Var3,KM.Var7(i))))< bestTaxa
            KM_closest=[];
            bestTaxa=taxa.Var1(find(strcmp(taxa.Var3,KM.Var7(i))));
            KM_closest(end+1)=KM.Var3(i);
        end
    end
    KM_mean(f,1)=mean(KM_closest,'all');
end

