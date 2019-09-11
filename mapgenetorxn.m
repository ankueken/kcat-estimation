function maprxnvalues=mapgenetorxn(GEM,genenames,genevalues,modellogic)

%**********Function to map individual gene values to reactions*************
% Arguments

%GEM: Genome-scale model, COBRA structure
%genenames: cell of char entries containing the gene names corresponding to
%the supplied data values. They must be written in the same format as the
%gene names provided in the GEM
%genevalues: a vector containing the expression values, must be in the same
%order as genenames
%modellogic: Logical; modellogic=1 => A OR B, A+B; 
%modellogic=0 => A OR B, max(A,B)

%Semid�n, October, 2014
%**************************************************************************
    
%Maps gene value to gene name in model genes list
%genenames must be a cell of char entries
    
modelgenes=GEM.genes;
modelrules=GEM.rules;

if nargin<4,
    modellogic=0;
end
if isempty(modellogic),
    modellogic=0;
end

%Eliminates the last 2 characters in Recon1 genes to match the exact 
%Entrez gene name

% for i=1:length(modelgenes),
%     modelgenes{i}=modelgenes{i}(1:end-2);
% end

modeldata=[];loc=[];geneloc=[];
for n=1:length(modelgenes),
    for p=1:length(genenames),
    loc(p)=strcmpi(modelgenes{n},genenames{p});
    end
    try
        geneloc(n)=find(loc==1,1);
        modeldata(n)=genevalues(geneloc(n));
    catch
        if modellogic==1,
        modeldata(n)=0;
        elseif modellogic==0,
        modeldata(n)=nan;
        end 
    end       
end

%Map gene values to reactions according to model rules

modelgenedata={};
errorinrule=[];

for n=1:length(modelgenes)
    modelgenedata{n,1}=sprintf('x(%d)',n);
    modelgenedata{n,2}=sprintf('(x(%d))',n);
    modelgenedata{n,3}=sprintf('%d',modeldata(n));
end

    %Asigns x(i) to corresponding value

mrm=modelrules;
for n=1:length(modelgenes),
    mrm=strrep(mrm,modelgenedata(n,2),modelgenedata(n,3));   
    mrm=strrep(mrm,modelgenedata(n,1),modelgenedata(n,3)); 
end

for n=1:length(mrm)
    
    %Selects substring within last pair of brackets
    
while isempty(strfind(mrm{n},'('))==0,

    lastbracket=strfind(mrm{n},'(');
    firstbracket=strfind(mrm{n},')');
    lastsubstr=mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1)));

if isempty(strfind(lastsubstr,'&'))==0,
    lastsubstr=strrep(lastsubstr,'&',',');
    lastsubstr=strrep(lastsubstr,'(','[');lastsubstr=strrep(lastsubstr,')',']');lastsubstr=strcat('min','(',lastsubstr,')');
    
elseif isempty(strfind(lastsubstr,'|'))==0,
    if modellogic==0,
        lastsubstr=strrep(lastsubstr,'|',',');
        lastsubstr=strrep(lastsubstr,'(','[');lastsubstr=strrep(lastsubstr,')',']');lastsubstr=strcat('max','(',lastsubstr,')');
    elseif modellogic==1,
        lastsubstr=strrep(lastsubstr,'|','+');
    end
end

     %Evaluates expression within last brackets and substitutes value in logical rule

    try
       lastsubeval{n}=eval(lastsubstr); 
       mrm{n}=strrep(mrm{n},mrm{n}(max(lastbracket):firstbracket(find(firstbracket>max(lastbracket),1))),sprintf('%d',lastsubeval{n}));
    catch
        lastsubeval{n}=0;
        mrm{n}=strrep(mrm{n},'(','');mrm{n}=strrep(mrm{n},')','');
    end
     
end

     %Perfomrs substitution when no brackets left 

if isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0)&&(isempty(strfind(mrm{n},'|'))==1),
   mrm{n}=strrep(mrm{n},'&',',');mrm{n}=strcat('min','(','[',mrm{n},']',')');   
    
elseif isempty(strfind(mrm{n},'('))==1 &&(isempty(strfind(mrm{n},'|'))==0)&&(isempty(strfind(mrm{n},'&'))==1),
    if modellogic==0,
        mrm{n}=strrep(mrm{n},'|',',');mrm{n}=strcat('max','(','[',mrm{n},']',')');
    elseif modellogic==1,
        mrm{n}=strrep(mrm{n},'|','+');
    end
elseif isempty(strfind(mrm{n},'('))==1 && (isempty(strfind(mrm{n},'&'))==0)&&(isempty(strfind(mrm{n},'|'))==0),
    mrm{n}=strrep(mrm{n},'&',',');isor=strfind(mrm{n},'|');
    if length(isor)==1,
    mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1)),'min','(','[',mrm{n}(isor(1)+1:length(mrm{n})),']',')');
    elseif length(isor)>1,
        mrm{n}=strcat('min','(','[',mrm{n}(1:isor(1)-1),']',')',mrm{n}(isor(1):length(mrm{n})));
        for q=2:length(isor),
            if q==length(isor),
               isor=strfind(mrm{n},'|');
               mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q)),'min','(','[',mrm{n}(isor(q)+1:length(mrm{n})),']',')');
            else
              isor=strfind(mrm{n},'|');
              mrm{n}=strcat(mrm{n}(1:isor(q-1)),'min','(','[',mrm{n}(isor(q-1)+1:isor(q)-1),']',')',mrm{n}(isor(q):length(mrm{n})));
            end            
        end
        
    end

end
end
  
for n=1:length(mrm),
    if isempty(mrm{n})==0,
        if modellogic==0,
        mrm{n}=strrep(mrm{n},'|',',');mrm{n}=strcat('max','(','[',mrm{n},']',')');
        elseif modellogic==1,
        mrm{n}=strrep(mrm{n},'|','+');
        end
        try
        mrm{n}=strrep(mrm{n},mrm{n},sprintf('%d',eval(mrm{n})));
        catch
            mrm{n}=modelrules{n};
        end
    end
end

    %Pinpoints mistakes in model.rules   

errormrm=[];
for t=1:length(mrm),
   errormrm(t)=isempty(strfind(mrm{t},'x'));
end
 numerrors=find(errormrm==0); 
if numerrors>0,
   disp('');
   fprintf('%g reaction rules have been found to be incomplete',length(numerrors));
   disp('');
   disp('Please, check following reaction rules: ');
   disp('');
   mrm{[numerrors]}
    
end

maprxndata=[];
for i=1:length(mrm),
    if isempty(mrm{i})==1,
        mrm{i}='NaN';
    end
   if ~isempty(str2num(mrm{i}))
   maprxndata(i)=str2num(mrm{i});
   end
end
maprxndata(isnan(maprxndata))=NaN;
maprxnvalues=maprxndata';

end
        
    
    