
[name,path]=uigetfile('H:\3.Juvenile reference brain\registration to templete\data\2021_ind\cluter_durCond\n_pc 5--n_cluster 3\dr_type 2--clut_type 1\*.txt');
%[name,path]=uigetfile(fullfile('F:\DUlab\FC_analyse\Ca_actvities_analyze\whole_brain_mapping_20220527\','-logfile.txt'));
a=textread([path,name],'%s');
ind=[];kk=1;aa=string;
for ii=1:length(a)
    if contains(a{ii},'!!!!')
        ind(kk)=ii;
        aa(kk,1)=[a{ii}];
        kk=kk+1;
    end
end