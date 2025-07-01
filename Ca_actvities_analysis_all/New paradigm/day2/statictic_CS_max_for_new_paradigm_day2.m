set(0,'defaultfigurecolor','w');
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 36;
%1:trial num;2:start trial;3:end trial
trial.spon_bef=[6*3 1 6*3];
trial.test = [6 trial.spon_bef(3)+1 trial.spon_bef(3)+6];
trial.total =trial.spon_bef(1)+trial.test(1);
fs_ca=0.4;
fs_behavior=60;

[inputname,inputpath]=uigetfile('G:\data\.xlsx','location');
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
[actname,actpath]=uigetfile('G:\data\.mat','activities_aft_process');
load([behaviorpath,behaviorname]);
load([actpath,actname]);
outputpath=uigetdir('G:\data','outputpath');


[~,sheet,~]=xlsfinfo([inputpath,inputname]);
statistic_CS_max=[];
for kk=1:size(sheet,2)
disp(sheet{kk});
outpath=checkpath([outputpath '\' sheet{kk}]);
num=xlsread([inputpath,inputname],sheet{kk});
%num=xlsread(['H:\20180623_fish1\分区\分区_20180623_fish1.xlsx'],'Olfactory');

acti=[];a_event={};a_event_rep=[];CS_max=[];
for ii=1:size(num,1)
    acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
    a_event{ii,1}=activities_event_preCS{num(ii,1),1}{num(ii,2),1}';
    a_event_rep=[a_event_rep a_event{ii,1}.ind];
end
 %figure,plot(acti(:,3))
for ii=1:size(num,1)
    aa=reshape(acti(:,ii)',frame.per_cycle,trial.total)';
    CS_max(:,ii)=max(aa(:,frame.cs_start:frame.cs_end),[],2);
    if ~isempty(find(CS_max(:,ii)>1))
        CS_max(:,ii)=ones(1,60)*100;
    end
end
%CS_max_test=mean(CS_max([trail.test(2):trial.test(3)],ii),2);
%CS_max(:,3)=[];CS_max(:,16)=[];CS_max(:,16)=[];
CS_max(:,find(CS_max(1,:)==100))=[];
%figure,plot(CS_max)
statistic_CS_max{kk,1}=CS_max;
statistic_CS_max{kk,2}=acti;
%statistic_CS_max{kk,3}=CS_max_test;
statistic_CS_max{kk,3}=sheet{kk};
end
save([outputpath '\statistic CS_max'],'statistic_CS_max');