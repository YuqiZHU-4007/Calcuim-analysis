set(0,'defaultfigurecolor','w');
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 34;

trial.total = 60;
trial.hab = 15;
trial.acq = 30;
trial.test = 15;
fs_ca=0.4;
fs_behavior=60;

[inputname,inputpath]=uigetfile('G:\data\.xlsx','location');
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
[actname,actpath]=uigetfile('G:\data\.mat','activities_aft_process');
load([behaviorpath,behaviorname]);
load([actpath,actname]);
outputpath=uigetdir('G:\data','outputpath');

[~,sheet,~]=xlsfinfo([inputpath,inputname]);

statistic_CS_max=[];cs_res_mean=[];

for kk=1:size(sheet,2)
disp(sheet{kk});
outpath=checkpath([outputpath '\' sheet{kk}]);
num=xlsread([inputpath,inputname],sheet{kk});
%num=xlsread(['H:\20180623_fish1\分区\分区_20180623_fish1.xlsx'],'Olfactory');

acti=[];a_event={};a_event_rep=[];CS_max=[];acti_durcs=[];
for ii=1:size(num,1)
    acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
   % figure,plot(activities_new{num(ii,1),1}(51:100,num(ii,2)))
end
for ii=1:size(num,1)
    aa=reshape(acti(:,ii)',frame.per_cycle,trial.total)';
    CS_max(:,ii)=max(aa(:,frame.cs_start:frame.cs_end),[],2);
    CS_max(trial.hab+1:trial.hab+trial.acq,ii)=max(aa(trial.hab+1:trial.hab+trial.acq,frame.cs_start:frame.us_start-1),[],2);
   
    %aaa=normalize(acti(:,ii),1); 
    aaa=acti(:,ii); 
    aaaa=reshape(aaa,frame.per_cycle,trial.total)';
    acti_durcs(:,ii)=reshape(aaaa(trial.hab+1:trial.hab+trial.acq,frame.cs_start:frame.us_start-1)',[],1);%trial.hab+1:trial.hab+trial.acq,
    %figure,plot(acti(:,ii));
    if ~isempty(find(CS_max(:,ii)>1))
        CS_max(:,ii)=ones(1,60)*100;
    end
end
CS_max(:,find(CS_max(1,:)==100))=[];
%figure,plot(CS_max)
statistic_CS_max{kk,1}=CS_max;
statistic_CS_max{kk,2}=mean(CS_max,2);
statistic_CS_max{kk,3}=acti;
statistic_CS_max{kk,4}=sheet{kk};
cs_res_mean(:,kk)=mean(acti_durcs,2);
% CS_max_hab_test=CS_max([1:trial.hab trial.hab+trial.acq+1:trial.total],ii);
% CS_max_acq=CS_max(trial.hab+1:trial.hab+trial.acq,ii);
end
save([outputpath '\statistic CS_max'],'statistic_CS_max');
figure,plot(cs_res_mean);

data = cs_res_mean;
[R,P] = corrcoef(data);
P = 1 - P;
clims = [-1 1];
figure,f = imagesc(R,clims);
set(gca,'YTick',statistic_CS_max{4});