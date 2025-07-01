% a=activities_spon(:,:,2);
% A=reshape(activities_c,[],frame.per_cycle,trial.total);
% b=reshape (A(:,:,7:9),[],frame.per_cycle*3);
% aa=cat(2,b,a);bb=compute_dfdf(cat(2,aa(:,1:120),aa,aa(:,end-119:end)),241);
% 
% figure,subplot(3,1,1);plot(mean(aa,1));
% subplot(3,1,2);plot(mean(bb,1));
% cat_act_spon1_dfdf=bb(:,121:end-120);
% subplot(3,1,3);plot(mean(cat_act_spon1_dfdf,1));
% 

clc;clear all;close all;
global fs
global time
global frame
global frameb
global trial
global re_startpoint
global startpoint
global y_3sd

set(0,'defaultfigurecolor','w');
[actname,actpath]=uigetfile('J:\newfile\.mat','activities new');
[loctname,locpath]=uigetfile([actpath '.xls'],'location');
[behaviorname,behaviorpath]=uigetfile([actpath '.mat'],'behavior');
[envname,envpath]=uigetfile([actpath '.mat'],'env');
outpath=uigetdir(actpath,'outputpath');


[~,~,frame,~,trial,~,~,T_non_spon,T_spon]=setpara_spon();
load([actpath,actname]); load([envpath,envname]);
activities=activities_cutback2;
activities_r=[activities(:,1),activities,activities(:,end)];
activities_spon=zeros(size(activities_r,1),size(T_spon{1},2),trial.acq_block_num+2);
if size(activities_r,2)~=trial.total*frame.per_cycle
    warning(['wrong in addomisssion','---------act length:',num2str(size(activities_r,2))]);
end
if activities_r>trial.total*frame.per_cycle
    is_addiomission=1;
    frame_omission=(trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+1:...
        (trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+frame.per_cycle;
    activities_c=activities_r(:,setdiff([1:size(activities_r,2)],frame_omission));
else
    is_addiomission=0;activities_c=activities_r;
end

for ii=1:trial.acq_block_num+2
    activities_spon(:,:,ii)=activities_c(:,T_spon{ii});
end

activities_c=activities_c(:,T_non_spon);
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara_cut_spon([behaviorpath,behaviorname]);

c=reshape(7:30,3,[]);
cat_act_spon_dfdf=nan(88918,975,8);
for ii=2:9
a=activities_spon(:,:,ii);
A=reshape(activities_c,[],frame.per_cycle,trial.total);
b=reshape (A(:,:,c(:,ii-1)),[],frame.per_cycle*3);
aa=cat(2,b,a);bb=compute_dfdf(cat(2,aa(:,1:120),aa,aa(:,end-119:end)),241);
cat_act_spon_dfdf(:,:,ii-1)=bb(:,121:end-120);

figure,subplot(3,1,1);plot(mean(aa,1));
subplot(3,1,2);plot(mean(bb,1));
subplot(3,1,3);plot(mean(cat_act_spon_dfdf(:,:,ii-1),1));
end

%% act spon dfdf
clc;clear all;close all;

set(0,'defaultfigurecolor','w');
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        actpath=Path{batchi}{fishi};
        load( [actpath,'\act_spon.mat']);
        activities_spon_dfdf=nan(size(activities_spon,1),size(activities_spon,2),size(activities_spon,3));
        for ii=1:size(activities_spon,3)
            aa=activities_spon(:,:,ii);
            bb=compute_dfdf(cat(2,aa(:,1:120),aa,aa(:,end-119:end)),241);
            activities_spon_dfdf(:,:,ii)=bb(:,121:end-120);
        end
        save( [actpath,'\act_spon_dfdf.mat'],'activities_spon_dfdf');
    end
end
