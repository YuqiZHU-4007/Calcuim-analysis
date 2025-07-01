%test

%1,regress seed

%2.distance

%3.amp
batchi=1;fishi=9;
path=Path{batchi};
load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align');
load(fullfile(path{fishi}, '/para.mat'));
nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);

A=activities_dfdf_align;

frame_indCS=1:frame.cs_end;
frame_indUS=frame.us_start:frame.per_cycle;
sessions=size(A,2);
regressor1=[];regressor2=[];regressor3=[];regress1=zeros(5,sessions*frame.per_cycle); regress2=zeros(5,(sessions-2)*frame.per_cycle);regress3=zeros(5,(sessions)*frame.per_cycle);
figure,
for ii=1:sessions
    regressor1(:,ii)=zeros(1,frame.per_cycle);
    regressor1(frame.cs_start:frame.cs_end,ii)=(ii-1)/sessions*2;
    subplot(1,sessions,ii),plot(regressor1(:,ii));
    ylim([0 2])
end
figure,
for ii=1:sessions-2
    regressor2(:,ii)=zeros(1,frame.per_cycle);
    regressor2(frame.us_start:frame.us_start+2,ii)=(sessions-2-ii)/(sessions-2)*2;
    subplot(1,sessions,ii),plot(regressor2(:,ii));
    ylim([0 2])
end
figure,
for ii=1:sessions
    regressor3(:,ii)=zeros(1,frame.per_cycle);
    regressor3(frame.cs_start:frame.cs_end,ii)=(ii-1)/sessions;
    if ii>=2 && ii<=sessions-1
        regressor3(frame.us_start:frame.us_start+2,ii)=(sessions-ii+1)/(sessions)*2;
    end
    subplot(1,sessions,ii),plot(regressor3(:,ii));
    ylim([0 2])
end
regress1(1,:)=reshape(regressor1,1,[]);
regress2(1,:)=reshape(regressor2,1,[]);
regress3(1,:)=reshape(regressor3,1,[]);
close all;

[regressors,~,~,~] = GetMotorRegressor(regress2);output=[regressors(1,1).im]./max(regressors(1,1).im);
output= reshape(output,frame.per_cycle,sessions-2);output= output(frame_indUS,:);output= reshape(output,[],1);
AA=reshape(A(frame_indUS,2:7,:),length(frame_indUS)*(sessions-2),[]);
a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
[stimcorr,~] = MotorSourceCorrelation(AA', output',[]);
US_related_stimcorr{batchi,fishi}=stimcorr;

[regressors,~,~,~] = GetMotorRegressor(regress1);output=[regressors(1,1).im]./max(regressors(1,1).im);
output= reshape(output,frame.per_cycle,sessions);output= output(frame_indCS,:);output= reshape(output,[],1);
%figure,plot(output)
AA=reshape(A(frame_indCS,:,:),length(frame_indCS)*(sessions),[]);
a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
[stimcorr,~] = MotorSourceCorrelation(AA',output',[]);
CS_related_stimcorr{batchi,fishi}=stimcorr;

[regressors,~,~,~] = GetMotorRegressor(regress3);output=[regressors(1,1).im]./max(regressors(1,1).im);
output= reshape(output,frame.per_cycle,sessions);output= output(:,2:sessions-1);output= reshape(output,[],1)';
%figure,plot(output)
AA=reshape(A(:,2:7,:),frame.per_cycle*(sessions-2),[]);
a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
[stimcorr,~] = MotorSourceCorrelation(AA',output,[]);
CSUS_related_stimcorr{batchi,fishi}=stimcorr;

stimcorr=CS_related_stimcorr{batchi,fishi};
thr=0.5;%thrfun2(stimcorr);
ind_regress=find(stimcorr>thr);
% stimcorr=US_related_stimcorr{batchi,fishi};
% thr=0.2;%thrfun2(stimcorr);
% ind_regress2=find(stimcorr>thr);
% ind_regress=intersect(ind_regress,ind_regress2);

frame_indCS=frame.cs_start:frame.cs_end;
frame_indUS=frame.us_start:frame.us_start+5;
a=squeeze(nanmean(A(frame_indUS,:,:),1));
b= squeeze(nanmean(A(frame_indCS,:,:),1));
B1=(a-b)>0 & abs(a-b)>b*0.3;B2=(a-b)<0 & abs(a-b)>b*0.1;
ind_thr=find(sum(B1(2:4,:),1)>=2 & sum(B2(5:7,:),1)>=2);

ind=intersect(ind_thr,ind_regress);

CS_related_act=A(:,:,ind);CS_related_act=reshape(CS_related_act,size(CS_related_act,1)*size(CS_related_act,2),[]);
figure('position',[1,1,1200,800]),
tiledlayout(2,2);nexttile([1 2]);
patch([[frame.cs_start:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
    [frame.cs_end:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
    [frame.cs_end:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
    [frame.cs_start:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]']',...
    [-1*ones(size(activities_dfdf_align,2),1), -1*ones(size(activities_dfdf_align,2),1) 1*ones(size(activities_dfdf_align,2),1) 1*ones(size(activities_dfdf_align,2),1)]',...
    'r','facealpha',0.2,'edgealpha',0);hold on;
shadedErrorBar(1:size(CS_related_act,1),nanmean(CS_related_act,2),nanstd(CS_related_act,[],2),'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
plot(mean(CS_related_act,2),'r','linewidth',2);%legend('regressor','corr > thr');
title(num2str(length(ind)))
%figure('position',[1,1,800,400]),
nexttile,scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),0.5,[0.5 0.5 0.5]);hold on;
scatter3(supervolxeli(ind,1),supervolxeli(ind,2),supervolxeli(ind,3),16,'r','filled');grid off;view([0,-90]);
nexttile,scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),0.5,[0.5 0.5 0.5]);hold on;
scatter3(supervolxeli(ind,1),supervolxeli(ind,2),supervolxeli(ind,3),16,'r','filled');grid off;view([-90,0]);










