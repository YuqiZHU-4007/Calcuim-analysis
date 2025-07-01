%onset time
clc;clear all;
isspon2=true; istraining=true;

input=struct;
[actname,actpath]=uigetfile('E:\A_Data_lightsheet\Data_huc\20190905\*.mat','act');
[envname,envpath]=uigetfile([actpath '\*.mat'],'env');
if isspon2
    [actname2,actpath2]=uigetfile([actpath '\*.mat'],'act2');
    a2=load([actpath2,actname2]);
end
if istraining
    [actname3,actpath3]=uigetfile([actpath '\*.mat'],'act_training_aft_process');
    a3=load([actpath3,actname3]);
end
load([envpath,envname]);
a1=load([actpath,actname]);
input.spacial_location=env.supervoxel(:,1:3);

winsize=121;
activities=a1.activities;
dfdf=compute_dfdf(activities,winsize);
dfdf_m=dfdf;%smoothdata(dfdf,2,'movmean',3);
figure;
n=randi(size(activities,1),[1,1]);
subplot(3,1,1),plot(activities(n,:)');xlim([1 size(activities,2)]);
subplot(3,1,2),plot(dfdf(n,:)');xlim([1 size(activities,2)]);
subplot(3,1,3),plot(dfdf_m(n,:)');xlim([1 size(activities,2)]);
input.act=dfdf(:,(winsize-1)/2+1:end-(winsize-1)/2);

if isspon2
    activities=a2.activities;
    dfdf=compute_dfdf(activities,winsize);
    dfdf_m=dfdf;%smoothdata(dfdf,2,'movmean',3);
    figure;
    n=randi(size(activities,1),[1,1]);
    subplot(3,1,1),plot(activities(n,:)');xlim([1 size(activities,2)]);
    subplot(3,1,2),plot(dfdf(n,:)');xlim([1 size(activities,2)]);
    subplot(3,1,3),plot(dfdf_m(n,:)');xlim([1 size(activities,2)]);
    input.act2=dfdf(:,(winsize-1)/2+1:end-(winsize-1)/2);
end
if istraining
    activities=[a3.activities_preCS_dfdf_aftcorrect]';
    %dfdf=compute_dfdf(activities,winsize);
    figure;
    n=randi(size(activities,1),[1,1]);
    plot(activities(n,:)');xlim([1 size(activities,2)]);
    input.act_training=activities;
end
clear dfdf;clear dfdf_m;
%% spon
A=input.act;
figure,imagesc(A',[0 0.05]);colorbar;
A=normalize(A,2,'zscore');
[cof,pca1,latent,tsquare]  = pca(A');
figure,plot(cumsum(latent)./sum(latent))

A=input.act2;
figure,imagesc(A',[0 0.05]);colorbar;
A=normalize(A,2,'zscore');
[cof,pca2,latent,tsquare]  = pca(A');
figure,plot(cumsum(latent)./sum(latent))

 figure,plot3(pca1(:,1),pca1(:,2),pca1(:,3),'r');hold on;
 plot3(pca2(:,1),pca2(:,2),pca2(:,3),'b');hold on;
 xlabel('pc1'); ylabel('pc2'); zlabel('pc3');
 figure,plot(pca1(:,1))

%% training
A=reshape(input.act_training',frame.per_cycle,trial.total,[]);figure;plot(A(:,30,n))
%A=reshape(input.act_training',frame.per_cycle,trial.total,[]);
A=A(:,trial.acq(2):trial.acq(3),:);
M=[];label_roi=[];label_trial=[];
for ii=1:size(A,2)
    M=cat(3,M,A(:,ii,:));
    label_roi=cat(2,label_roi,[1:size(A,3)]);
    label_trial=cat(2,label_trial,ii*ones(1,size(A,3)));
end
M=reshape(M,frame.per_cycle,[])';figure,plot(M(find(label_roi==n & label_trial==24),:));

base=M(:,frame.cs_start-15:frame.cs_start-1); m=mean(base(:));sd=std(base(:));
onset=findonset(M(:,frame.cs_start:frame.us_start)',m+3*sd);figure,plot(onset)

figure('name','M'),imagesc(M,[0 0.05]);colorbar;hold on;
scatter(onset+(frame.cs_start-1),[length(onset):-1:1]);


