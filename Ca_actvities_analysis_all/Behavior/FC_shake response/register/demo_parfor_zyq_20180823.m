clear,clc;
tic
set(0,'defaultfigurecolor','w')  
io.regPath='I:\register\imageJ.jar';
io.ijPath='I:\Fiji_64bits.app\ij.jar';
javaaddpath(io.regPath);
javaaddpath(io.ijPath);
a=imageJ.align();

%和模板(每个trial第一帧)
sourceFolder=uigetdir('Z:\','source folder');
%targetFolder=uigetdir('Z:\','target folder');
outputname=uigetdir('Z:\','output folder');
%sourceFolder='H:\20180509\fish2\roi-fin-total';
%targetFolder='H:\20180509\fish2\roi-fin-total'; 
%outputName='H:\20180509\fish2\roi_fin_reg_result\'; 
%a.alignMultiImage(source,targetFolder,outputName);
threadNum=20;%线程数
%a.alignMultiImageParallel(source,targetFolder,outputName,threadNum);
framenum=1:72000;
trialnum=60;%
totframe_percycle=1200;

%模板为每个trial第一张
parfor ii=framenum
    sourcenum=ceil(ii/totframe_percycle)*totframe_percycle;%fix((ii-1)/totframe_percycle)*totframe_percycle+1;%%当ii=n*totframe_percycle时，trial数算的不对
    source=[sourceFolder '\' num2str(sourcenum) '.tif'];
    target=[sourceFolder '\' num2str(ii) '.tif'];
    pp{ii}=a.getAlignPoints(source,target);
    %delta_r2((sourcenum-1)/totframe_percycle,ii-sourcenum+1)=2*pp{ii}(1,2);
    %delta_c2((sourcenum-1)/totframe_percycle,ii-sourcenum+1)=2*pp{ii}(1,1);
    delta_r(ii)=2*pp{ii}(1,2);delta_c(ii)=2*pp{ii}(1,1);
end
%reshape
delta_r2=reshape(delta_r,[totframe_percycle,trialnum])';
delta_c2=reshape(delta_c,[totframe_percycle,trialnum])';
toc

h1=figure;subplot(2,1,1),plot(delta_r);title('r');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_r) max(delta_r)],'color','r');hold on;
txt=num2str([1:trialnum]');
text([1:trialnum]'*totframe_percycle-800,(min(delta_r)+max(delta_r))/2*ones(trialnum,1),txt);
subplot(2,1,2),plot(delta_c);title('c');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_c) max(delta_c)],'color','r');hold on;
text([1:trialnum]'*totframe_percycle-800,(min(delta_c)+max(delta_c))/2*ones(trialnum,1),txt);
saveas(h1,[outputname  '\corr_1stframe' '.jpg']);

sourcenum=0;
%和前一帧(每个trial第一帧和自己)
for ii=framenum
    if mod(ii,totframe_percycle)==1
        sourcenum=ii;
    else
        sourcenum=ii-1;
    end
    source=[sourceFolder '\' num2str(sourcenum) '.tif'];
    target=[sourceFolder '\' num2str(ii) '.tif'];
    pp_bef1{ii}=fluctuation.getAlignPoints(source,target);
    delta_r_bef1(ii)=2*pp_bef1{ii}(1,2);
    delta_c_bef1(ii)=2*pp_bef1{ii}(1,1);
    %delta_r_bef12(ii,jj)=2*pp_bef1{ii}(1,2);
    %delta_c_bef12(ii,jj)=2*pp_bef1{ii}(1,1);
end
%reshape
delta_r_bef12=reshape(delta_r_bef1,[totframe_percycle,trialnum])';
delta_c_bef12=reshape(delta_c_bef1,[totframe_percycle,trialnum])';
toc

h2=figure;subplot(2,1,1),plot(delta_r_bef1);title('r');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_r_bef1) max(delta_r_bef1)],'color','r');hold on;
txt=num2str([1:trialnum]');
text([1:trialnum]'*totframe_percycle-800,(min(delta_r_bef1)+max(delta_r_bef1))/2*ones(trialnum,1),txt);
subplot(2,1,2),plot(delta_c_bef1);title('c');hold on;
line([[1:trialnum]'*totframe_percycle [1:trialnum]'*totframe_percycle],[min(delta_c_bef1) max(delta_c_bef1)],'color','r');hold on;
text([1:trialnum]'*totframe_percycle-800,(min(delta_c_bef1)+max(delta_c_bef1))/2*ones(trialnum,1),txt);
saveas(h2,[outputname  '\corr_before_1frame' '.jpg']);

clc;clear all;close all;
[inputname,inputpath]=uigetfile('G:\data\.mat','behavior');
load([inputpath,inputname]);
totframe_percycle=1200;
disp(['trial num:' num2str(length(delta_c_bef1)/totframe_percycle)])
disp(['hab:15']);disp(['acq:' num2str(length(delta_c_bef1)/totframe_percycle-15-6)]);disp(['test:6']);
x=delta_r_bef1;figure,plot(x);hold on;plot(y_3sd,'r')
a=input('input max thr of baseline(can be larger than it)\n');
[y_3sd,startpoint_sd,thr]=findpeak_use3sd(x,32,a,totframe_percycle);
[y_t,startpoint_t,thr]=findpeak_usetendency(x,thr,32);
% 
trials=[0 66 69];
x=delta_r_bef1;
bin=[32 32];
a=[0.5 1];
[y_3sd,startpoint_sd,thr]=findpeak_use3sd_seg(x,bin,a,trials,totframe_percycle);%分段，分成HAB,ACQ,TEST,适用于基线变化的情况

%reshape
re_startpoint_sd=[];
re_startpoint_sd(:,1)=ceil(startpoint_sd/totframe_percycle);
re_startpoint_sd(:,2)=mod(startpoint_sd,totframe_percycle);
re_startpoint_sd(re_startpoint_sd(:,2)==0,:)=[];

re_startpoint_t(:,1)=ceil(startpoint_t/totframe_percycle);
re_startpoint_t(:,2)=mod(startpoint_t,totframe_percycle);
re_startpoint_t(re_startpoint_t(:,2)==0,:)=[];

% figure,plot(y_3sd);hold on
% scatter(startpoint_sd,y_3sd(startpoint_sd));
% hold on
% plot([1:1200:72000;1:1200:72000],[-200 100]);

save([inputpath inputname],'delta_r_bef1','delta_r_bef12','delta_c_bef1','delta_c_bef12','pp_bef1',...
     'y_3sd','startpoint_sd','y_t','startpoint_t',...
     're_startpoint_sd','re_startpoint_t','-v7.3');
% save([outputName '\agin_result.mat'], 'delta_r', 'delta_r2','delta_c', 'delta_c2', ...
%     'delta_r_bef1','delta_r_bef12','delta_c_bef1','delta_c_bef12','pp','pp_bef1',...
%      'y_3sd','startpoint_sd','y_t','startpoint_t',...
%      're_startpoint_sd','re_startpoint_t','-v7.3');