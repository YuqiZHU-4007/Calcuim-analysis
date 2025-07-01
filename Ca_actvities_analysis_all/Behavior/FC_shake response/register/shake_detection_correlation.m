clear all; close all; clc;
%%
tic
io.regPath='Z:\data\fear_conditioning\matlab_analysis\behavior\shake_detect\correlation\register\imageJ.jar';
io.ijPath='D:\Fiji.app\ij.jar';
javaaddpath(io.regPath);
javaaddpath(io.ijPath); 
fluctuation=imageJ.align();

frame.per_cycle = 1200;
frame.cs_start = 601;
frame.cs_end = 889;
frame.us_start = 793;
frame.cs_dur = frame.cs_end-frame.cs_start;
frame.fps = 60;
frame.bin = 32;
frame.csus_interval = frame.us_start-frame.cs_start;

trial.hab = 6;
trial.acq = 30;
trial.test = 6;
trial.total = trial.hab + trial.acq + trial.test;
trial.perblock = 6;
trial.block = fix(trial.acq/trial.perblock);
framenum = trial.total*frame.per_cycle;

%ï¿½ï¿½Ä£ï¿½ï¿½(Ã¿ï¿½ï¿½trialï¿½ï¿½Ò»Ö¡)
sourceFolder=uigetdir('Z:\','source folder');
%targetFolder=uigetdir('Z:\','target folder');
outputName=uigetdir('Z:\','output folder');
%sourceFolder='H:\20180509\fish2\roi-fin-total';
%targetFolder='H:\20180509\fish2\roi-fin-total'; 
%outputName='H:\20180509\fish2\roi_fin_reg_result\'; 
%a.alignMultiImage(source,targetFolder,outputName);
threadNum=20;%ï¿½ß³ï¿½ï¿½ï¿½
%a.alignMultiImageParallel(source,targetFolder,outputName,threadNum);

%Ä£ï¿½ï¿½ÎªÃ¿ï¿½ï¿½trialï¿½ï¿½Ò»ï¿½ï¿½
% parfor ii=framenum
%     sourcenum=fix((ii-1)/frame.per_cycle)*frame.per_cycle+1;%%ï¿½ï¿½ii=n*frame.per_cycleÊ±ï¿½ï¿½trialï¿½ï¿½ï¿½ï¿½Ä²ï¿½ï¿½ï¿?
%     source=[sourceFolder '\' num2str(sourcenum) '.tif'];
%     target=[sourceFolder '\' num2str(ii) '.tif'];
%     pp{ii}=a.getAlignPoints(source,target);
%     %delta_r2((sourcenum-1)/frame.per_cycle,ii-sourcenum+1)=2*pp{ii}(1,2);
%     %delta_c2((sourcenum-1)/frame.per_cycle,ii-sourcenum+1)=2*pp{ii}(1,1);
%     delta_r(ii)=2*pp{ii}(1,2);delta_c(ii)=2*pp{ii}(1,1);
% end
% %reshape
% delta_r2=reshape(delta_r,[frame.per_cycle,trial.total])';
% delta_c2=reshape(delta_c,[frame.per_cycle,trial.total])';
% toc
% 
% h1=figure;subplot(2,1,1),plot(delta_r);title('r');hold on;
% line([[1:trial.total]'*frame.per_cycle [1:trial.total]'*frame.per_cycle],[min(delta_r) max(delta_r)],'color','r');hold on;
% txt=num2str([1:trial.total]');
% text([1:trial.total]'*frame.per_cycle-800,(min(delta_r)+max(delta_r))/2*ones(trial.total,1),txt);
% subplot(2,1,2),plot(delta_c);title('c');hold on;
% line([[1:trial.total]'*frame.per_cycle [1:trial.total]'*frame.per_cycle],[min(delta_c) max(delta_c)],'color','r');hold on;
% text([1:trial.total]'*frame.per_cycle-800,(min(delta_c)+max(delta_c))/2*ones(trial.total,1),txt);
% saveas(h1,[outputName  '\corr_1stframe' '.jpg']);

sourcenum=0;
%ï¿½ï¿½Ç°Ò»Ö¡(Ã¿ï¿½ï¿½trialï¿½ï¿½Ò»Ö¡ï¿½ï¿½ï¿½Ô¼ï¿½)
for ii=1:(trial.total*frame.per_cycle)
    if mod(ii,frame.per_cycle)==1
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
delta_r_bef12=reshape(delta_r_bef1,[frame.per_cycle,trial.total])';
delta_c_bef12=reshape(delta_c_bef1,[frame.per_cycle,trial.total])';
toc

h2=figure;subplot(2,1,1),plot(delta_r_bef1);title('r');hold on;
line([[1:trial.total]'*frame.per_cycle [1:trial.total]'*frame.per_cycle],[min(delta_r_bef1) max(delta_r_bef1)],'color','r');hold on;
txt=num2str([1:trial.total]');
text([1:trial.total]'*frame.per_cycle-800,(min(delta_r_bef1)+max(delta_r_bef1))/2*ones(trial.total,1),txt);
subplot(2,1,2),plot(delta_c_bef1);title('c');hold on;
line([[1:trial.total]'*frame.per_cycle [1:trial.total]'*frame.per_cycle],[min(delta_c_bef1) max(delta_c_bef1)],'color','r');hold on;
text([1:trial.total]'*frame.per_cycle-800,(min(delta_c_bef1)+max(delta_c_bef1))/2*ones(trial.total,1),txt);
% saveas(h2,[outputName  '\corr_before_1frame' '.jpg']);

startpoint_sd=[];
startpoint_t=[];
re_startpoint_sd=[];
re_startpoint_t=[];
thr=[];

% %ï¿½ï¿½spontaneousï¿½ï¿½ï¿½ï¿½Òª
% delta_r_bef1 = delta_r_bef1(1,21601:75600);

x=delta_r_bef1;figure,plot(x);
    hold on
    %ï¿½ï¿½Ã¿ï¿½ï¿½trialï¿½Ö¿ï¿½
    x1 = [frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total];
    x1 = repmat(x1,2,1)';
    y1 = [min(x(:)) max(x(:))];
    plot(x1,y1,'g');
    %ï¿½ï¿½ï¿½Ã¿ï¿½ï¿½acqï¿½ï¿½ï¿½ï¿½Ê¼ï¿½Í½ï¿½ï¿½ï¿½
    x2 = [frame.per_cycle*trial.hab:frame.per_cycle*trial.perblock:frame.per_cycle*(trial.acq+trial.hab-1)];
    x2 = repmat(x2,2,1)';
    y2 = [min(x(:)) max(x(:))];
    plot(x2,y2,'k','LineWidth',2);
    x3 = [frame.per_cycle*(trial.hab+1):frame.per_cycle*trial.perblock:frame.per_cycle*(trial.acq+trial.hab)];
    x3 = repmat(x3,2,1)';
    y3 = [min(x(:)) max(x(:))];
    plot(x3,y3,'k','LineWidth',2);
    %ï¿½ï¿½ï¿½CSï¿½ï¿½US
    x4 = [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total];
    x4 = repmat(x4,2,1)';
    y4 = [min(x(:)) max(x(:))];
    plot(x4,y4,'b','LineStyle','--','LineWidth',1);
    x5 = [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total];
    x5 = repmat(x5,2,1)';
    y5 = [min(x(:)) max(x(:))];
    plot(x5,y5,'b','LineStyle','--','LineWidth',1);
    x6 = [(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)];
    x6 = repmat(x6,2,1)';
    y6 = [min(x(:)) max(x(:))];
    plot(x6,y6,'r','LineStyle','--','LineWidth',1);
    %ï¿½ï¿½ï¿½ï¿½ï¿½trial
    xt = [20:frame.per_cycle:frame.per_cycle*trial.total];
    yt = -1*ones(size(xt));
    txt = num2cell(1:trial.total);%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ÒªÊ¹ï¿½ï¿½cellï¿½ï¿½Ê½ï¿½ï¿½ï¿½Ð£ï¿½ï¿½ï¿½ï¿½Ê¹ï¿½ï¿½chrï¿½ï¿½Ê½ï¿½Ä»ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Òª×ªï¿½Ã³ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½
    text(xt,yt,txt)
    hold off
fluctuation=input('input max thr of baseline(can be larger than it)\n');
[y_3sd,startpoint_sd,endpoint_sd,thr_sd]=findpeak_use3sd(x,frame.bin,fluctuation,frame.per_cycle);
[y_t,startpoint_t,thr_t]=findpeak_usetendency(x);

% trials=[0 trial.hab (trial.hab+trial.acq) (trial.total)];
% x=delta_r_bef1;
% bin=[30 30 30];
% a=[20 2 5];
% [y_3sd,startpoint_sd,thr_sd]=findpeak_use3sd_seg(x,bin,a,trials,frame.per_cycle);%ï¿½Ö¶Î£ï¿½ï¿½Ö³ï¿½HAB,ACQ,TEST,ï¿½ï¿½ï¿½ï¿½ï¿½Ú»ï¿½ï¿½ß±ä»¯ï¿½ï¿½ï¿½ï¿½ï¿?

% trial=[0 33 57 81]; %after201810
% x=delta_r_bef1;
% bin=[32 32 32];
% fluctuation=[2 2 2];
% [y,startpoint_sd,thr]=findpeak_use3sd_seg(x,[32 32 32],fluctuation,trials,frame.per_cycle);%ï¿½Ö¶Î£ï¿½ï¿½Ö³ï¿½HAB,ACQ,TEST,ï¿½ï¿½ï¿½ï¿½ï¿½Ú»ï¿½ï¿½ß±ä»¯ï¿½ï¿½ï¿½ï¿½ï¿?
% [y_t,startpoint_t,thr]=findpeak_usetendency(x);

%ÕÒµ½ÏÂÒ»´ÎÎ²°Í¶¯µÄÆðÊ¼µãÓëÉÏ´ÎÎ²°Í¶¯µÄ½áÊøµã£¬Èç¹ûÐ¡ÓÚ3¸öµã£¬¾ÍËãÊÇÒ»´ÎÎ²°Í¶¯
nextStart_minus_preEnd = [startpoint_sd,0] - [0,endpoint_sd];
[locs_start_end] = find(nextStart_minus_preEnd(2:end-1)<4);
startpoint_sd(locs_start_end+1) = [];
endpoint_sd(locs_start_end) = [];
%reshapeÆðÊ¼µãºÍ½áÊøµã
re_start_end_point(:,1) = fix(startpoint_sd/frame.per_cycle)+1;
re_start_end_point(:,2) = mod(startpoint_sd,frame.per_cycle);
%     re_start_end_point(re_start_end_point(:,2)==0,:) = [];
re_endpoint_sd(:,1) = fix(endpoint_sd/frame.per_cycle)+1;
re_endpoint_sd(:,2) = mod(endpoint_sd,frame.per_cycle);
%     re_endpoint_sd(re_endpoint_sd(:,2)==0,:) = [];
re_start_end_point(:,3) = re_endpoint_sd(:,2);
[x_zero,y_zero] = find(re_start_end_point==0);
re_start_end_point(x_zero,:) = [];


% re_startpoint_t(:,1)=fix(startpoint_t/frame.per_cycle)+1;
% re_startpoint_t(:,2)=mod(startpoint_t,frame.per_cycle);
% re_startpoint_t(re_startpoint_t(:,2)==0,:)=[];


% save([outputName '\agin_result.mat'], 'delta_r', 'delta_r2','delta_c', 'delta_c2', ...
%     'delta_r_bef1','delta_r_bef12','delta_c_bef1','delta_c_bef12','pp','pp_bef1',...
%      'y_3sd','startpoint_sd','y_t','startpoint_t',...
%      're_startpoint_sd','re_startpoint_t','-v7.3');
  save([outputName '\shake.mat'],'trial','frame','delta_r_bef1','delta_r_bef12','delta_c_bef1','delta_c_bef12','pp_bef1',...
     'y_3sd','startpoint_sd','endpoint_sd','re_start_end_point','fluctuation');