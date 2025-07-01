clc;clear all;close all;
[inputname,inputpath]=uigetfile('G:\data\.mat','behavior');
load([inputpath,inputname]);

set(0,'defaultfigurecolor','w');
[fs_ca,fs_behavior,frame,frameb,trial]=setpara;
re_startpoint=re_startpoint_sd;
startpoint=startpoint_sd;
y_3sd=y_3sd;%y_t
% frame.per_cycle = 100;
% frame.cs_start = 26;
% frame.cs_end = 38;
% frame.us_start = 34;

% totframe_percycle=1200;
% frameb.per_cycle =totframe_percycle;
% frameb.cs_start = 10*fs_behavior+1;
% frameb.cs_end = 14.8*fs_behavior+1;
% frameb.us_start =  13.2*fs_behavior+1;

% %1:trial num;2:start trial;3:end trial
% a=0;trial.spon_bef=[a min(a,1) a];
% a=15;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
% trial.acq_block_num=5;
% trial.acq_block_trial=6;
% trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
% a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
% a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
% trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);

% fs_ca=0.2;
% fs_behavior=60;

% re_startpoint=re_startpoint_sd;startpoint=startpoint_sd;
% y_3sd=y_3sd;

method='probability';

h=figure;
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.8]) ;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w') 
%subplot(3,1,1),
title(inputname,'Interpreter','none')
set(axes,'position',[0.054,0.547,0.90,0.38]);
y=delta_r_bef1;
xlim([1 frameb.per_cycle*trial.total]);ylim([min(y) max(y)]);xlabel('frames');
patch([[frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.test(3)-trial.hab(2)+1,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.2);hold on %CS on-off
line([frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3));frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3))],...
    [min(y) max(y)],'color','k','linewidth',1.05,'linestyle','--');hold on;
plot(y,'b');hold on;plot(y_3sd,'r');

%subplot(3,1,2),
set(axes,'position',[0.054,0.20,0.90,0.28]);box on;
xlabel('trial number');ylabel('time(s)');
x=1:trial.total;
patch([trial.hab(2)  trial.test(3) trial.test(3) trial.hab(2)],([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs_ca,'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
line([trial.acq(2) trial.acq(3)],([frame.us_start frame.us_start]-frame.cs_start)*fs_ca,'color','r');hold on;
xlim([1 trial.total]);ylim(([1 frame.per_cycle+1]-frame.cs_start)*fs_ca);
onset=[];onset=re_startpoint;%行为
scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior,12,[0 0 0],'filled');hold on;

%subplot(3,1,3),
set(axes,'position',[0.054,0.08,0.90,0.11]);box on;
onset=[];onset=re_startpoint;%行为
x=onset(find(onset(:,2)<frameb.cs_start) ,1);%& onset(:,1)>3
% [C,ia,ic] = unique(x);
% a_counts = accumarray(ic,1);
% value_counts = [C, a_counts]
% histogram('Categories',{'Bef CS(4.8s)','CS','Aft CS(4.8s)'},'BinCounts',[xx(1),xx(2),xx(3)]/sum(xx));
histogram(x,'facecolor','k','Normalization',method,'BinWidth',3);hold on%
set(gca,'XTick',3:3:trial.total);xlabel('trial number');xlim([1 trial.total]);

saveas(h,[inputpath '\'  'behavior onset' '.tif']);
savefig(h,[inputpath '\'  'behavior onset']);

%统计图
hh=figure;
binnum=frameb.per_cycle/(0.4*fs_behavior);bin2num=4.8*fs_behavior/32;higher_than_edge=0.3;
subplot(2,4,1),
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2) & re_startpoint(:,1)<=trial.hab(3)),2);
x=(x-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs_ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
%rotate(h2,[0 0 1],180) 
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs_ca);title('Habituation');
ylabel(method);

subplot(2,4,2)
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.test(2) & re_startpoint(:,1)<=trial.test(3)),2);
x=(x-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs_ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs_ca);title('Test');

subplot(2,4,3)
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.acq(2) & re_startpoint(:,1)<=trial.acq(3)),2);
x=(x-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs_ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
line(([frame.us_start frame.us_start]-frame.cs_start)*fs_ca,[min(h2.Values)  max(h2.Values)],'color','r');
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs_ca);title('Acqusition');

subplot(2,4,5)
plot_shakenum(re_startpoint,[trial.hab(1)-9 trial.hab(2)+10 trial.hab(3)],frameb,fs_behavior,bin2num);ylabel('mov num/sum');
subplot(2,4,6);
plot_shakenum(re_startpoint,trial.test,frameb,fs_behavior,bin2num)
subplot(2,4,7);
plot_shakenum(re_startpoint,trial.acq,frameb,fs_behavior,bin2num)

subplot(2,4,8);
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2)+10 & re_startpoint(:,1)<=trial.hab(3)),2);xx=[];
xx(1)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs_behavior;%(frameb.cs_end-frameb.cs_start);
x=re_startpoint(find(re_startpoint(:,1)>=trial.test(2) & re_startpoint(:,1)<=trial.test(3)),2);
xx(2)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs_behavior;%(frameb.cs_end-frameb.cs_start);
histogram('Categories',{'Bef learning','Aft learning'},'BinCounts',[xx(1),xx(2)]./[(trial.hab(1)-9)*bin2num trial.test(1)*bin2num]);
f= table([xx(1);(trial.hab(1)-9)*bin2num-xx(1)],[xx(2);trial.test(1)*bin2num-xx(2)],'VariableNames',{'Beflearning','Aftlearning'},'RowNames',{'MOV','NoMOV'});f
higher_than_edge=max([xx(1),xx(2)]./[(trial.hab(1)-9)*bin2num trial.test(1)*bin2num])*0.03;
[h,p,stats]=plot_fishertest_p(f,xx(1)/((trial.hab(1)-9)*bin2num),xx(2)/(trial.test(1)*bin2num),categorical({'Bef learning','Aft learning'}),higher_than_edge);
%ylabel('mov num/sum');

saveas(hh,[inputpath '\'  'behavior statistics' '.tif']);
savefig(hh,[inputpath '\'  'behavior statistics']);