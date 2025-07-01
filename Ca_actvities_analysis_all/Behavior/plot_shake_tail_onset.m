clc;clear all;close all;
[inputname,inputpath]=uigetfile('G:\data\.mat','behavior_fin');
load([inputpath,inputname]);

fishnum=inputname(1:14);
y_fin=y_3sd;y_raw_fin=delta_r_bef1;
re_startpoint_fin=re_startpoint_sd;startpoint_fin=startpoint_sd;

[inputname,inputpath]=uigetfile('G:\data\.mat','behavior_tail');
load([inputpath,inputname]);
y_tail=y_3sd;y_raw_tail=delta_r_bef1;
re_startpoint_tail=re_startpoint_sd;startpoint_tail=startpoint_sd;

[fs_ca,fs_behavior,frame,frameb,trial]=setpara;
[interfin,intertail,inter]=findintersect(re_startpoint_fin,re_startpoint_tail,62);
inter_fin=(interfin(:,1)-1)*(frameb.per_cycle)+interfin(:,2);
inter_tail=(intertail(:,1)-1)*(frameb.per_cycle)+intertail(:,2);

for k=1
hh=figure;
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.8]) ;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w') 
set(axes,'position',[0.054,0.70,0.90,0.25],'xtick',[],'xticklabel',[]);
y=y_fin;
xlim([1 frameb.per_cycle*trial.total]);ylim([min(y) max(y)]);ylabel('Fin','FontSize',16);
patch([[frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.test(3)-trial.hab(2)+1,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.2);hold on %CS on-off
line([frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3));frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3))],...
    [min(y) max(y)],'color','k','linewidth',1.05,'linestyle','--');hold on;
plot(y_raw_fin,'b');hold on;plot(y,'r');hold on
scatter(startpoint_fin,y(startpoint_fin),12,[0 0 0],'filled');
scatter(inter_fin,y(inter_fin),12,[1 0 1],'filled');
title(fishnum,'Interpreter','none','FontSize',25);
%legend('CS','US','rotation angle','detected event');

set(axes,'position',[0.054,0.42,0.90,0.25],'xtick',[1:frameb.per_cycle*3:frameb.per_cycle*trial.total]);
y=y_tail;
xlim([1 frameb.per_cycle*trial.total]);ylim([min(y) max(y)]);xlabel('Frames','FontSize',16);ylabel('Tail','FontSize',16);
patch([[frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.test(3)-trial.hab(2)+1,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.2);hold on %CS on-off
line([frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3));frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3))],...
    [min(y) max(y)],'color','k','linewidth',1.05,'linestyle','--');hold on;
plot(y_raw_tail,'b');hold on;plot(y,'r');hold on
scatter(startpoint_tail,y(startpoint_tail),12,[0 0 1]);
scatter(inter_tail,y(inter_tail),12,[1 0 1],'filled');
% title(['tail'],'FontSize',25);

set(axes,'position',[0.054,0.08,0.7,0.25]);box on;
xlabel('trial number','FontSize',16);ylabel('time(s)','FontSize',16);
x=1:trial.total;
patch([trial.hab(2)  trial.test(3) trial.test(3) trial.hab(2)],([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs_ca,'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
line([trial.acq(2) trial.acq(3)],([frame.us_start frame.us_start]-frame.cs_start)*fs_ca,'color','r');hold on;
xlim([1 trial.total]);ylim(([1 frame.per_cycle+1]-frame.cs_start)*fs_ca);
onset=[];[~,ind]=setdiff(startpoint_fin,inter_fin);onset=re_startpoint_fin(ind,:);%行为
scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior,12,[0 0 0],'filled');hold on;
onset=[];[~,ind]=setdiff(startpoint_tail,inter_tail);onset=re_startpoint_tail(ind,:);%行为
scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior,12,[0 0 1]);hold on;
onset=[];onset=interfin;%行为
scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs_ca)*fs_behavior-1)/fs_behavior,30,[1 0 1],'*');hold on;
legend('CS','US','fin nonoverlap','tail nonoverlap','overlap','Location','northeastoutside')
% set(axes,'position',[0.68,0.08,0.10,0.22]);
% histogram('Categories',{'fin','tail'},'BinCounts',[length(startpoint_fin) length(startpoint_tail)],'BarWidth',0.7);
% ylim([0 max(length(startpoint_fin),length(startpoint_tail))+0.1*max(length(startpoint_fin),length(startpoint_tail))]);
% text(categorical({'fin'}),length(startpoint_fin)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(length(startpoint_fin)));
% text(categorical({'tail'}),length(startpoint_tail)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(length(startpoint_tail)));

set(axes,'position',[0.8,0.08,0.18,0.25]);
y = [length(startpoint_fin)-length(inter),length(inter); length(startpoint_tail)-length(inter),length(inter)];
bar(categorical({'fin','tail'}),y,'stacked');ylabel('event number','FontSize',16)
ylim([0 max(length(startpoint_fin),length(startpoint_tail))+0.1*max(length(startpoint_fin),length(startpoint_tail))]);
text(categorical({'fin'}),length(startpoint_fin)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(length(startpoint_fin)),'FontSize',12);
text(categorical({'fin'}),y(1,1)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(y(1,1)),'FontSize',12);
text(categorical({'tail'}),length(startpoint_tail)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(length(startpoint_tail)),'FontSize',12);
text(categorical({'tail'}),y(2,1)+0.05*max(length(startpoint_fin),length(startpoint_tail)),num2str(y(2,1)),'FontSize',12);
legend('nonoverlap','overlap','Location','northeastoutside')

saveas(hh,[inputpath '\'  'onset_fin_and_tail' '.tif']);
savefig(hh,[inputpath '\'  'onset_fin_and_tail']);
end

%%统计所有的fin和tail结果
y_fin = y(1,:);
y_fin_all = [y_fin_all;y_fin];
y_tail = y(2,:);
y_tail_all = [y_tail_all;y_tail];
