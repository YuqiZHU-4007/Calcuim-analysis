function plot_behavior_onset(y,y_3sd,fs,frame,frameb,trial,re_startpoint)
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.8]) ;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w') 
linewidth=1.2;scattersize=20;method='probability';

%1
set(axes,'position',[0.054,0.547,0.90,0.38]);
%CS bar
color=[0.5 0.5 0.5];
patch1=patch([[frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_end+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]'...
    [frameb.cs_start+trial.spon_bef(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*trial.test(3)]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.test(3)-trial.hab(2)+1,1)',...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
line([frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3));frameb.us_start+trial.hab(3)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.acq(3))],...
    [min(y) max(y)],'color',[1 0.5 0],'linestyle','--','linewidth',linewidth);hold on;
plot(y,'b');hold on;plot(y_3sd,'r');
xlim([1 frameb.per_cycle*trial.total]);ylim([min(y) max(y)]);%xlabel('frames');
set(gca,'FontSize',13,'XTick',[]);box on;

%subplot(3,1,2),
set(axes,'position',[0.054,0.20,0.90,0.28]);box on;
%CS bar
color=[0.5 0.5 0.5];
p1=patch([trial.hab(2)  trial.test(3) trial.test(3) trial.hab(2)],([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,...
    color,'edgecolor',color,'facealpha',0.2,'edgealpha',0.25);hold on;
%US bar
l1=line([trial.acq(2) trial.acq(3)],([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,'color','r','linewidth',linewidth);hold on;
onset=[];onset=re_startpoint;%行为
s1=scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs.ca)*fs.behavior-1)/fs.behavior,scattersize,[0 0 0],'^','filled');hold on;
xlim([1 trial.total]);ylim(([1 frame.per_cycle+1]-frame.cs_start)*fs.ca);
ylabel('Time(s)');
set(gca,'FontSize',13);box on;

%subplot(3,1,3),
set(axes,'position',[0.054,0.08,0.90,0.11]);box on;
onset=[];onset=re_startpoint;%行为
x=onset(find(onset(:,2)<frameb.cs_start) ,1);%& onset(:,1)>3
histogram(x,'facecolor','k','Normalization',method,'BinWidth',3);hold on%
set(gca,'XTick',3:3:trial.total);xlabel('Trial number');ylabel(method);
xlim([1 trial.total]);
set(gca,'FontSize',13);box on;


