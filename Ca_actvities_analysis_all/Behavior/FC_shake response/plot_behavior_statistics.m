function plot_behavior_statistics(fs,frame,frameb,trial,re_startpoint)

method='probability';
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.8]) ;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w') 
%统计图
binnum=frameb.per_cycle/(0.4*fs.behavior);bin2num=4.8*fs.behavior/32;higher_than_edge=0.3;
subplot(2,4,1),
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2) & re_startpoint(:,1)<=trial.hab(3)),2);
x=(x-((frame.cs_start-1)*fs.ca)*fs.behavior-1)/fs.behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
%rotate(h2,[0 0 1],180) 
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs.ca);title('Habituation');
ylabel(method);set(gca,'FontSize',13);

subplot(2,4,2)
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.test(2) & re_startpoint(:,1)<=trial.test(3)),2);
x=(x-((frame.cs_start-1)*fs.ca)*fs.behavior-1)/fs.behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs.ca);title('Test');
set(gca,'FontSize',13);

subplot(2,4,3)
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.acq(2) & re_startpoint(:,1)<=trial.acq(3)),2);
x=(x-((frame.cs_start-1)*fs.ca)*fs.behavior-1)/fs.behavior;
h2=histogram(x,binnum,'facecolor','k','Normalization',method);hold on
patch(([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,[min(h2.Values)  max(h2.Values) max(h2.Values) min(h2.Values)],'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[min(h2.Values)  max(h2.Values)],'color','r');
xlim(([1 frame.per_cycle+1]-frame.cs_start)*fs.ca);title('Acqusition');
set(gca,'FontSize',13);

subplot(2,4,5)
plot_shakenum(re_startpoint,[trial.hab(1)-9 trial.hab(2)+10 trial.hab(3)],frameb,fs.behavior,bin2num);ylabel('mov num/sum');
subplot(2,4,6);
plot_shakenum(re_startpoint,trial.test,frameb,fs.behavior,bin2num)
subplot(2,4,7);
plot_shakenum(re_startpoint,trial.acq,frameb,fs.behavior,bin2num)

subplot(2,4,8);
onset=[];onset=re_startpoint;%行为
x=re_startpoint(find(re_startpoint(:,1)>=trial.hab(2)+10 & re_startpoint(:,1)<=trial.hab(3)),2);xx=[];
xx(1)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs.behavior;%(frameb.cs_end-frameb.cs_start);
x=re_startpoint(find(re_startpoint(:,1)>=trial.test(2) & re_startpoint(:,1)<=trial.test(3)),2);
xx(2)=length(find(x>=frameb.cs_start & x<frameb.cs_end));%/4.8*fs.behavior;%(frameb.cs_end-frameb.cs_start);
histogram('Categories',{'Bef learning','Aft learning'},'BinCounts',[xx(1),xx(2)]./[(trial.hab(1)-9)*bin2num trial.test(1)*bin2num]);
f= table([xx(1);(trial.hab(1)-9)*bin2num-xx(1)],[xx(2);trial.test(1)*bin2num-xx(2)],'VariableNames',{'Beflearning','Aftlearning'},'RowNames',{'MOV','NoMOV'});f
higher_than_edge=max([xx(1),xx(2)]./[(trial.hab(1)-9)*bin2num trial.test(1)*bin2num])*0.03;
[h,p,stats]=plot_fishertest_p(f,xx(1)/((trial.hab(1)-9)*bin2num),xx(2)/(trial.test(1)*bin2num),categorical({'Bef learning','Aft learning'}),higher_than_edge);
set(gca,'FontSize',13);
%ylabel('mov num/sum');

