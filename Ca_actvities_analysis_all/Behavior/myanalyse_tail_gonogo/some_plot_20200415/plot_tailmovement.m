function plot_tailmovement(y,fs,frameb,trial,startpoint)
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.8]) ;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w') 
linewidth=1.2;scattersize=12;method='probability';
total_trial=trial.hab(1)+trial.acq(1)+trial.test(1);
re_startpoint(:,1)=fix(startpoint/frameb.per_cycle)+1;
re_startpoint(:,2)=mod(startpoint,frameb.per_cycle);

%1
set(axes,'position',[0.054,0.73,0.90,0.25]);
%CS bar
color=[0.5 0.5 0.5];
patch1=patch([[frameb.cs_start:frameb.per_cycle:frameb.per_cycle*total_trial]'...
    [frameb.cs_end:frameb.per_cycle:frameb.per_cycle*total_trial]'...
   [frameb.cs_end:frameb.per_cycle:frameb.per_cycle*total_trial]'...
     [frameb.cs_start:frameb.per_cycle:frameb.per_cycle*total_trial]']',...
    repmat([min(y) min(y) max(y) max(y)],total_trial,1)',...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
if trial.acq_block_num~=0
line([frameb.us_start+trial.hab(1)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.hab(1)+trial.acq(1));frameb.us_start+trial.hab(1)*frameb.per_cycle:frameb.per_cycle:frameb.per_cycle*(trial.hab(1)+trial.acq(1))],...
    [min(y) max(y)],'color','r','linestyle','--','linewidth',linewidth);hold on;
end
plot(y,'b');hold on;
%scatter(startpoint,y(startpoint),scattersize,[0 0 0],'o','filled');hold on;
xlim([1 frameb.per_cycle*total_trial]);ylim([min(y) max(y)]);%xlabel('frames');
set(gca,'FontSize',13,'XTick',[]);box on;

%subplot(3,1,2),
set(axes,'position',[0.25,0.08,0.4,0.6]);box on;axis ij;
clr=jet(trial.acq_block_num+2);clr_r=kron(clr,ones(trial.acq_block_trial,1));clr_r=[repmat(clr_r(1,:),trial.hab(1)-trial.acq_block_trial,1);clr_r;repmat(clr_r(end,:),trial.hab(1)-trial.acq_block_trial,1)];
re_y=reshape(y,[frameb.per_cycle,total_trial]);
[h_all, incre_all] = sepplot(([1:frameb.per_cycle]-frameb.cs_start)/fs.behavior, re_y,clr_r);
%CS bar
p1=patch(([frameb.cs_start frameb.cs_start frameb.cs_end frameb.cs_end]-frameb.cs_start)/fs.behavior,[min(incre_all(:)) max(incre_all(:)) max(incre_all(:)) min(incre_all(:))],...
    color,'edgecolor',color,'facealpha',0.2,'edgealpha',0.25);hold on;
[h_all, incre_all] = sepplot(([1:frameb.per_cycle]-frameb.cs_start)/fs.behavior, re_y,clr_r);
%US bar
if trial.acq_block_num~=0
l1=line(([frameb.us_start frameb.us_start]-frameb.cs_start)/fs.behavior,[min(incre_all([trial.hab(1)+1:trial.hab(1)+trial.acq(1)],1)) max(incre_all([trial.hab(1)+1:trial.hab(1)+trial.acq(1)],1))],...
    'LineStyle','--','color','r','linewidth',linewidth);hold on;
end
%scatter(startpoint,y(startpoint),scattersize,[0 0 0],'o','filled');hold on;
xlim(([1 frameb.per_cycle+1]-frameb.cs_start)/fs.behavior);%ylim([min(incre_all(:)) max(incre_all(:))]);
xlabel('Time(s)');
set(gca,'FontSize',13,'yticklabel',[]);box on;
%'Colorbar'
set(axes,'position',[0.655,0.08,0.01,0.6]);box on;l=0;
for zz=1:size(incre_all,1)
    ll=incre_all(zz,1)-(incre_all(2,1)-incre_all(1,1))/2;
    y = [l ll ll l];
    x = [0 0 1 1];
    l=ll;
    patch(x,y,clr_r(zz,:),'edgealpha',0);hold on; 
end
set(gca,'visible','off','yticklabel',[],'xticklabel',[]);
ylim([min(incre_all(:)) max(incre_all(:))]); xlim([0 1]);

set(axes,'position',[0.054,0.08,0.02,0.6]);box on;
x=re_startpoint(find(re_startpoint(:,2)<frameb.cs_start),1);
histogram(x,'facecolor','k','Normalization',method,'BinWidth',1,'Orientation','horizontal');hold on%
set(gca,'XTick',3:3:total_trial);set(gca,'FontSize',13);box on;
ylim([1 total_trial]);ylabel('# Trial');axis ij;


%subplot(3,1,3),
set(axes,'position',[0.078,0.08,0.16,0.6]);box on;axis ij;
%CS bar
p1=patch(([frameb.cs_start frameb.cs_start frameb.cs_end frameb.cs_end]-frameb.cs_start)/fs.behavior,[1 total_trial total_trial 1],...
    color,'edgecolor',color,'facealpha',0.2,'edgealpha',0.25);hold on;
%US bar
l1=line(([frameb.us_start frameb.us_start]-frameb.cs_start)/fs.behavior,[trial.hab(1)+1 trial.hab(1)+trial.acq(1)],'color','r','linewidth',linewidth);hold on;
onset=[];onset=re_startpoint;%ÐÐÎª
s1=scatter((onset(:,2)-(frameb.cs_start-1))/fs.behavior,onset(:,1),scattersize,[0 0 0],'s','filled');hold on;
ylim([1 total_trial]);xlim(([1 frameb.per_cycle+1]-frameb.cs_start)/fs.behavior);
xlabel('Time(s)');
set(gca,'FontSize',13,'yticklabel',[]);box on;






