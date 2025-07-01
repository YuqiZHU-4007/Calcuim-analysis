function [rawacti, area] = plot_mean_hab_test_control(acti,fs,trial,frame,frameb,startpoint,cut_move,ind_cut_trial)

linewidth=3;legend_location='northwest';scattersize=50;
color=[0.5 0.5 0.5];%CS bar color
set(gcf,'Units','normalized', 'position',[0.05,0.05,1.2,0.65]) ;
% set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.65]) ;
% set(gcf,'Units','normalized', 'position',get(0,'ScreenSize')) ;
if trial.hab(1)<=trial.test(1)
    indcut1=0;
    indcut2=trial.test(1)-trial.hab(1);
else
    indcut1=trial.hab(1)-trial.test(1);
    indcut2=0;
end
zz=acti([(trial.hab(2)-1)*frame.per_cycle+1:trial.hab(3)*frame.per_cycle trial.acq(3)*frame.per_cycle+1:trial.test(3)*frame.per_cycle]);


subplot(2,4,[1,2,5,6]),
rawacti=[];
edge=[100 -100];
rawacti(:,1)=mean(reshape(acti((trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3)-indcut1)*frame.per_cycle),frame.per_cycle,[]),2);%hab
rawacti(:,2)=mean(reshape(acti((trial.test(2)-1)*frame.per_cycle+1:(trial.test(3)-indcut2)*frame.per_cycle),frame.per_cycle,[]),2);%test
if max(rawacti(frame.cs_start:frame.cs_end,1))>max(rawacti(frame.cs_start:frame.cs_end,2))
    ttest_tail='right';
else
    ttest_tail='left';
end
[h,p]=ttest2(rawacti(frame.cs_start:frame.cs_end,1),rawacti(frame.cs_start:frame.cs_end,2),'Tail',ttest_tail);
aa=acti((trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3)-indcut1)*frame.per_cycle);
aa=[aa acti((trial.test(2)-1)*frame.per_cycle+1:(trial.test(3)-indcut2)*frame.per_cycle)];
edge=[min(aa(:)) max(aa(:))];edge=[edge(1) edge(2)]+[-0.6 1.1]*min((edge(2)-edge(1)),0.05);

% edge=[min(zz) max(zz)];edge=[edge(1) edge(2)]+[-0.25 0.25]*(edge(2)-edge(1));
% edge=[min(zz) max(zz)];edge=[edge(1) edge(2)]+[-0.08 0.08]*(edge(2)-edge(1));
%CS bar
p1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2)]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2)]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2)]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2)]']'*fs.ca,...
    repmat([edge(1) edge(1) edge(2) edge(2)],(trial.test(1)-indcut2),1)',...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on 
%hab
p2=plot(([1:(trial.test(1)-indcut2)*frame.per_cycle])*fs.ca,...
    acti((trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3)-indcut1)*frame.per_cycle),'b','linewidth',linewidth); hold on;
%test
p3=plot(([1:(trial.test(1)-indcut2)*frame.per_cycle])*fs.ca,...
    acti((trial.test(2)-1)*frame.per_cycle+1:(trial.test(3)-indcut2)*frame.per_cycle),'r','linewidth',linewidth); hold on;
%trial间隔
p4=line([frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2);frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.test(1)-indcut2)]*fs.ca,...
    [edge(1) edge(2)],'color',[0.5 0.5 0.5],'linestyle','--');hold on;
%hab行为
onset=[];onset=startpoint-(trial.hab(2)-1)*frameb.per_cycle;onset=onset(find(onset<=(trial.hab(3)-indcut1)*frameb.per_cycle));onset=onset(find(onset>0));
s2=scatter(onset/fs.behavior,ones(size(onset))*(min(edge)+0.1*(max(edge)-min(edge))),scattersize,'b','^','filled');hold on;
%test行为
onset=[];onset=startpoint-(trial.test(2)-1)*frameb.per_cycle;onset=onset(find(onset<=(trial.test(3)-indcut2)*frameb.per_cycle));onset=onset(find(onset>0));
s1=scatter(onset/fs.behavior,ones(size(onset))*(min(edge)+0.06*(max(edge)-min(edge))),scattersize,'r','^','filled');hold on;
xlim([1 (trial.test(1)-indcut2)*frame.per_cycle]*fs.ca); ylim([min(edge) max(edge)]);
ylabel('△F/F0');xlabel('Time(s)');legend('off')
set(gca,'FontSize',20);box off;

subplot(2,4,[3,7]),
% rawacti=[];
% edge=[100 -100];
% rawacti(:,1)=mean(reshape(acti((trial.hab(2)+1)*frame.per_cycle+1:(trial.hab(3)-indcut+2)*frame.per_cycle),frame.per_cycle,[]),2);%hab
% rawacti(:,2)=mean(reshape(acti((trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle),frame.per_cycle,[]),2);%test
% if max(rawacti(frame.cs_start:frame.cs_end,1))>max(rawacti(frame.cs_start:frame.cs_end,2))
%     ttest_tail='right';
% else
%     ttest_tail='left';
% end
% [h,p]=ttest2(rawacti(frame.cs_start:frame.cs_end,1),rawacti(frame.cs_start:frame.cs_end,2),'Tail',ttest_tail);
%CS bar
edge=[min(rawacti(:)) max(rawacti(:))];edge=[edge(1) edge(2)]+[-0.6 1.1]*(edge(2)-edge(1));
patch1=patch(([frame.cs_start...
    frame.cs_end...
    frame.cs_end'...
    frame.cs_start]-frame.cs_start)*fs.ca,...
    [edge(1) edge(1) edge(2) edge(2)],...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on 
%mean_hab_test
p2=plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,rawacti(:,1),'b','linewidth',linewidth); hold on;
p3=plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,rawacti(:,2),'r','linewidth',linewidth);hold on;
%显著性
lineheight=max(max(rawacti(frame.cs_start:frame.cs_end,:)));
l1=line(([frame.cs_start frame.cs_end]-frame.cs_start)*fs.ca,[lineheight lineheight]+0.1*(max(edge)-min(edge)),...
        'color','k','linewidth',2);hold on;
text(((frame.cs_start+frame.cs_end)/2-frame.cs_start)*fs.ca,...
     lineheight+0.125*(max(edge)-min(edge)),...
    num2str(p),'fontsize',14,'HorizontalAlignment','center');hold on;
title('Mean','FontSize',50);xlabel('Time(s)');

xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca); ylim([edge(1) edge(2)]);
legend([patch1,p2,p3,s2,s1],{'cs','Habituation','Test','Hab. shake onset','Test shake onset'},'position',[[0.78 0.7 0.12 0.225]]);
set(gca,'FontSize',20);box off;

subplot(2,4,8)
[~,area]=calculate_integtate_dfdf_main_control(acti,cut_move,trial,frame,fs,ind_cut_trial);
%area_CS_hab_tst = sum(rawacti(frame.cs_start:(frame.cs_end-1),:));%trapz(frame.cs_start:(frame.cs_end-1),rawacti(frame.cs_start:(frame.cs_end-1),:));
x = [1 2];
plot(x,area.CS_hab_tst,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r')
set(gca,'FontSize',20,'xtick',x,'xticklabel',{'Hab';'Test'},'TickLength',[0.01 0.01]);box off;
y_min = min(area.CS_hab_tst)-0.2*abs(min(area.CS_hab_tst));
y_max = max(area.CS_hab_tst)+0.2*abs(min(area.CS_hab_tst));
if y_max==0 && y_min==0
    ylim([-0.2 0.2]);
else
ylim([y_min y_max]);
end
xlim([0.6 2.4]);
title('Area of CS','FontSize',20);
ylabel('Integrated △F/F0');



