function [area] = plot_meantrace_acq(acti,fs,trial,frame,cut_move,ind_cut_trial)

zz=acti((trial.hab(3))*frame.per_cycle+1:trial.acq(3)*frame.per_cycle);
% set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
%set(gcf,'Units','normalized', 'position',get(0,'ScreenSize')) ;
linewidth=3;legend_location='northwest';
color=[0.5 0.5 0.5];%CS bar color
ColorOrder=[0,0,0;...
    0.3,0.75,0.93;...
    0,0.45,0.74;...
    0,0,1;...
    0.93,0.69,0.13;...
    0.85,0.33,0.10;...
    1,0,0;...
     0.64,0.08,0.18;...
    0.49,0.18,0.56;...
    1,0,1;];
%% subplot(12,2,[1:2:13])
% subplot(12,1,1:8)
% edge=[100 -100];
% %CS bar
% patch1=patch(([frame.cs_start...
%     frame.cs_end...
%     frame.cs_end'...
%     frame.cs_start]-frame.cs_start)*fs.ca,...
%     [edge(1) edge(1) edge(2) edge(2)],...
%     color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
% %US line
% l1=line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[edge(1) edge(2)],'color',[1 0 0],'linestyle','--','linewidth',linewidth);hold on
% area_CS_acq_1st_block = [];
% for ii=1:trial.acq_block_num
%     ind1=trial.acq_block_trial*(ii-1)+1+trial.hab(3);
%     a=acti((ind1-1)*frame.per_cycle+1:ind1*frame.per_cycle);
%     edge(1)=min(edge(1),min(a));edge(2)=max(edge(2),max(a));
%     plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,a,'linewidth',linewidth);hold on
%     area_CS_acq_1st_block = [area_CS_acq_1st_block,sum(a(frame.cs_start:(frame.us_start-1)))];
% end
% % edge=[edge(1) edge(2)]+[-0.25 0.25]*(edge(2)-edge(1));
% edge=[edge(1) edge(2)]+[-0.06 0.06]*(edge(2)-edge(1));
% legend('CS','US','1','2','3','4','5','location',legend_location);title('1st trial of each block','FontSize',50);
% xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
% ylabel('¡÷F/F0');xlabel('Time(s)');
% set(gca,'FontSize',20);box off;
% 
% % subplot(12,2,[19:2:23])
% subplot(12,1,10:12)
% x = (1:1:length(area_CS_acq_1st_block));
% plot(x,area_CS_acq_1st_block,'-o',...
%     'LineWidth',3,...
%     'MarkerEdgeColor','r',...
%     'MarkerFaceColor','r')
% set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
% title('Area of CS-US in each block','FontSize',20);
% ylabel('Integrated ¡÷F/F0');xlabel('Blocks');

subplot(12,1,1:7,'ColorOrder',ColorOrder)
edge=[-100 100];
%CS bar
patch1=patch(([frame.cs_start...
    frame.cs_end...
    frame.cs_end'...
    frame.cs_start]-frame.cs_start)*fs.ca,...
    [edge(1) edge(1) edge(2) edge(2)],...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.5,'facealpha',0.5);hold on %CS on-off
%US line
l1=line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[edge(1) edge(2)],'color',[1 0 0],'linestyle','--','linewidth',linewidth);hold on

% area_CS_acq_block=[];
% for ii=1:trial.acq_block_num
%     ind1=trial.acq_block_trial*(ii-1)+1+trial.hab(3);
%     ind2=trial.acq_block_trial*(ii-1)+1+trial.hab(3)+trial.acq_block_trial-1;
%     %subplot(trial.acq_block_num,1,ii),
%     %plot(aa((ind1-1)*frame.per_cycle+1:ind1*frame.per_cycle));hold on
%     a=mean(reshape(acti((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)',1);
%     edge(1)=min(edge(1),min(a));edge(2)=max(edge(2),max(a));
%     plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,a,...
%         'linewidth',linewidth);hold on
%     area_CS_acq_block = [area_CS_acq_block,sum(a(frame.cs_start:(frame.us_start-1)))];
% end

[rawacti,area]=calculate_integtate_dfdf_main(acti,cut_move,trial,frame,fs,ind_cut_trial);
edge(1)=max(edge(1),min(rawacti.acq_mean_blocks{1,1}(:)));edge(2)=min(edge(2),max(rawacti.acq_mean_blocks{1,1}(:)));
if edge(1)==edge(2)
    edge=[edge(1) edge(2)]+[-0.02 0.02];
else
edge=[edge(1) edge(2)]+[-0.06 0.06]*(edge(2)-edge(1));
end
a=rawacti.acq_mean_blocks{:};
plot(([1:frame.per_cycle]-frame.cs_start)*fs.ca,a,...
    'linewidth',linewidth);hold on
% edge=[edge(1) edge(2)]+[-0.25 0.25]*(edge(2)-edge(1));
legend(['CS';'US';num2str([1:trial.acq_block_num]','%02d')],'location',legend_location);
title('Mean across trials of each block','FontSize',50);
xlim(([1,frame.per_cycle+1]-frame.cs_start)*fs.ca);ylim([edge(1) edge(2)]);
ylabel('¡÷F/F0');xlabel('Time(s)');
set(gca,'FontSize',20);box off;

subplot(12,1,10:12,'ColorOrder',ColorOrder)
x = (1:1:length(area.CS_acq_block));
edge(1)=max(-100,min(area.CS_acq_block(:)));edge(2)=min(100,max(area.CS_acq_block(:)));
if edge(1)==edge(2)
    edge=[edge(1) edge(2)]+[-0.02 0.02];
else
edge=[edge(1) edge(2)]+[-0.06 0.06]*(edge(2)-edge(1));
end
plot(x,area.CS_acq_block,'-o',...
    'LineWidth',3,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r');
ylim([edge(1) edge(2)]);
set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;
title('Area of CS-US in each block','FontSize',20);
ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
% 
% 
