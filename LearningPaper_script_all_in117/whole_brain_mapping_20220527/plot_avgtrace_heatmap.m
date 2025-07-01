function h=plot_avgtrace_heatmap(M_cluster,act_all_acq_block,frame,frame_ind,trial_ind_acq,startpoint,fs,frameb)

h=figure('position',[1,41,1920,962]);
a=mean(act_all_acq_block,2);
y=[-0.5,1];roundn([min(a(:)) max(a(:))],-1);
t=tiledlayout(5,1);%t.TileSpacing = 'compact';%t.Padding = 'compact';
ax1=nexttile([1,1]);plot(mean(act_all_acq_block,2),'k','linewidth',2);xlim([1 frame.per_cycle*length(trial_ind_acq)]);ylim(y);set(gca,'fontsize',16);hold on;hold on
%CS bar
patch([repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_acq)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_acq)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_acq)],1,1);...
    repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_acq)],1,1)],...
    [repmat(y(1),length(trial_ind_acq),1)';...
    repmat(y(1),length(trial_ind_acq),1)';...
    repmat(y(2),length(trial_ind_acq),1)';...
    repmat(y(2),length(trial_ind_acq),1)'],...
    'r','FaceAlpha',0.2,'EdgeColor','r','EdgeAlpha',0.2);
%US line
l3=plot(repmat([(frame.us_start):frame.per_cycle:frame.per_cycle*length(trial_ind_acq)],2,1)',[y(1) y(2)],'k','LineStyle','--','linewidth',1);hold on
%ÐÐÎª
%behavior

onset=startpoint;onset=onset-(trial_ind_acq(1)-1)*frameb.per_cycle+1;onset=onset/fs.behavior/fs.ca;
%(min(act_all_acq_block(:))+0.1*(max(act_all_acq_block(:))-min(act_all_acq_block(:))))
s2=scatter(onset,ones(size(onset))*y(1),14,[0 0 0],'^','filled');hold on;

a=mean(M_cluster,1);
ax2=nexttile([4,1]);y=[0 1];%roundn([min(a(:)) max(a(:))],-1);
imagesc(M_cluster,y);colormap('jet');colorbar;hold on;
line(repmat([1:length(frame_ind):length(frame_ind)*length(trial_ind_acq)],2,1),repmat([0 size(M_cluster,1)],length(trial_ind_acq),1)','color','y','linestyle','--','linewidth',1.5);
%linkaxes([ax1,ax2],'x');xticklabels(ax1,{});
set(gca,'fontsize',16);colorbar('Ticks',y,'Ticklabels',y);
ylabel([num2str(size(M_cluster,1)),' Neurons']);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;

end
