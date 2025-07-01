function h=plot_avgtrace_heatmap2(M_cluster_hab,act_all_hab,frame,frame_ind_hab,trial_ind_hab,startpoint,fs,frameb,M_cluster_test,act_all_test,frame_ind_test,trial_ind_test)

h=figure('position',[1,41,1920,962]);
a=mean(act_all_hab,2);
y=[-0.75,3];roundn([min(a(:)) max(a(:))],-1);if y(2)<=y(1) y(2)=y(1)+0.1;end
t=tiledlayout(5,2);%t.TileSpacing = 'compact';%t.Padding = 'compact';
ax1=nexttile([1,1]);plot(mean(act_all_hab,2),'k','linewidth',2);xlim([1 frame.per_cycle*length(trial_ind_hab)]);ylim(y);set(gca,'fontsize',16);hold on;hold on
%CS bar
patch([repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_hab)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_hab)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_hab)],1,1);...
    repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_hab)],1,1)],...
    [repmat(y(1),length(trial_ind_hab),1)';...
    repmat(y(1),length(trial_ind_hab),1)';...
    repmat(y(2),length(trial_ind_hab),1)';...
    repmat(y(2),length(trial_ind_hab),1)'],...
    'r','FaceAlpha',0.2,'EdgeColor','r','EdgeAlpha',0.2);
title('Pre Cond')
%行为
%behavior
onset=startpoint;onset=onset-(trial_ind_hab(1)-1)*frameb.per_cycle+1;onset=onset/fs.behavior/fs.ca;
%(min(act_all_acq_block(:))+0.1*(max(act_all_acq_block(:))-min(act_all_acq_block(:))))
s2=scatter(onset,ones(size(onset))*y(1),14,[0 0 0],'^','filled');hold on;

%%
ax1=nexttile([1,1]);plot(mean(act_all_test,2),'k','linewidth',2);xlim([1 frame.per_cycle*length(trial_ind_test)]);ylim(y);set(gca,'fontsize',16);hold on;hold on
%CS bar
patch([repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_test)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_test)],1,1);...
    repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*length(trial_ind_test)],1,1);...
    repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind_test)],1,1)],...
    [repmat(y(1),length(trial_ind_test),1)';...
    repmat(y(1),length(trial_ind_test),1)';...
    repmat(y(2),length(trial_ind_test),1)';...
    repmat(y(2),length(trial_ind_test),1)'],...
    'r','FaceAlpha',0.2,'EdgeColor','r','EdgeAlpha',0.2);
title('Post Cond')
%行为
%behavior
onset=startpoint;onset=onset-(trial_ind_test(1)-1)*frameb.per_cycle+1;onset=onset/fs.behavior/fs.ca;
%(min(act_all_acq_block(:))+0.1*(max(act_all_acq_block(:))-min(act_all_acq_block(:))))
s2=scatter(onset,ones(size(onset))*y(1),14,[0 0 0],'^','filled');hold on;

%%
a=mean(M_cluster_test,1);
y=roundn([min(a(:)) max(a(:))],0);if y(2)<=y(1) y(2)=y(1)+0.1;end
ax2=nexttile([4,1]);imagesc(M_cluster_hab,y);colormap('jet');colorbar;hold on;
line(repmat([1:length(frame_ind_hab):length(frame_ind_hab)*length(trial_ind_hab)],2,1),repmat([0 size(M_cluster_hab,1)],length(trial_ind_hab),1)','color','y','linestyle','--','linewidth',1.5);
%linkaxes([ax1,ax2],'x');xticklabels(ax1,{});
set(gca,'fontsize',16);colorbar('Ticks',y,'Ticklabels',y);
ylabel([num2str(size(M_cluster_hab,1)),' Neurons']);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;
%%
%%
a=mean(M_cluster_test,1);
y=roundn([min(a(:)) max(a(:))],0);if y(2)<=y(1) y(2)=y(1)+0.1;end
ax2=nexttile([4,1]);imagesc(M_cluster_test,y);colormap('jet');colorbar;hold on;
line(repmat([1:length(frame_ind_test):length(frame_ind_test)*length(trial_ind_test)],2,1),repmat([0 size(M_cluster_test,1)],length(trial_ind_test),1)','color','y','linestyle','--','linewidth',1.5);
%linkaxes([ax1,ax2],'x');xticklabels(ax1,{});
set(gca,'fontsize',16);colorbar('Ticks',y,'Ticklabels',y);
ylabel([num2str(size(M_cluster_test,1)),' Neurons']);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;
end
