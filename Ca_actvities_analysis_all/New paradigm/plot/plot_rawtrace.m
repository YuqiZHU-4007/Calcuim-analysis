function plot_rawtrace(acti,X,event,fs,frame,trial,startpoint)

if nargin<3
    event=[];
end
linewidth=1.2;scattersize=20;
aa=acti;
%trial间间隔
l1=plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',linewidth);hold on
%acq block间隔
l2=plot(repmat([frame.per_cycle*trial.hab(3):frame.per_cycle*trial.acq_block_trial:frame.per_cycle*(trial.acq(3)-1)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'k--','linewidth',linewidth);hold on
%CS bar
color=[0.5 0.5 0.5];
patch1=patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
    [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
    [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
    [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']'*fs.ca,...
    repmat([min(aa(:)) min(aa(:)) max(aa(:)) max(aa(:))],trial.test(3)-trial.hab(2)+1,1)',...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
% plot(repmat([frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
%US line
l3=plot(repmat([(frame.us_start+trial.hab(3)*frame.per_cycle):frame.per_cycle:frame.per_cycle*(trial.acq(3))],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'r','LineStyle','-','linewidth',linewidth);hold on
%behavior
onset=[];onset=startpoint;
s2=scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),scattersize,[0 0 0],'^','filled');hold on;
%calcuim event
if ~isempty(event)
    s1=scatter(event.ind*fs.ca,aa(event.ind),scattersize,'k','filled');hold on;
else
    s1=[];
end
%画原始trace
if size(aa,2)==1
    p1=plot([1:frame.per_cycle*trial.total]*fs.ca,aa,'b','linewidth',1.2);hold on;
else
    shadedErrorBar([1:size(aa,1)]*fs.ca,aa',{@(x) mean(x,1),@(x) std(x,[],1)},'lineprops','r','transparent',1,'patchSaturation',0.1); hold on
    aa=mean(aa,2);
    s2=scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),scattersize,[0 0 0],'^','filled');hold on;
end
%l1=line([frame.per_cycle frame.per_cycle]*fs.ca,[min(aa(:)) max(aa(:))],'color',[0.5 0.5 0.5],'linestyle','--','linewidth',linewidth,'visible','off');hold on
legend([patch1,l3(1),s2,s1],{'CS','US','Behavior event onset','Ca event onset'},...
    'position',[0.460393411255403,0.199579498409775,0.157894733348531,0.117283947346796]);
xlim(X*fs.ca);hold on;ylim([min(aa(:)) max(aa(:))]);
set(gca,'FontSize',13);

%以前画图code
% h2=figure;
% aa=acti(:,ii);
% set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
% subplot(3,1,1);
% plot([1:frame.per_cycle*trial.total]*fs.ca,aa);hold on;xlim([trial.spon_bef(3)*frame.per_cycle+1 trial.hab(3)*frame.per_cycle]*fs.ca);hold on;ylim([min(aa(:)) max(aa(:))]);
% plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'g');hold on
% plot(repmat([frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([(trial.hab(3)*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq(3))],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'r','LineStyle','--');hold on
% scatter(event.ind*fs.ca,aa(event.ind),9,'k','filled');hold on;
% %behavior
% onset=[];onset=startpoint;
% scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),11,[1 0.5 0],'filled');hold on;
% title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Hab']);xlabel('time(s)');ylabel(lab);
%
% subplot(3,1,2);
% plot([1:frame.per_cycle*trial.total]*fs.ca,aa,'r');hold on; xlim([trial.hab(3)*frame.per_cycle+1 (trial.acq(3))*frame.per_cycle]*fs.ca);ylim([min(aa(:)) max(aa(:))]);
% plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'g');hold on
% plot(repmat([frame.per_cycle*trial.hab(3):frame.per_cycle*trial.acq_block_trial:frame.per_cycle*(trial.acq(3)-1)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'k');hold on
% %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
% plot(repmat([frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([(trial.hab(3)*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq(3))],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'b','LineStyle','--');hold on
% scatter(event.ind*fs.ca,aa(event.ind),9,'k','filled');hold on;
% %behavior
% onset=[];onset=startpoint;
% scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),11,'b','filled');hold on;
% title(['Acq']);xlabel('time(s)');ylabel(lab);
%
% subplot(3,1,3);
% plot([1:frame.per_cycle*trial.total]*fs.ca,aa);hold on;xlim([(trial.acq(3))*frame.per_cycle+1 (trial.test(3))*frame.per_cycle]*fs.ca);ylim([min(aa(:)) max(aa(:))]);
% plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'g');hold on
% plot(repmat([frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
% plot(repmat([(trial.hab(3)*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq(3))],2,1)'*fs.ca,[min(aa(:)) max(aa(:))],'r','LineStyle','--');hold on
% scatter(event.ind*fs.ca,aa(event.ind),9,'k','filled');hold on;
% %behavior
% onset=[];onset=startpoint;
% scatter(onset/fs.behavior,ones(size(onset))*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),11,[1 0.5 0],'filled');hold on;
% title(['Test']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
% saveas(h2,[outpath '\raw trace\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
% close(h2);