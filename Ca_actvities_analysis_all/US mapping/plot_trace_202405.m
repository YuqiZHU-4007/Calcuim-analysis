function h=plot_trace_202405(data,yli,clrmap,s)
load('H:\1.Test US\5.fear conditioning behavioral data\20210709\fish2\para.mat');

% if linevisible
%     l3=plot(repmat([(frame.us_start):frame.per_cycle:frame.per_cycle*size(A,2)],2,1)',[min(yli(:)) max(yli(:))],'r','LineStyle','-','linewidth',1.2);hold on
% end
%CS bar
color=[0.5 0.5 0.5];
A=1;
patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*size(A,2)]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*size(A,2)]']',...
    repmat([min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],size(A,2),1)',...
    color,'edgecolor',color,'facecolor',color,'edgealpha',0.2,'facealpha',0.25);hold on
p=[];
for kk=1:size(data,2)
    m=squeeze(mean(data(:,kk,:),3,'omitnan'));
    shadedErrorBar([1:length(m)],squeeze(data(:,kk,:))',{@(x) mean(x,1,'omitnan'),@(x) std(x)*(1/sqrt(size(x,2)))},'lineprops',{clrmap(kk,:,:)},'transparent',1,'patchSaturation',0.1); hold on
    p(kk)=plot(m,'linewidth',1.2,'color',clrmap(kk,:,:));hold on
end
ylim(yli);
%legend([p],s,'location','northwest')