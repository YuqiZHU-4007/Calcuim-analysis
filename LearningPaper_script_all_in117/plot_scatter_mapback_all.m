function h=plot_scatter_mapback_all(h,supervoxel_plot2,neuron_id_in_regioni,cIX,gIX,clr,seg)
if isempty(h)
    h=figure('position',[313,74,1099,719]);
end
t=tiledlayout(h,3,2);t.TileSpacing = 'compact';t.Padding = 'compact';
ax1=nexttile(t,[2 1]);    scatter3(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,1),supervoxel_plot2(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
for jj=unique(gIX)'
    ind=find(gIX==jj);
    scatter3(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),1),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),3),14,repmat(clr(jj,:),length(ind),1),'filled');hold on;
end
axis equal;view([60 40]);xlabel('X');ylabel('Y');zlabel('Z');grid on;title(seg);set(gca,'ZDir','reverse','FontWeight','bold','linewidth',2);box on;
%set(gca,'xtick',[],'ytick',[],'ztick',[]);

ax1=nexttile(t,[2 1]);    scatter(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,1),14,[0.5 0.5 0.5],'filled');hold on;
for jj=unique(gIX)'
    ind=find(gIX==jj);
    scatter(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),1),14,repmat(clr(jj,:),length(ind),1),'filled');hold on;
end
axis equal;xlabel('X');ylabel('Y');grid off;
a=get(gca,'YLim');set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'ytick',[(a(2)-a(1))/2+a(1) a(2)],'yticklabel',{'Medial','Lateral'});
set(gca,'FontWeight','bold','linewidth',2);box on;view([90 90])

ax1=nexttile(t,[1 1]);    scatter(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
for jj=unique(gIX)'
    ind=find(gIX==jj);
    scatter(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),3),14,repmat(clr(jj,:),length(ind),1),'filled');hold on;
end
axis equal;view([180 90]);xlabel('X');ylabel('Z');grid on;
set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'XDir','reverse','ytick',[get(gca,'YLim')],'yticklabel',{'Dorsal','Ventral'});
set(gca,'FontWeight','bold','linewidth',2);box on;

ax1=nexttile(t,[1 1]); 
[fraction,edges]=histcounts(gIX,'Normalization','probability','BinWidth',1,'BinEdge',[0.5:1:max(gIX)+0.5]);
hh=bar(edges(1:end-1)+0.5,fraction,'FaceColor','flat');
for k =unique(gIX)' %1:size(fraction_in_region_in_clust,2)
    hh.CData(k,:) = clr(k,:);
end
text(edges(1:end-1)+0.5,fraction,num2str(fraction'.*length(gIX)),'HorizontalAlignment','center','VerticalAlignment','bottom');
text(edges(end)+0.5,0,['N=',num2str(size(neuron_id_in_regioni,1))],'HorizontalAlignment','center','VerticalAlignment','bottom');
ylabel('Neurons(%)');set(gca,'fontsize',14,'FontWeight','bold','linewidth',2,'xtick',1:max(gIX));box off;
length(find(gIX==1));
