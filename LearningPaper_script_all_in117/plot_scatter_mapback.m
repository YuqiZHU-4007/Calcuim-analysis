function h=plot_scatter_mapback(h,supervoxel_plot2,neuron_id_in_regioni,cIX,ind,seg)
if isempty(h)
h=figure('position',[313,74,1099,719]);
end
t=tiledlayout(h,3,2);t.TileSpacing = 'compact';t.Padding = 'compact';
ax1=nexttile(t,[2 1]);
scatter3(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,1),supervoxel_plot2(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
scatter3(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),1),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),3),14,'r','filled');
axis equal;view([60 40]);xlabel('X');ylabel('Y');zlabel('Z');grid on;title(seg,'Fontsize',12,'FontWeight','bold');
set(gca,'ZDir','reverse','FontWeight','bold','linewidth',2);box on;
%set(gca,'xtick',[],'ytick',[],'ztick',[]);

ax1=nexttile(t,[2 1]);
scatter(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,1),14,[0.5 0.5 0.5],'filled');hold on;
scatter(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),1),14,'r','filled');
axis equal;;xlabel('X');ylabel('Y');grid off;
a=get(gca,'YLim');set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'ytick',[(a(2)-a(1))/2+a(1) a(2)],'yticklabel',{'Medial','Lateral'});
set(gca,'FontWeight','bold','linewidth',2);box on;view([90 90])

ax1=nexttile(t,[1 1]);
scatter(supervoxel_plot2(neuron_id_in_regioni,2),supervoxel_plot2(neuron_id_in_regioni,3),14,[0.5 0.5 0.5],'filled');hold on;
scatter(supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),2),supervoxel_plot2(neuron_id_in_regioni(cIX(ind)),3),14,'r','filled');
axis equal;view([180 90]);xlabel('X');ylabel('Z');grid on;
legend({['n=',num2str(size(neuron_id_in_regioni,1))],['n=',num2str(length(ind))]},'Location','northeast','Fontsize',10)
set(gca,'xtick',[get(gca,'XLim')],'xticklabel',{'Rostral','Caudal'},'XDir','reverse','ytick',[get(gca,'YLim')],'yticklabel',{'Dorsal','Ventral'});
set(gca,'FontWeight','bold','linewidth',2);box on;