function plot_fraction_regions(fraction_in_region_in_clust,Lable,cmap)

figure('position',[200,200,1300,800]),
if ~isempty(fraction_in_region_in_clust)
    c = categorical(Lable)';y=fraction_in_region_in_clust;
    b=bar(c,y,'stacked','FaceColor','flat','XDataMode','manual','xdata',[1:length(Lable)]);
    set(gca,'xtick',[],'xticklabels',Lable,'XTickLabelRotation',45,'fontsize',12);
    kk=1;
    for k =1:size(fraction_in_region_in_clust,2) %1:size(fraction_in_region_in_clust,2)
        b(k).CData = cmap(kk,:);
        kk=kk+1;
    end
    legend([b(1:size(fraction_in_region_in_clust,2))],num2str([1:size(fraction_in_region_in_clust,2)]'));
    ylim([0 1]);%
else
    disp('is empty!!')
end

