    %count
    kk=1;hh=[];c={};
    for ii=unique(index_all(cIX_iii))';
        a=find(index_all(cIX_iii)==ii);
        h=histogram(gIX_iii(a),'binedge',[1:max(unique(gIX_iii))+1]);
        h.Values;
        hh(:,kk)=h.Values/length(find(index_all==ii));
        c{kk}=[num2str(ii) 'ex. fish'];kk=kk+1;
    end
    col=hsv(6);kk=1;
    figure,h=barh(hh,'stacked','FaceColor','flat');
    for ii=length(h)
        h(kk).CData = col(kk,:);kk=kk+1;
    end
    legend(c)
    ylim([0.5 max(unique(gIX_iii))+0.5]);ylabel('# Cluster');xlabel('Fraction');set(gca,'fontsize',15,'xcolor','w','ycolor','w');
  hh([4,7,10],[2 3 4 6])
    %%
   
l_fraction_in_region_in_clust=mean(fraction_in_region_in_clust(:,:,[1,2,3,4]),3,'omitnan');sum(l_fraction_in_region_in_clust,2,'omitnan');
% a=squeeze(sum(fraction_in_region_in_clust,1,'omitnan'));
% figure,h=bar([1:17],a(1:17,:),'stacked','xdata',[1:17]);set(gca,'fontsize',15,'xcolor','w','ycolor','w')

Lable={'l pallium rostral'; 'r pallium rostral';'l pallium middle';'r pallium middle';'l pallium lateral';'r pallium lateral';...
    'l habenula';'r habenula'; 'longitude';...
    'l tectum rostral';'r tectum rostral';'l tectum middle'; 'r tectum middle';'l tectum caudal';'r tectum caudal';...
    'l tegmentum';'r tegmentum';...
    'l cerebellum';'r cerebellum';...
    'l hindbrain';'r hindbrain';...
    }
c = categorical(Lable)';
cmap=clrmap_iii(num_clust,:);
y=l_fraction_in_region_in_clust;
figure,
%b=barh(c,y,'stacked','FaceColor','flat','YDataMode','manual','XDataMode','manual','xdata',[1:length( Lable)]);hold on;
b=barh([1:21],y,'stacked','FaceColor','flat');hold on;
kk=1;
for k =unique(gIX_iii)' %1:size(fraction_in_region_in_clust,2)
    b(k).CData = cmap(kk,:);
    kk=kk+1;
end
legend([b(unique(gIX_iii))],num2str(unique(gIX)),'location','bestoutside');
xlim([0 1]);set(gca,'xcolor','w','ycolor','w');ylim([0.5 21.5])
set(gca,'ytick',c,'XTickLabelRotation',45,'fontsize',12); axis ij;
%text([1:length(c)],num_in_region_in_clust{2}(id,1)+0.1*num_in_region_in_clust{2}(id,1),num2str(fraction_in_region_in_clust(id,unique(gIX)),'%.2f\n'));

a=fraction_in_region_in_clust(:,[3,14,4,7,10],1);%a(:,[11 12 13])=[];
a=squeeze(fraction_in_region_in_clust(:,10,:));%a(:,[11 12 13])=[];
figure,scatter(a)