figure,dendrogram(hfreq);title (['Cluster-Hierarchy:' segpara ',' seg]);%%%树状图
figure,dendrogram(hfreq,0);title (['Cluster-Hierarchy:' segpara ',' seg]);
%figure,scatter(matrix,15,Color(cfreq,:),'filled');title 'Cluster-Hierarchy';
%figure,plot(matrix(:,1),para.curv_max_max_loc,15,'color',Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
%figure,scatter(para.marker(:,2),para.marker(:,1),15,Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
figure,%帧数分布图
title(['result of hierarchy-' segpara ',' seg]);
line([sti_start sti_start],[0 43],'Color','r');hold on
line([sti_start+duration sti_start+duration],[0 43],'Color','r');hold on
xlabel('frames');ylabel('trial number');
for i=1:size(para.marker,1)
    scatter(t_dur_cell{para.marker(i,1),1}(para.marker(i,2),1),para.marker(i,1),15,Color(indh(i),:),'filled');hold on;
end
hold off;
%interval=10;

figure,%%%分类统计图
C = categorical(indh,(1:Ncluster),{'1','2','3'});%,'4','5','6','7','8','9'
histogram(C,'BarWidth',0.5,'Normalization','probability');title(['class statistics' segpara ',' seg]);

figure,%%原数数据分类结果
title(['result of hierarchy-' segpara ',' seg]);
xlabel('frames');ylabel('rawdata of ang_10');
for i=1:Ncluster
    a=para.marker(find(indh==i),:);%(trial,movnumber)
    for j=1:size(a,1)
        t=t_dur_cell{a(j,1),1};%t=t(find(a(:,1)==a(j,1)),:);
        plot(taildeg_everytrial(a(j,1),t(a(j,2),1):t(a(j,2),3)),'color',Color(i,:));hold on;%/10+interval*(i-1)
    end
end
hold off;