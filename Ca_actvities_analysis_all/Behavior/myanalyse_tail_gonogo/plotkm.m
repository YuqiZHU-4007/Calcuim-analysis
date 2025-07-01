%Ncluster1=Ncluster+1
Ncluster1=4;
matrix=[norm.a1 norm.a3 norm.b2 norm.b3 norm.c3 norm.d3 norm.g3];
%matrix=[norm.a1-norm.a2 norm.b1(:,a:b) norm.b2 norm.b3 norm.c3 norm.d3 norm.d4 norm.g1(:,a:b) norm.g2(:,a:b) norm.g3 norm.g4];%norm.f2  norm.f1(:,1:norm.f4)  
matrixx=[norm.a1 norm.a3  norm.b2 norm.b3 norm.c1 norm.c3 norm.d1 norm.d3  norm.d4  norm.g3 norm.g4];%norm.f2 norm.g1(:,1:norm.f4)
[indk1,center]=kmeans(matrix,Ncluster1,'distance','sqEuclidean','start','cluster','display','final');
%indk=indsequence(indk1);%figure,histogram(indh1);%figure,histogram(indk1);
indk=indk1;
figure,
histogram(indk);title (['Cluster-kmeans:' segpara ',' seg]);

ind=indk;
if(Ncluster1>12)
    Color=linspecer(Ncluster1);
else
    Color=linspecer(Ncluster1,'qualitative');
end

%saveas(gcf,[savepath,'compare-hist','.png']);
% figure;title (['Cluster-Kmeans' segpara ',' seg]);
% for i=1:Ncluster
%     scatter(matrix1(ind==i,1),matrix1(ind==i,2),15,Color(ind(i),:),'filled');
%     hold on
% end
% % plot(c(:,1),c(:,2),'kx',...
% %      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Cluster 3');%,'Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8','Cluster 9'
% hold off
% saveas(gcf,[savepath,'c21','.png']);

figure,%帧数分布图
title(['result of kmeans-' segpara ',' seg]);
line([sti_start sti_start],[0 43],'Color','r');hold on
line([sti_start+duration sti_start+duration],[0 43],'Color','r');hold on
xlabel('frames');ylabel('trial number');
for i=1:size(para.marker,1)
    scatter(t_dur_cell{para.marker(i,1),1}(para.marker(i,2),1),para.marker(i,1),15,Color(ind(i),:),'filled');hold on;
    grid on;
end
hold off;
%interval=10;
saveas(gcf,[savepath,'c22','.png']);

% figure,%%%分类统计图
% C = categorical(ind,(1:Ncluster),{'1','2','3'});%,'4','5','6','7','8','9'
% histogram(C,'BarWidth',0.5,'Normalization','probability');title(['class statistics' segpara ',' seg]);
% saveas(gcf,[savepath,'c23','.png']);

figure,%%原数数据分类结果
for i=1:Ncluster1
    a=para.marker(find(ind==i),:);%(trial,movnumber)
    c=find(ind==i);
    c=c(1);
    for j=1:size(a,1)
        t=t_dur_cell{a(j,1),1};%t=t(find(a(:,1)==a(j,1)),:);
        plot(taildeg_everytrial(a(j,1),t(a(j,2),1):t(a(j,2),3))/15+10*(i-1),'color',Color(ind(c),:));hold on;%/10+interval*(i-1)
    end
end
title(['result of kmeans-' segpara ',' seg]);
xlabel('frames');ylabel('rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c24','.png']);

%%%均值对比
summ=zeros;
figure,
for i=1:Ncluster1
    a=para.marker(find(ind==i),:);%(trial,movnumber)
    c=find(ind==i);
    c=c(1);
    for j=1:size(a,1)
        t=t_dur_cell{a(j,1),1};sum1=zeros(1,max(t(:,2)));
        %t=t(find(a(:,1)==a(j,1)),:);
        sum1=sum1(1:t(a(j,2),2))+taildeg_everytrial(a(j,1),t(a(j,2),1):t(a(j,2),3));
    end
    summ(i,1:length(sum1))=sum1;
    plot(summ(i,:)/size(a,1),'color',Color(ind(c),:));hold on;%/10+interval*(i-1)
end
title(['Mean result of kmeans-' segpara]);
xlabel('frames');ylabel('mean rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c25','.png']);

%%%%组内组间显著性差异
figure,
a={};hin=zeros;hout=zeros;hout1=zeros;
for i=1:Ncluster1
     a{i}=find(ind==i);
     if length(a{i})==1
         hin(i)=0;
     else
     hin(i)=mean(pdist(matrix(a{i},:),'euclidean'));
     end
end
hin(:,find(hin==0))=mean(hin);
for i=1:Ncluster1
     aa=a{i};
     for k=1:Ncluster1
         bb=0;
         for j=1:size(aa,1)
             m=[matrix(aa(j),:);matrix(a{k},:)];
             b=pdist(m,'euclidean');b=b(1:length(a{k}));
             bb=bb+sum(b);
         end
         bb=bb/(size(aa,1)*length(a{k}));
         hout(i,k)=bb;
     end
     hout1(i)=mean(hout(i,:));
end
for i=1:Ncluster1
    scatter(1,hin(i),15,'r','filled');hold on;
    scatter(3,hout1(i),15,'b','filled');hold on;%Color(i,:)
    line([1 3],[hin(i) hout1(i)]);hold on;
end
legend('within','between');
[h,p,ci,stats]=ttest2(hin,hout1,0.05);
if h==1
    drawbrace([1 max([hin hout1])+0.3 ],[3 max([hin hout1])+0.3],10, 'Color', 'k');hold on;
    %line([1 3],[max([hin hout1])+0.3 max([hin hout1])+0.3],'Color','k','LineWidth',2);hold on;
    scatter(2,max([hin hout1])+0.5,15,'k','*');hold on;
end
hold off;
xlim([0 4]);axis equal;
title(['difference between and within groups(k-means)' segpara ]);
saveas(gcf,[savepath,'c27','.png']);

% %%%%%%%%%%%%%%%%对比
% % figure,subplot(3,1,1);scatter(matrix1(:,1),matrix1(:,2),15,'r','filled');
% % subplot(3,1,2);scatter(matrix1(:,1),matrix1(:,2),15,Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
% % subplot(3,1,3);scatter(matrix1(:,1),matrix1(:,2),15,Color(indk,:),'filled');title(['Cluster-Kmeans:' segpara]);
% % saveas(gcf,[savepath,'compare','.png']);