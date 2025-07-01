%%%%%%%%画图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%savepath=[pathname 'figcluster'];
h=pdist(matrix,'euclidean');%%%计算各行向量的距离
hfreq=linkage(h,'average');
Ncluster=input('输入类个数\n');%%%
indh1=cluster(hfreq,'maxclust',Ncluster);%
indh=indsequence(indh1);
figure,
histogram(indh);title (['Cluster-Hierarchy:' segpara ',' seg]);

figure,dendrogram(hfreq);title (['Cluster-Hierarchy:' segpara ',' seg]);%%%树状图
ind=indh;
if(Ncluster>12)
    Color=linspecer(Ncluster);
else
    Color=linspecer(Ncluster,'qualitative');
end
%%%%Cluster-Hierarchy

%figure,scatter(matrix,15,Color(cfreq,:),'filled');title 'Cluster-Hierarchy';
%figure,subplot(2,1,1),scatter(matrix1(:,1),ones(1,size(matrix1,1)),15,'r','filled');title(['Cluster-Hierarchy:' segpara]);
% subplot(2,1,2),scatter(matrix1(:,1),ones(1,size(matrix1,1)),15,Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
% saveas(gcf,[savepath,'c13','.png']);

%figure,scatter(para.marker(:,2),para.marker(:,1),15,Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
figure,%帧数分布图
title(['result of hierarchy-' segpara ',' seg]);
line([sti_start sti_start],[0 43],'Color','r');hold on
line([sti_start+duration sti_start+duration],[0 43],'Color','r');hold on
xlabel('frames');ylabel('trial number');
for i=1:size(para.marker,1)
    scatter(t_dur_cell{para.marker(i,1),1}(para.marker(i,2),1),para.marker(i,1),15,Color(ind(i),:),'filled');hold on;
    grid on;
end
hold off;
saveas(gcf,[savepath,'c14','.png']);
%interval=10;

% figure,%%%分类统计图
% C = categorical(ind,(1:Ncluster),{'1','2','3'});%,'4','5','6','7','8','9'
% histogram(C,'BarWidth',0.5,'Normalization','probability');title(['class statistics' segpara ',' seg]);
% saveas(gcf,[savepath,'c15','.png']);

figure,%%原数数据分类结果
for i=1:Ncluster
    a=para.marker(find(ind==i),:);%(trial,movnumber)
    c=find(ind==i);
    c=c(1);
    for j=1:size(a,1)
        t=t_dur_cell{a(j,1),1};%t=t(find(a(:,1)==a(j,1)),:);
        plot(taildeg_everytrial(a(j,1),t(a(j,2),1):t(a(j,2),3))/15+15*(i-1),'color',Color(ind(c),:));hold on;%/10+interval*(i-1)
    end
end
title(['result of hierarchy-' segpara ',' seg]);
xlabel('frames');ylabel('rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c16','.png']);

%%%均值对比
summ=zeros;
figure,
for i=1:Ncluster
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
title(['Mean result of hierarchy-' segpara]);
xlabel('frames');ylabel('mean rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c17','.png']);

%%%%组内组间显著性差异
figure,
a={};hin=zeros;hout=zeros;hout1=zeros;
for i=1:Ncluster
     a{i}=find(ind==i);
     if length(a{i})==1
         hin(i)=0;
     else
     hin(i)=mean(pdist(matrix(a{i},:),'euclidean'));
     end
end
hin(:,find(hin==0))=mean(hin);
for i=1:Ncluster
     aa=a{i};
     for k=1:Ncluster
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
for i=1:Ncluster
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
title(['difference between and within groups(hierarchy)' segpara]);
saveas(gcf,[savepath,'c18','.png']);