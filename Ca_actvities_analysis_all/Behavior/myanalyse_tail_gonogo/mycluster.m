%%%%%%%%%%%%计算距离后进行标准化！！！
%1、计算特征点
%2、标准化
%3、计算距离，聚类
clear all;
close all;
global t_dur_cell
global alltailpos_everytrial
global taildeg_everytrial
global para
global norm
global fs;
global N;%%
global trailnum
global sti_start
global duration
global seg
global Ncluster
seg=' fish1,9.29';
fs=240;
N=12000;
trailnum=42;
sti_start=7105;duration=240*1;

set(0, 'DefaultFigureColor', 'white');
%%%
disp('Pick files include "alltailpos";"alltheta";"t_dur";"backim"');
filename={};
[filename, pathname] = uigetfile( {'*.mat'}, 'Pick a file','MultiSelect','on');%'backim.mat';'alltailpos.mat';'alltheta.mat'
if size(filename,2)~=1
    for i = 1: length(filename)
        load(strcat(pathname,filename{i}));
    end
else
    error('you forget choosing multi files');
    %load([pathname filename]);
end
savepath = pathname;

% for i=1:trailnum
%     t_dur_cell{i,1}(find(t_dur_cell{i,1}),:)
% end
alltailpos1=zeros(11,2,length(alltheta));
alltailpos1(:,:,movind) =alltailpos;
if size(alltailpos1,3)<12000*trailnum
    alltailpos1(:,:,length(alltailpos1)+1:12000*trailnum)=0;
end
for i=1:trailnum
    alltailpos_everytrial{i,1}=alltailpos1(:,:,1+(i-1)*12000:12000+(i-1)*12000);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%计算paraments(+左右对称，尾动中部，尾部)
getpara;
%%%%聚类变量%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getnorm;
% figure,histogram(norm.a1-norm.a2,50);title('ang10r-ang10l');
% figure,histogram(norm.b2,50);title('ang10mean');
% figure,histogram(norm.b3,50);title('ang10max');
% figure,histogram(norm.c3,50);title('angpeaknum');%jia peakmax and loc 
% figure,histogram(norm.d3,50);title('freqmax');%%jia loc
% figure,histogram(norm.g3,50);title('curvmaxmax1');%%jia loc
% figure,scatter(norm.a1,norm.a2,15,'r','filled');xlabel('ang10r');ylabel('ang10l');
% figure,scatter(norm.b2,norm.b3,15,'r','filled');xlabel('ang10mean');ylabel('ang10max');
% figure,scatter(max(norm.c1,[],2),norm.c3,15,'r','filled');xlabel('angpeakmax');ylabel('angpeaknum');
% figure,scatter(norm.d4,norm.d3,15,'r','filled');xlabel('freqmaxloc');ylabel('freqmax');
% figure,scatter(norm.g4,norm.g3,15,'r','filled');xlabel('curvmaxmax1loc');ylabel('curvmaxmax1');

a=1;b=30;
matrix=[norm.a1 norm.a3 norm.b2 norm.b3 norm.c3 norm.d3 norm.g3];
%matrix=[norm.a1-norm.a2 norm.b1(:,a:b) norm.b2 norm.b3 norm.c3 norm.d3 norm.d4 norm.g1(:,a:b) norm.g2(:,a:b) norm.g3 norm.g4];%norm.f2  norm.f1(:,1:norm.f4)  
matrixx=[norm.a1 norm.a3 norm.b1(:,a:b) norm.b2 norm.b3 norm.c1 norm.c3 norm.d1 norm.d3  norm.d4 norm.g1(:,a:b) norm.g2(:,a:b) norm.g3 norm.g4];%norm.f2 norm.g1(:,1:norm.f4)
matrix1=[para.ang10r para.ang10l ];%画图用ones(length(para.ang10r-para.ang10l),1)
segpara='para-6';%%%%%%%%%%%%
[number,row]=size(matrix);
%%%my h cluster
h=pdist(matrix,'euclidean');%%%计算各行向量的距离
hfreq=linkage(h,'average');
figure,dendrogram(hfreq);title (['Cluster-Hierarchy:' segpara ',' seg]);%%%树状图
saveas(gcf,[savepath,'c11','.png']);
figure,dendrogram(hfreq,0);title (['Cluster-Hierarchy:' segpara ',' seg]);
saveas(gcf,[savepath,'c12','.png']);
Ncluster=input('输入类个数\n');%%%
indh1=cluster(hfreq,'maxclust',Ncluster);%
%figure,histogram(indh);
%%%k-means
[indk1,center]=kmeans(matrix,Ncluster,'distance','sqEuclidean','start','cluster','display','final');

indh=indsequence(indh1);indk=indsequence(indk1);%figure,histogram(indh1);%figure,histogram(indk1);
figure,
subplot(2,1,1),
histogram(indh);title (['Cluster-Hierarchy:' segpara ',' seg]);
subplot(2,1,2),histogram(indk);title (['Cluster-kmeans:' segpara ',' seg])
saveas(gcf,[savepath,'compare-hist','.png']);

%%%%%%%%画图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%savepath=[pathname 'figcluster'];
figure,histogram(h);title(['distance:' segpara]);
saveas(gcf,[savepath,'pdist','.png']);

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


%%%%%%%%%%%%%%%%%%%%%%k-means%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ind=indk;

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
for i=1:Ncluster
    a=para.marker(find(ind==i),:);%(trial,movnumber)
    c=find(ind==i);
    c=c(1);
    for j=1:size(a,1)
        t=t_dur_cell{a(j,1),1};%t=t(find(a(:,1)==a(j,1)),:);
        plot(taildeg_everytrial(a(j,1),t(a(j,2),1):t(a(j,2),3))/15+5*(i-1),'color',Color(ind(c),:));hold on;%/10+interval*(i-1)
    end
end
title(['result of kmeans-' segpara ',' seg]);
xlabel('frames');ylabel('rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c24','.png']);

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
title(['Mean result of kmeans-' segpara]);
xlabel('frames');ylabel('mean rawdata of ang_10');
hold off;
saveas(gcf,[savepath,'c25','.png']);

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
title(['difference between and within groups(k-means)' segpara ]);
saveas(gcf,[savepath,'c27','.png']);

% %%%%%%%%%%%%%%%%对比
% % figure,subplot(3,1,1);scatter(matrix1(:,1),matrix1(:,2),15,'r','filled');
% % subplot(3,1,2);scatter(matrix1(:,1),matrix1(:,2),15,Color(indh,:),'filled');title(['Cluster-Hierarchy:' segpara]);
% % subplot(3,1,3);scatter(matrix1(:,1),matrix1(:,2),15,Color(indk,:),'filled');title(['Cluster-Kmeans:' segpara]);
% % saveas(gcf,[savepath,'compare','.png']);

function ind=indsequence(ind)
global Ncluster;
for i=1:Ncluster
    l(i)=length((find(ind==i)));
    loc(i,1:l(i))=find(ind==i);
end
l2=sort(l);
for i=1:length(l)
    locc=loc(i,:);locc(:,find(locc==0))=[];
    if length(find(l2==l(i)))~=1
        a=find(l2==l(i));
        b=find(l==l2(a(1)));
        for j=1:length(a)            
            locc=loc(b(j),:);locc(:,find(locc==0))=[];
            ind(locc,:)=a(j);
        end
    else
    ind(locc,:)=find(l2==l(i));
    end
end
end

function a=mynorm(a)
a=(a-min(a))./(max(a)-min(a));
end