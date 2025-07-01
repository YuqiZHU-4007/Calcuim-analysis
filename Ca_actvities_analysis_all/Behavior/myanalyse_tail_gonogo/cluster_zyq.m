%%%%%%%%%%%%计算距离后进行标准化！！！
%1、计算特征点
%2、标准化
%3、计算距离，聚类
clear all;
global t_dur_cell
global taildeg_everytrial
global fs;
global N;%%
global trailnum
global sti_start
global duration
global seg
global Ncluster
seg='9.29,fish1';
fs=240;
N=12000;
trailnum=42;
sti_start=7105;duration=240*1;
Ncluster=2;%input('输入类个数\n');
disp('Pick files include "alltailpos";"alltheta";"backim"');
[filename, pathname] = uigetfile( {'*.mat'}, 'Pick a file','MultiSelect','on');%'backim.mat';'alltailpos.mat';'alltheta.mat'
load([pathname filename]);
savepath = pathname;
% alltailpos1=zeros(11,2,length(alltheta));
% alltailpos1(:,:,movind) =alltailpos;
% if size(alltailpos1,3)<12000*trailnum
%     alltailpos1(:,:,length(alltailpos1)+1:12000*trailnum)=0;
% end
% for i=1:trailnum
%     alltailpos_everytrial{i,1}=alltailpos1(:,:,1+(i-1)*12000:12000+(i-1)*12000);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%计算paraments(+左右对称，尾动中部，尾部)
getpara;
k=1;
for i=1:trailnum
    a=t_dur_cell{i,1};
    test=zeros;
    for j=1:size(a,1)
        para.marker(k,:)=[i,j];
        para.rawdata(k,1:a(j,2))=taildeg_everytrial(i,a(j,1):a(j,3));
        para.dur(k,:)=a(j,2);
        para.freq(k,1:a(j,2))=abs(fft(para.rawdata(k,1:a(j,2))))*2./a(j,2);
        k=k+1;
        %plot((0:floor(a(j,2)/2)-1)/fs,fourier(1:floor(a(j,2)/2)));hold on;
        %scatter(locs(loc),para.frenquency(i,j));hold on;
    end
end
[para.freqmax,para.freqmaxloc,para.marker]=findfreqmax(taildeg_everytrial,t_dur_cell,trailnum);
%%%%聚类变量
matrix=[para.freqmax para.freqmaxloc];
[number,row]=size(matrix);
%%%my h cluster
h=pdist(matrix,'euclidean');%%%计算各行向量的距离
%h=zscore(h);%%
hfreq=linkage(h,'average');
cfreq=cluster(hfreq,'maxclust',Ncluster);%
if(Ncluster>12)
    Color=linspecer(Ncluster);
else
    Color=linspecer(Ncluster,'qualitative');
end
figure,dendrogram(hfreq);title 'Cluster-Hierarchy';
figure,dendrogram(hfreq,0);title 'Cluster-Hierarchy';
figure,scatter(para.marker(:,2),para.marker(:,1),15,Color(cfreq,:),'filled');title 'Cluster-Hierarchy';
% figure,
% for i=1:Ncluster
%     for j=1:number
%         if(cfreq(j)==i)
%             hold on
%             plot(para.marker(j,2),para.marker(j,1),'o','MarkerFaceColor',Color(i,:),'MarkerEdgeColor',Color(i,:))%dis(j,1),dis(j,2)
%         end
%     end
% end
%%%画分布图，多少被分为一类
%%%k-means
[idx,c]=kmeans(matrix,Ncluster,'distance','sqEuclidean','start','cluster');%,'display','final'
figure,
for i=1:Ncluster
scatter(para.marker(idx==i,2),para.marker(idx==i,1),15,Color(i,:),'filled');
hold on
end
% plot(c(:,1),c(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2',...
%     'Location','NW')
title 'Cluster-Kmeans'
hold off
