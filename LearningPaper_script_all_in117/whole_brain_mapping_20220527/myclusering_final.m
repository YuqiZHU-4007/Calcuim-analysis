function [clusterind,explained]=myclusering_final(A,M,cluster_type,clusternum,n_pc)
% A:n X t; n:neuron number; t:time series
% Allen_2018: Z-score,PCA(top 20 PCs), shared nearest neighbor graph-based clustering
% Vaziri_2020: Z-score, t-SNE(2D), Hierarchical clustering
% Orger_2014: Z-score,PCA(top 10 PCs),k-means (n=25)
%1、z-score
%2、降维pca/t-sne
%3、cluster

%% 3、autocluster/Hierarchical clustering/k++/shared nearest neighbor graph-based clustering
if cluster_type==1
    Y = pdist(M);Y=squareform(Y);
    Z = linkage(Y,'ward','euclidean' );figure,dendrogram(Z);c = cophenet(Z,Y);
else
    Y=[];Z=[];
end
switch cluster_type
    case 1 %Hierarchical clustering
        clusterind=cluster(Z,'maxclust',clusternum);
    case 2
        [clusterind,~,~ ]= kmeans(M,clusternum,'Distance','sqeuclidean','Replicates',1000,'MaxIter',1000);
end

explained=[];kk=1;ss=[];
if ~isempty(A)
    for ii=unique(clusterind)'
        m=A(:,find(clusterind==ii));
        [~,score,~,~,a,~] = pca(m);
        explained(1:min(n_pc,size(a,1)),kk)=a(1:min(n_pc,size(a,1)));
        ss(kk)=size(m,2);
        kk=kk+1;
    end
end


%% plot test
% %raw result
% load('H:\1.Test US\2.Tail free——Data from 117\20210805\fish1\para.mat');
% ref_win=ceil(0/fs.ca+1):frame.cs_start-1;%取时间段的baseline作为参考判断event
% time=([1:frame.per_cycle]-frame.cs_start)*fs.ca;
% area_win_hab=frame.cs_start:frame.cs_end-1;
% stimCS=zeros(1,frame.per_cycle*trial.total(1));
% for ii=1:size(stimCS,2)/frame.per_cycle
%     stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=0;%23
% end
% stimUS=ones(1,frame.per_cycle*trial.acq(1))*3;
% for ii=trial.acq(2):trial.acq(3)
%     stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=0;%16
% end
%
% figure('position',[680,558,946,420]),subplot(1,2,1),gscatter(M(:,1),M(:,2));xlabel('Tsne-1');ylabel('Tsne-1');set(gca,'fontsize',16);
% subplot(1,2,2), gscatter(M(:,1),M(:,2),clusterind);xlabel('Tsne-1');ylabel('Tsne-1');set(gca,'fontsize',16);
% %         [~,a]=sort(clusterind);
% %         figure,imagesc(A(:,a)',[0,8]);colormap('hot');colorbar;
% cIX=1:size(A,2);clrmap = GetColormap('hsv_new',max(clusterind));
% [h]=pushbutton_popupplot_Callback(A',cIX,clusterind,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
% %cut small cluster
% [cluster_ind_cut,cluster_ind_cut_remat,ind_cut,ind_save,value_cut,value_save]=cut_small_cluster(clusterind,10);
% disp(['saved neuron num. ',num2str(length(cluster_ind_cut)),'-saved cluster num. ',num2str(length(value_save))])
% figure, gscatter(M(ind_save,1),M(ind_save,2),cluster_ind_cut);hold on;scatter(M(ind_cut,1),M(ind_cut,2),12,[0.5 0.5 0.5],'filled');
% cIX=1:size(A(:,ind_save),2);
% [h]=pushbutton_popupplot_Callback(A(:,ind_save)',cIX,clusterind(ind_save),clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);

