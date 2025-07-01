function [clusternum,clusterind,M,Y,explained_all,n_cluster]=myclusering(A,dimention_reducntion_type,cluster_type,n_pc)
% A:n X t; n:neuron number; t:time series
% Allen_2018: Z-score,PCA(top 20 PCs), shared nearest neighbor graph-based clustering
% Vaziri_2020: Z-score, t-SNE(2D), Hierarchical clustering
% Orger_2014: Z-score,PCA(top 10 PCs),k-means (n=25)
%1、z-score
%2、降维pca/t-sne
%3、cluster
switch dimention_reducntion_type
    case 1
        PCnum=20;
        [coeff,score,latent,tsquared,explained,mu2] = pca(A','Rows','pairwise');
        explained_var = cumsum(latent)/sum(latent);a=find(explained_var>0.8, 1 );
        if explained_var(PCnum)>0.8 M=score(:,1:PCnum); else M=score(:,1:a); end
    case 2
        if size(A,2)<6000
            M = tsne(A','Algorithm','exact','Distance','euclidean','Perplexity',50);
            figure, gscatter(M(:,1),M(:,2))
        else
            % for p=ceil(size(A,2)/35):ceil(size(A,2)/100);
            M = tsne(A,'Algorithm','exact','Distance','euclidean','Perplexity',ceil(size(A,2)/100));
        end
    case 3
        M=A;
end
%figure,imagesc(M,[-1,1])
%% 3、autocluster/Hierarchical clustering/k++/shared nearest neighbor graph-based clustering
if cluster_type==1
    Y = pdist(M);Y=squareform(Y);
    Z = linkage(Y,'ward','euclidean' );figure,dendrogram(Z);c = cophenet(Z,Y);
else
    Y=[];Z=[];
end
t=[1,1];cluster_ind_cut=[];n_cluster=3;
while ~(~isempty(cluster_ind_cut) && (sum(t(:,2)>n_cluster)>=3)) && n_cluster<=50
    explained_all={};
    n=3;a=zeros(n,1);b=zeros(n,1);c=1;
    while ~(sum(a(c) | b(c))==length(c) && ~isempty(c)) && n<=min(size(A,2),size(A,2))
        disp(['--------------\\\\running:n_cluster:',num2str(n_cluster),'--n:',num2str(n)]);
        switch cluster_type
            case 1 %Hierarchical clustering
                T=cluster(Z,'maxclust',n);
            case 2
                [T,~,~ ]= kmeans(M,n,'Distance','sqeuclidean','Replicates',100,'MaxIter',1000);
        end
        explained=[];kk=1;ss=[];
        for ii=unique(T)'
            m=A(:,find(T==ii));
            [~,score,~,~,a,~] = pca(m);
            explained(1:min(n_pc,size(a,1)),kk)=a(1:min(n_pc,size(a,1)));
            ss(kk)=size(m,2);
            kk=kk+1;
        end
        a=explained(1,:)>50;b=sum(explained,1)>80; c=find(ss>n_cluster);
        explained_all{n}= explained;n=n+1;
    end
    clusternum=n;n_cluster_finnal=n_cluster;
    clusterind=cluster(Z,'maxclust',n);
    t=tabulate(clusterind);
    [cluster_ind_cut,~,~,~]=cut_small_cluster(clusterind,20);
    disp(['!!!!!!n_cluster:',num2str(n_cluster),'--n:',num2str(clusternum),'-clusterindcut',num2str(length(cluster_ind_cut))]);
    n_cluster=n_cluster+1;
end
1
end

%% plot test
% %raw result
% load('H:\1.Test US\2.Tail free——Data from 117\20210805\fish1\para.mat')
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
% figure, gscatter(M(:,1),M(:,2),clusterind);
% %         [~,a]=sort(clusterind);
% %         figure,imagesc(A(:,a)',[0,8]);colormap('hot');colorbar;
% cIX=1:size(A,2);clrmap = GetColormap('hsv_new',max(clusterind));
% [h]=pushbutton_popupplot_Callback(A',cIX,clusterind,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
% %cut small cluster
% [cluster_ind_cut,ind_cut,ind_save,value_cut]=cut_small_cluster(clusterind,n_cluster);disp(length(cluster_ind_cut))
% figure, gscatter(M(ind_save,1),M(ind_save,2),cluster_ind_cut);hold on;scatter(M(ind_cut,1),M(ind_cut,2),12,[0.5 0.5 0.5],'filled');

