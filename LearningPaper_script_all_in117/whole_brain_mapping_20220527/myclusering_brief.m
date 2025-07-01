function [clusternum,clusterind,M,explained_all,n_pc,n_cluster]=myclusering_brief(A,dimention_reducntion_type,cluster_type)
%zyq,20220903
% A:n X t; n:neuron number; t:time series
% Allen_2018: Z-score,PCA(top 20 PCs), shared nearest neighbor graph-based clustering
% Vaziri_2020: Z-score, t-SNE(2D), Hierarchical clustering
% Orger_2014: Z-score,PCA(top 10 PCs),k-means (n=25)
%1¡¢z-score
%2¡¢½µÎ¬pca/t-sne
%3¡¢cluster
switch dimention_reducntion_type
    case 1 %pca
        PCnum=20;
        [coeff,score,latent,tsquared,explained,mu2] = pca(A','Rows','pairwise');
        explained_var = cumsum(latent)/sum(latent);a=find(explained_var>0.8, 1 );
        if explained_var(PCnum)>0.8 M=score(:,1:PCnum); else M=score(:,1:a); end
    case 2 %tsne
        if ceil(size(A,2)/35)<1000
            M = tsne(A','Algorithm','exact','Distance','euclidean');
            figure, gscatter(M(:,1),M(:,2))
        else
            for p=ceil(size(A,2)/35):ceil(size(A,2)/100);                
                M = tsne(A,'Algorithm','exact','Distance','euclidean','Perplexity',p);
            end
        end
    case 3 %no-dc
        M=A;
end
%3¡¢autocluster/Hierarchical clustering/k++/shared nearest neighbor graph-based clustering
switch cluster_type
    case 1 %Hierarchical clustering
        explained_all={};
        Y = pdist(M);Y=squareform(Y);
        Z = linkage(Y,'ward','euclidean' );figure,dendrogram(Z);c = cophenet(Z,Y);
        n_pc=3;n_cluster=3;
        for n=3:size(A,2)
            T=cluster(Z,'maxclust',n);explained=[];kk=1;ss=[];
            for ii=unique(T)'
                m=A(:,find(T==ii));
                [~,score,~,~,a,~] = pca(m);
                explained(1:min(n_pc,size(a,1)),kk)=a(1:min(n_pc,size(a,1)));
                ss(kk)=size(m,2);
                kk=kk+1;
            end
            a=explained(1,:)>50;b=sum(explained,1)>80; c=find(ss>n_cluster);
            explained_all{n}= explained;
            if sum(a(c) | b(c))==length(c) && ~isempty(c)
                break;
            else
                continue;
            end
        end
        clusternum=n;
        clusterind=cluster(Z,'maxclust',n);

    case 2 %kmeans
        n_pc=3;n_cluster=3;
        for n=3:size(A,2)
            [T,~,~ ]= kmeans(M,n,'Distance','sqeuclidean','Replicates',10,'MaxIter',1000);explained=[];kk=1;ss=[];
            for ii=unique(T)'
                m=A(:,find(T==ii));
                [~,score,~,~,a,~] = pca(m);
                explained(1:min(n_pc,size(a,1)),kk)=a(1:min(n_pc,size(a,1)));
                ss(kk)=size(m,2);
                kk=kk+1;
            end
            a=explained(1,:)>50;b=sum(explained,1)>70; c=find(ss>n_cluster);
            explained_all{n}= explained;
            if sum(a(c) | b(c))==length(c) && ~isempty(c)
                break;
            else
                continue;
            end
        end
        clusternum=n;
        [clusterind,~,~ ]= kmeans(M,clusternum,'Distance','sqeuclidean','Replicates',10,'MaxIter',1000);
end
end