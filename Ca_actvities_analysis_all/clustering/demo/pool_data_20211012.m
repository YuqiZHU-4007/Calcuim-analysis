%% === variables=========================
% FuncHeteroMat = [ASW_fish,...%%==1)peak amplitude,n_rep x n_varicosity
%     ASSR_fish,...%%==2)synchronous activity strength
%     ASAR_fish,...%%==3)decay activity strength
%      freq_saturation_fish,...%%==4)n_rep x n_varicosity
%     freq_threshold_fish,...%%==5)n_rep x n_varicosity,每个varicosity钙反应起始频率；大于10%
%     freq_range_fish,...%%==6)n_rep x n_varicosity,每个varicosity钙反应tuning range
%     freq_half_fish,...%%==7)n_rep x n_varicosity, each varicosity's inflection point
%     decay_tau_fish,...%%==8)Ca_response_decay time constant
%     roi_size_fish,...%%%==10)ROI size==============
%     ];
% SpatioMat = [path length,branch order, brain region ];
%% === pool data=========================
clear; clc; close all;
[Pathname] = uipickfiles('out','ch'); %'num',3,选择同一条鱼所有‘fish_ramp_2_32Hz.mat’files, pick by order
file_num = size(Pathname);
fish_ID_sum = [];%%====
fish_ID_path_sum = [];%%====
imageID_sum = [];%%==钙成像图像的编号==
Dff0Mat_sum = [];
DFF0_FR_sum = [];
FuncHeteroMat_sum = [];
SpatioMat_sum = [];
roi_loc_fish_sum = [];
for n_fish = 1:file_num(1)
    filename = strtrim(Pathname(n_fish,:));
    load(filename);
    roi_num_fish = size(FuncHeteroMat,1);
    fish_ID_sum = [fish_ID_sum;repmat(n_fish,roi_num_fish,1)];
    imageID_sum = [imageID_sum;imageID_fish];%%==钙成像图像的编号==
    Dff0Mat_sum = [Dff0Mat_sum;Dff0Mat_fish'];
    DFF0_FR_sum = [DFF0_FR_sum;DFF0_FR_fish'];
    FuncHeteroMat_sum = [FuncHeteroMat_sum;FuncHeteroMat];
    SpatioMat_sum = [SpatioMat_sum;SpatioMat];
    roi_loc_fish_sum = [roi_loc_fish_sum;roi_loc_fish];
end
fish_ID_path_sum = Pathname;
clearvars -except fish_ID_sum fish_ID_path_sum imageID_sum = [] Dff0Mat_sum DFF0_FR_sum FuncHeteroMat_sum SpatioMat_sum ...
    roi_loc_fish_sum
save('FuncHetero_sum.mat');
%% === plot data=========================

[~,I] = sort(FuncHeteroMat_sum(:,1));
B =FuncHeteroMat_sum(I,:);

%%
%z-score
FuncHeteroMat_sum_n=normalize(FuncHeteroMat_sum,1);
%pca
[coeff,score,latent,~,explained,mu]=pca(FuncHeteroMat_sum_n);
explained_var=cumsum(latent)/sum(latent);
figure,bar(explained/100,'FaceColor',[0.5 0.5 0.5]);hold on;plot(explained_var,'-k', 'Marker','o','MarkerSize',5,'linewidth',2);hold on;
plot(min(find(explained_var>0.8)),explained_var(min(find(explained_var>0.8))),'-r', 'Marker','o','MarkerSize',5,'linewidth',2);hold off;
xlim([0 10]);xlabel('# PC','fontsize',15);ylabel('Variance explained (%)','fontsize',15);
figure,scatter3(score(:,1),score(:,2),score(:,3));
%%clustering
X=Dff0Mat_sum;%score(:,1:min(find(explained_var>0.8)));
[~,~,~,K]=kmeans_opt(X,20,0.8,3);
%hierarchy
Y = pdist(X,'euclidean');Z = linkage(Y,'average');figure,dendrogram(Z,0);clust_1 = cluster(Z,'maxclust',K);
%k-mean
clust_2 = kmeans(X,K,'distance','sqeuclidean','Replicates',15,'MaxIter',1000,'Display','final','EmptyAction','drop');
%DBSCAN
clust_3 = dbscan(X,0.5,10);n=length(unique(clust_3))
clust_3(find( clust_3<0))=n;
%plot
figure,kk=1;
for ii=1:3
    switch ii
        case 1
            clust=clust_1;str='hierarchy';
        case 2
            clust=clust_2;str='k-means';
        case 3
            clust=clust_3;str='DBSCAN';
    end
    n=length(unique(clust));
    clr=GetColormap('hsv_new',n);clr_idx=clr(clust,:);
    subplot(3,2,kk),
    scatter3(score(:,1),score(:,2),score(:,3),10,clr_idx,'filled');title([str ' Using Euclidean Distance Metric']);
    subplot(3,2,kk+1),silhouette(X,clust,'Euclidean');
    kk=kk+2;
end
for ii=2
    switch ii
        case 1
            clust=clust_1;str='hierarchy';
        case 2
            clust=clust_2;str='k-means';
        case 3
            clust=clust_3;str='DBSCAN';
    end
    figure,kk=1;
    for jj=unique(clust)'
        id=find(clust==jj);
        subplot(3,ceil(length(unique(clust))/3),kk),plot(Dff0Mat_sum(id,:)');kk=kk+1;
        ylim([0 10]);
    end
        figure,kk=1;
    for jj=unique(clust)'
        id=find(clust==jj);
        subplot(3,ceil(length(unique(clust))/3),kk),scatter3(SpatioMat_sum(id,1),SpatioMat_sum(id,2),SpatioMat_sum(id,3));kk=kk+1;
        ylim([0 10]);
    end
end









