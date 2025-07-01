%% load movement parament (theta & duration)
clc;clear all;close all;
path_summary='H:\1.Test US\2.Tail free¡ª¡ªData from 117\summary.mat';
load(path_summary);
n=55:83;
%Theta=nan(3,8,length(n));Dur=nan(3,8,length(n));env_loc=nan(8,3,2,length(n));
Theta=[];Dur=[];env_loc=[];
P=[];kk=1;index_UR_fishnum=[];alltheta_all=[];
for ii=n
    a=load(fullfile(Path(ii,:), 'Results_of_alltheta.mat'));
    if ~isempty(cell2mat(Data.env_loc(:,ii)'))
        Theta(:,:,kk)=cell2mat(Data.Theta_US_n(:,ii)')';
        Dur(:,:,kk)=cell2mat(Data.env_dur_us(:,ii)')';
        env_loc(:,:,:,kk)=cell2mat(Data.env_loc(:,ii)');
        for jj=1:size(env_loc,2)
            [p,x]=calcula_bout_para(a.alltheta,squeeze(env_loc(:,jj,1,kk)),squeeze(env_loc(:,jj,2,kk)));
            P(:,jj,:,kk)=squeeze(cell2mat(squeeze(struct2cell(p))));
            alltheta_all(:,jj,:,kk)=x';
        end
        index_UR_fishnum(ii)=kk;kk=kk+1;
    end
end
Path(find(index_UR_fishnum~=0));
size(alltheta_all)
alltheta_all_rep=reshape(alltheta_all,size(alltheta_all,1),[],size(alltheta_all,4));alltheta_all_rep=reshape(alltheta_all_rep,size(alltheta_all,1),[],1)';
%% PCA
size(P)
label=zeros(size(Theta));label(:,4:8,index_UR_fishnum([55,60]))=1;label=reshape(label,[],1);
%reshaoe
A=reshape(Theta,[],1);B=reshape(Dur,[],1);
Para=reshape(P,size(P,1),[],size(P,4));Para=reshape(Para,size(P,1),[],1);
Para=Para';
%Para=[Para;A';B']';
Para=normalize(Para,'zscore');
[coeff,score,latent,~,explained] = pca(Para);
explained_var = cumsum(latent)/sum(latent)
figure('position',[100,100,1200,430]),subplot(1,2,1),plot(explained_var,'-k', 'Marker','o','MarkerSize',5,'linewidth',2);hold on;
bar(explained/100,'FaceColor',[0.5 0.5 0.5]);
plot(min(find(explained_var>0.9)),explained_var(min(find(explained_var>0.9))),'-r', 'Marker','o','MarkerSize',5,'linewidth',2);hold off;
xlim([0 25]);xlabel('# PC','fontsize',15);ylabel('variance explained (%)','fontsize',15);
subplot(1,2,2),imagesc(coeff);colorbar;xlabel('# PC','fontsize',15);ylabel('Parameters','fontsize',15);

figure,scatter3(score(:,1),score(:,2),score(:,3));hold on; 
scatter3(score(find(label==1),1),score(find(label==1),2),score(find(label==1),3),'r','filled');hold on; 
xlabel('PC1');ylabel('PC2');zlabel('PC3');
figure,
kk=1;
for ii=1:4
    for jj=1:4
        subplot(4,4,kk),
        scatter(score(:,ii),score(:,jj));hold on;
        scatter(score(find(label==1),ii),score(find(label==1),jj),'r','filled');hold on;
        xlabel(['PC' num2str(ii)]);ylabel(['PC' num2str(jj)]);
        kk=kk+1;
    end
end
%% clustering
index_UR_2=find(~isnan(score(:,1)));
PC_num=2:4;
X=score(index_UR_2,PC_num);
range=1:20;
for i = 1:length(range),
    disp(['k-means k=' num2str(range(i))]);
    [ix,~,sum ]= kmeans(X,range(i),'Distance','sqeuclidean','Replicates',10,'Start','plus');
    silh = silhouette(X,ix,'sqeuclidean');
    Silh_mean(i) = mean(silh);
end
figure;plot(range,Silh_mean,'o-','color',[0.5 0.5 0.5]);

[~,ix] = max(Silh_mean);n=3;
col=hsv(n);GetColormap('greedy_hsv',n);
[idx_kmeans,C,sumd,D] = kmeans(X,n,'Distance','sqeuclidean','Replicates',10,'Start','plus');
figure,
plot3(C(:,1),C(:,2),C(:,3),'kx','Markersize',15);hold on;
scatter3(X(:,1),X(:,2),X(:,3),12,col(idx_kmeans,:),'filled');
xlabel(['PC' num2str(PC_num(1))]);ylabel(['PC' num2str(PC_num(2))]);zlabel(['PC' num2str(PC_num(3))]);
figure,
kk=1;
for ii=1:size(X,2)
    for jj=1:size(X,2)
        subplot(size(X,2),size(X,2),kk),
        gscatter(X(:,ii),X(:,jj),idx_kmeans);hold on;
        scatter(X(find(label(index_UR_2)==1),ii),X(find(label(index_UR_2)==1),jj),10,'k');hold on;
        xlabel(['PC' num2str(PC_num(ii))]);ylabel(['PC' num2str(PC_num(jj))]);
        kk=kk+1;legend off;
    end
end
legend('Cluster 1','Cluster 2','Cluster 3','Labeled')
unique(idx_kmeans(find(label(index_UR_2)==1)))
figure,
for ii=1:n
    a=1:size(alltheta_all_rep,2);
    b= mean(abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:)),'omitnan');
    errBar=std(abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:)),'omitnan')/sqrt(length(find(idx_kmeans==ii))-1);
    shadedErrorBar(a,b,errBar,'lineprops',col(ii,:));hold on;
end
figure,
for ii=1
    b=abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:))';
    plot(abs(alltheta_all_rep(index_UR_2(find(idx_kmeans==ii)),:))')
end
%% svm