%node definition
%zyq-20190911
%refer to Richard,2018,BioRxiv
% (1) spatial contiguity, (2) functional homogeneity,(3) inter-subject consistency.
clc;clear all;
isspon2=true; istraining=false;ismask=false;
isnew_8X3_data=true;

input=struct;
[actname,actpath]=uigetfile('H:\1.Test US\2.Tail free¡ª¡ªData from 117\20210514\fish2\*.mat','act');
[envname,envpath]=uigetfile([actpath '\*.mat'],'env');
input.savepath=checkpath([actpath 'graph_results']);
if isspon2
    if isnew_8X3_data
       a=load([actpath,actname]);
       a1=squeeze(a.activities_spon(:,:,1)); activities=a1;
       a2=squeeze(a.activities_spon(:,:,end));activities2=a2;
    else
    [actname2,actpath2]=uigetfile([actpath '\*.mat'],'act2');
    a1=load([actpath,actname]); activities=a1.activities;
    a2=load([actpath2,actname2]);activities2=a2.activities;
    end
end
if istraining
    [actname3,actpath3]=uigetfile([actpath '\*.mat'],'act_training_aft_process');
    a1=load([actpath,actname]); activities=a1.activities;
    a3=load([actpath3,actname3]);activities3=[a3.activities_preCS_dfdf_aftcorrect]';
end
load([envpath,envname]);
if ismask
    %% getmask
    temp=max(env.vol(:,:,1:24),[],3);
    figure,imshow(temp,[min(temp(:)) max(temp(:))]);
    mask=getmask_imfreehand(env.vol(:,:,1:24),5);
    figure,imshow(mask',[min(mask(:)) max(mask(:))]);
    save([input.savepath '\mask.mat'],'mask');
elseif  isregistotemp
    a=load('H:\3.Juvenile reference brain\registration to templete\ÄÔÇø·Ö¸î\segmentation_file_0525_DSregion_mask.mat');
else                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    a2=load([input.savepath '\mask.mat']);
end
if isregist_to_temp
    loc_in_temp=readtable(fullfile(actpath,'vol_env_spatialloc_warped_SyN_add_brainregion.csv'));
    input.spacial_location=table2array(loc_in_temp(:,1:3));
else
    input.spacial_location=env.supervoxel(:,1:3);
    %% adjust location:Z
    X=input.spacial_location;
    for ii=unique(X(:,3))'
        id=find(X(:,3)==ii);
        X(id,3)=(ii-1)*8/0.66+1;
    end
    input.spacial_location=X;
    % spacial_location=input.spacial_location;
    % save('E:\A_Data_lightsheet\Data_huc\20190905\fish2\spon1\20190905fish2spon_act.mat','activities','activities_dfdf','spacial_location','reg_mask','reg_name','-v7.3');
    
end

winsize=121;
dfdf=compute_dfdf(activities,winsize);
dfdf_m=dfdf;%smoothdata(dfdf,2,'movmean',3);
figure;
n=randi(size(activities,1),[1,1]);
subplot(3,1,1),plot(activities(n,:)');xlim([1 size(activities,2)]);
subplot(3,1,2),plot(dfdf(n,:)');xlim([1 size(activities,2)]);
subplot(3,1,3),plot(dfdf_m(n,:)');xlim([1 size(activities,2)]);
input.act=dfdf(:,(winsize-1)/2+1:end-(winsize-1)/2);

if isspon2
    dfdf=compute_dfdf(activities2,winsize);
    dfdf_m=dfdf;%smoothdata(dfdf,2,'movmean',3);
    figure;
    n=randi(size(activities2,1),[1,1]);
    subplot(3,1,1),plot(activities2(n,:)');xlim([1 size(activities2,2)]);
    subplot(3,1,2),plot(dfdf(n,:)');xlim([1 size(activities2,2)]);
    subplot(3,1,3),plot(dfdf_m(n,:)');xlim([1 size(activities2,2)]);
    input.act2=dfdf(:,(winsize-1)/2+1:end-(winsize-1)/2);
%     activities_dfdf=input.act2(:,1:3000);activities=activities(:,1:3000);
end
if istraining
    %dfdf=compute_dfdf(activities,winsize);
    figure;
    n=randi(size(activities3,1),[1,1]);
    plot(activities3(n,:)');xlim([1 size(activities3,2)]);
    input.act_training=activities3;
end
clear dfdf;clear dfdf_m;
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');

%% k-means of spacial coordinates
n=5000;
clrmap = GetColormap(clrmap_name,n);
X=input.spacial_location;
%X=normalize(X,1,'zscore');
[idx_spacial,Center_spacial,sumd,D] = kmeans(X,n,'Distance','sqeuclidean','Replicates',10);
within_clus=sqrt(sumd);m(1,1)=mean(within_clus);m(1,2)=std(within_clus);m(1,3)=length(within_clus);
inter_clus=pdist(Center_spacial,'euclidean')';m(1,4)=mean(inter_clus);m(1,5)=std(inter_clus);m(1,6)=length(inter_clus);
clr=clrmap(idx_spacial,:);
temp = pdist(Center_spacial,'euclidean');Dist = squareform(temp);
figure;
%X=input.spacial_location;
scatter3(X(:,1),X(:,2),X(:,3),12,clr,'filled');hold on
plot3(Center_spacial(:,1),Center_spacial(:,2),Center_spacial(:,3),'kx','MarkerSize',15,'LineWidth',3);
xlabel('X');ylabel('Y');zlabel('Z');axis equal;
title 'Cluster Assignments and Centroids';hold off;
figure;hist(idx_spacial,1:n);xlim([1 n]);

%% k-means of activities
n=200;
X=input.act;act_center=[];
for ii=unique(idx_spacial)'
    id=find(idx_spacial==ii);
    act_center(ii,:)=mean(X(id,:),1);
end
% X=activities;act_raw_center=[];
% for ii=unique(idx_spacial)'
%     id=find(idx_spacial==ii);
%     act_raw_center(ii,:)=mean(X(id,:),1);
% end
if isspon2
    X=input.act2;act_center2=[];
    for ii=unique(idx_spacial)'
        id=find(idx_spacial==ii);
        act_center2(ii,:)=mean(X(id,:),1);
    end
else
    act_center2=[];
end
if istraining
    X=input.act_training;act_center_training=[];
    for ii=unique(idx_spacial)'
        id=find(idx_spacial==ii);
        act_center_training(ii,:)=mean(X(id,:),1);
    end
else
    act_center_training=[];
end
Cor=corr(act_center');
figure,imagesc(Cor,[0.5 max(Cor(:))]);colormap('hot');colorbar;
figure,imagesc(Cor,[0 max(Cor(:))]);colormap('hot');colorbar;
[idx_functional,~,~,~] = kmeans(Cor',n,'Distance','correlation','Replicates',10);
figure;hist(idx_functional,[1:n]);xlim([1 n]);
figure;silhouette(Cor',idx_functional,'correlation');
[idxAs,idxAa_ind]=sort(idx_functional);
figure;imagesc(Cor(idxAa_ind,idxAa_ind),[0.5 max(Cor(:))]);colormap('hot');colorbar;
figure;imagesc(Cor(idxAa_ind,idxAa_ind),[0 1]);colormap('hot');colorbar;
%figure;Z = linkage(Cor(idxAa_ind,idxAa_ind)');dendrogram(Z,0)
%% adjust clusters
idx_functional_adj=idx_functional;
thr_divide=100;thr_merge=70;divided=ones(length(unique(idx_functional_adj')),1);a=[];
for rep=1
    divided_idx=[];ii=0;
    while ii<=max(idx_functional_adj)-1;%unique(idxAadj')
        ii=ii+1;
        id=find(idx_functional_adj==ii);n=1;a(ii,1)=length(id);
        %figure;scatter3(C(id,1),C(id,2),C(id,3),12,'filled');hold on
        [~,Cn,sumd,~] = kmeans(Center_spacial(id,:),1,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
        %idxn=ones(size(id));
        if max(sumd)>=thr_divide
            [idxn,Cn,sumd,~] = kmeans(Center_spacial(id,:),2,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
            [counts,~]=hist(idxn,[1:max(idxn)]);
            d_inter_clusters=squareform(pdist(Cn,'euclidean'));
            [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
            adj_merge=find(~isequal(merge_idx,merge_idy));
            while  max(sumd)>=thr_divide && isempty(adj_merge) && min(counts)>=2; %mad>=thr1 && length(id)>=3 && squareform(pdist(C(id,:),'euclidean'))<=thr2
                n=n+1;
                if n>length(id)
                    warning([num2str(ii) 'not divided well!']);break;
                else
                    [idxn,Cn,~,~] = kmeans(Center_spacial(id,:),n,'Distance','sqeuclidean','Replicates',10);
                    d_inter_clusters=squareform(pdist(Cn,'euclidean'));
                    [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
                    if ~isempty(merge_idx)
                        for jj=1:length(merge_idx)
                            id1=find(idxn==merge_idx(jj));id2=find(idxn==merge_idy(jj));
                            if ~isequal(id1,id2)
                                idxn([id1,id2])=min(merge_idx(jj),merge_idy(jj));
                            else
                                break;
                            end
                        end
                    end
                    [Cn,sumd] = FindCentroid_Direct(idxn,Center_spacial(id,:));sumd=sqrt(sumd);
                    d_inter_clusters=squareform(pdist(Cn,'euclidean'));
                    [merge_idx,merge_idy]=find(d_inter_clusters<=thr_merge);
                    adj_merge=find(~isequal(merge_idx,merge_idy));
                    [counts,~]=hist(idxn,[1:max(idxn)]);
                end
            end %while
            divided_idx=id;%union(divided_idx,id);
            adj=setdiff(1:length(idx_functional),union(id,divided_idx));
            adj2=find(idx_functional_adj(adj)> ii);
            %min(idxAadj(adj(adj2)))
            idx_functional_adj(adj(adj2))=idx_functional_adj(adj(adj2))+max(idxn)-1;
            %min(idxAadj(adj(adj2)))
            idx_functional_adj(id)=idxn+idx_functional_adj(id)-1;
            divided(ii,:)=n;
            if length(unique(idx_functional_adj))~=length([1:max(idx_functional_adj)]')
                disp true1
                disp (num2str(ii))
                break;
            end
        else
            divided_idx=id;%union(divided_idx,id);
            divided(ii,:)=n;
            %continue;
            if length(unique(idx_functional_adj))~=length([1:max(idx_functional_adj)]')
                disp false1
                disp (num2str(ii));
                break;
            end
        end
        %         figure;
        %         subplot(2,1,1);hist(idxA,1:max(idxA));xlim([1 max(idxA)])
        %         subplot(2,1,2);hist(idxAadj,1:max(idxAadj));xlim([1 max(idxAadj)])
        %         d=squareform(pdist(C(id,:),'euclidean'));
        %         [mad,~]=max(d(:));[maix,maiy]=find(d==mad);mai=unique([maix,maiy]);
        %         [miix,miiy]=find(d<=thr2);mii=unique([miix,miiy])
        %         figure;
        %         subplot(2,1,1);hist(idxA,1:max(idxA));
        %         subplot(2,1,2);hist(idxAadj,1:max(idxAadj));
    end
end
figure;
subplot(2,1,1);hist(idx_functional,1:max(idx_functional));xlim([1 max(idx_functional)]); max(idx_functional)
subplot(2,1,2);hist(idx_functional_adj,1:max(idx_functional_adj));xlim([1 max(idx_functional_adj)]);max(idx_functional_adj)
[Center_final,~] = FindCentroid_Direct(idx_functional_adj,Center_spacial);
figure;
clrmap = GetColormap(clrmap_name,length(unique(idx_functional_adj)));clr=clrmap(idx_functional_adj,:);
X=Center_spacial;
%scatter3(X(:,1),X(:,2),X(:,3),12,clr,'filled');hold on
plot3(Center_final(:,1),Center_final(:,2),Center_final(:,3),'k*','MarkerSize',5,'LineWidth',1.5);
xlabel('X');ylabel('Y');zlabel('Z');axis equal;
title 'Cluster Assignments and Centroids';hold off;
set(gca,'visible','off');
h=figure;idx=[];n=50:80;
for ii=n
    idx=[idx; find(idx_functional_adj==ii)];
end
clrmap = GetColormap(clrmap_name,length(n));clr=clrmap(idx_functional_adj(idx)-min(n)+1,:);
X=Center_spacial(idx,:);
[Cn,~] = FindCentroid_Direct(idx_functional_adj(idx),Center_spacial(idx,:));
scatter3(X(:,1),X(:,2),X(:,3),12,clr,'filled');hold on
plot3(Cn(:,1),Cn(:,2),Cn(:,3),'kx','MarkerSize',15,'LineWidth',3);
%scatter3(Cn(:,1),Cn(:,2),Cn(:,3),100,clrmap,'x','LineWidth',2);hold on
xlabel('X');ylabel('Y');zlabel('Z');%axis equal;
title 'Cluster Assignments and Centroids';hold off;
saveas(h,[checkpath(input.savepath) 'node_definition_zoom.fig']);
%% save
graph.nodes=Center_final;%[graph.nodes_counts,~]= hist(idx_functional_adj,1:max(idx_functional_adj));
graph.node_mean_act=[];
for ii=1:size(graph.nodes,1)
    graph.node_mean_act(ii,:)=mean(act_center(idx_functional_adj==ii,:),1);
end
if isspon2
    graph.node_mean_act2=[];
    for ii=1:size(graph.nodes,1)
        graph.node_mean_act2(ii,:)=mean(act_center2(idx_functional_adj==ii,:),1);
    end
else
    graph.node_mean_act2=[];
end
if istraining
    graph.node_mean_act_training=[];
    for ii=1:size(graph.nodes,1)
        graph.node_mean_act_training(ii,:)=mean(act_center_training(idx_functional_adj==ii,:),1);
    end
else
    graph.node_mean_act_training=[];
end
% graph.node_label=zeros(length(idx_spacial),1);
% for ii=unique(idx_functional_adj)'
%     ind=find(idx_functional_adj==ii);
%     for jj=unique(ind)'
%         graph.node_label(find(idx_spacial==jj))=ii;
%     end
% end
graph.node_label=zeros(length(idx_spacial),1);
for ii=1:length(idx_functional_adj)
  graph.node_label(find(idx_spacial==ii))=idx_functional_adj(ii);
end
% graph.node_mean_act_raw=[];
% for ii=1:size(graph.nodes,1)
%     graph.node_mean_act_raw(ii,:)=mean(act_raw_center(idxAadj==ii,:),1);
% end
graph.idx_spacial=idx_spacial;
graph.Center_spacial=Center_spacial;
graph.idx_functional=idx_functional;
graph.idx_functional_adj=idx_functional_adj;
%% label node
label=unique(mask);graph.node_label_region=zeros(size(graph.nodes,1),1);
node=sub2ind(size(mask),ceil(graph.nodes(:,2))+1,ceil(graph.nodes(:,1))+1);
for ii=1:length(label)
    reg=find(mask==label(ii));
    Lia = ismember(node,reg);
    graph.node_label_region(Lia == 1)=label(ii);
end
save([input.savepath '\graph.mat'],'graph');

% idxAadj=idxA;
% thr1=200;thr2=60;thres_merge=60;divided=ones(length(unique(idxAadj')),1);
% for rep=1
%     divided_idx=[];
%     for ii=unique(idxA)'
%         id=find(idxA==ii);n=1;
%         %figure;scatter3(C(id,1),C(id,2),C(id,3),12,'filled');hold on
%         [~,Cn,sumd,~] = kmeans(C(id,:),1,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
%         idxn=ones(size(id));
%         if max(sumd)>=thr1
%             [~,Cn,sumd,~] = kmeans(C(id,:),2,'Distance','sqeuclidean','Replicates',10);sumd=sqrt(sumd);
%             d_inter_clusters=squareform(pdist(Cn,'euclidean'));
%             [merge_idx,merge_idy]=find(d_inter_clusters<=thr2);
%             adj_merge=find(~isequal(merge_idx,merge_idy));
%             while  max(sumd)>=thr1 && isempty(adj_merge); %mad>=thr1 && length(id)>=3 && squareform(pdist(C(id,:),'euclidean'))<=thr2
%                 n=n+1;
%                 if n>length(id)
%                     warning([num2str(ii) 'not divided well!']);break;
%                 else
%                     [idxn,Cn,~,~] = kmeans(C(id,:),n,'Distance','sqeuclidean','Replicates',10);
%                     d_inter_clusters=squareform(pdist(Cn,'euclidean'));
%                     [merge_idx,merge_idy]=find(d_inter_clusters<=thr2);
%                     if ~isempty(merge_idx)
%                         for jj=1:length(merge_idx)
%                             id1=find(idxn==merge_idx(jj));id2=find(idxn==merge_idy(jj));
%                             if ~isequal(id1,id2)
%                                 idxn([id1,id2])=min(merge_idx(jj),merge_idy(jj));
%                             else
%                                 break;
%                             end
%                         end
%                     end
%                     [Cn,sumd,~] = FindCentroid_Direct(idxn,C(id,:));sumd=sqrt(sumd);
%                     d_inter_clusters=squareform(pdist(Cn,'euclidean'));
%                     [merge_idx,merge_idy]=find(d_inter_clusters<=thr2);
%                     adj_merge=find(~isequal(merge_idx,merge_idy));
%                 end
%             end %while
%             idxAadj(id)=idxn+max(idxAadj)-1;
%             divided(ii,:)=n;
%         else
%            continue;
%         end
%         %         d=squareform(pdist(C(id,:),'euclidean'));
%         %         [mad,~]=max(d(:));[maix,maiy]=find(d==mad);mai=unique([maix,maiy]);
%         %         [miix,miiy]=find(d<=thr2);mii=unique([miix,miiy])
% %         figure;
% %         subplot(2,1,1);hist(idxA,1:max(idxA));
% %         subplot(2,1,2);hist(idxAadj,1:max(idxAadj));
%     end
% end
% figure;
% subplot(2,1,1);hist(idxA,1:max(idxA));
% subplot(2,1,2);hist(idxAadj,1:max(idxAadj));
% figure;
% clrmap = GetColormap(clrmap_name,length(unique(idxAadj)));clr=clrmap(idxAadj,:);
% X=C;
% scatter3(X(:,1),X(:,2),X(:,3),12,clr,'filled');hold on
% xlabel('X');ylabel('Y');zlabel('Z');%axis equal;