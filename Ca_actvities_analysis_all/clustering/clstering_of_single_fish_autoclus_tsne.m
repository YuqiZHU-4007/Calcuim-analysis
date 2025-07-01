addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
addpath(genpath('F:\DUlab\Matlab\bhtsne'));
addpath(genpath('F:\DUlab\Matlab\tSNE_matlab'));

act_all=reshape(graph.node_mean_act_training',frame.per_cycle,trial.total,[]);
act_all_acq=act_all(:,trial.acq(2):trial.acq(3),:);
act_all_acq_block_mean=mean(reshape(act_all_acq,frame.per_cycle,trial.acq_block_trial,trial.acq_block_num,[]),2);
act_all_acq_block_mean=reshape(act_all_acq_block_mean,frame.per_cycle,trial.acq_block_num,[]);
figure,plot(graph.node_mean_act_training(5,trial.hab(3)*frame.per_cycle+1:trial.acq(3)*frame.per_cycle))
figure,plot(act_all_acq(:,:,5))
figure,plot(act_all_acq_block_mean(:,:,5))

ref_win=ceil(0/fs.ca+1):frame.cs_start-1;%取时间段的baseline作为参考判断event
time=([1:frame.per_cycle]-frame.cs_start)*fs.ca;
area_win_hab=frame.cs_start:frame.cs_end-1;
stimCS=zeros(1,frame.per_cycle*trial.total(1));
for ii=1:size(stimCS,2)/frame.per_cycle
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=0;%23
end
stimUS=ones(1,frame.per_cycle*trial.acq(1))*3;
for ii=trial.acq(2):trial.acq(3)
    stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=0;%16
end
%% Acq.
frame_ind=1:frame.us_start-2;%%%%%%%%%%%%%%%%%%%!!!!!!!
%act_all_acq_block_mean(:,:,1)=[];
M_0=reshape(act_all_acq_block_mean(frame_ind,:,:),length(frame_ind)*size(act_all_acq_block_mean,2),[])';
M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
base=M_0(:,frame.cs_start-15:frame.cs_start-1); m=mean(base(:));sd=std(base(:));
M=M_0;

% M_0=reshape(M_0',length(frame_ind),size(act_all_acq_block_mean,2),[]);
% frame_ind=frame.cs_start:frame.us_start-2;
% M_0= M_0(frame_ind,:,:);
% figure('name','M_0'),plot(M_0(:,:,5))

% M=[];label_roi=[];label_block=[];
% for ii=1:size(M_0,2)
%     M=cat(3,M,M_0(:,ii,:));
%     label_roi=cat(2,label_roi,[1:size(M_0,3)]);
%     label_block=cat(2,label_block,ii*ones(1,size(M_0,3)));
% end
% M=reshape(M,length(frame_ind),[])';
% figure('name','M'),plot(M')
% onset=findonset(M',m+1*sd);
[gIXx,~,numK,~]=find_best_k_in_range(M,1:20);
length(unique(gIXx))
cIXx=(1:size(M,1))';
cIX_reg = (1:size(M,1))';%ind;
if isempty(numK)
    numK=20;
end
%% clusering
cIX=cIXx;gIX=gIXx;
masterthres=0.8;
para=struct;
para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',10);
for ii=1:20
    [cIX,gIX] = AutoClustering(cIX,gIX,M,cIX_reg,1);%1,para,1,masterthres
end
length(unique(gIX))
result_cluster_raw.cIX=cIX;result_cluster_raw.gIX=gIX;%!!!!!!!!!!!!!!!!!!!!!!!!
gIX=[gIX;(max(gIX)+1)*ones(size(setdiff(cIX_reg,cIX)))];
cIX=[cIX;setdiff(cIX_reg,cIX)];
%% t-SNE
numDims = 2; pcaDims = 40; perplexity = 200;
mappedA = tsne(M, [], numDims,[],perplexity ); % no_dims = 2

eucD = pdist(mappedA,'euclidean');
clustTreeEuc = linkage(eucD,'average');
cophenet(clustTreeEuc,eucD)
figure('name','hireachy_tree_CS');[h,nodes] = dendrogram(clustTreeEuc,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
gIX = cluster(clustTreeEuc,'criterion','distance','cutoff',4);
%sort
mean_resp_clus=[];
for ii=1:max(gIX)
    mean_resp_clus(ii)=mean(mean(M(find(gIX==ii),:)));
end
[~,id]=sort(mean_resp_clus,'descend');
gIX_s= gIX;
for ii=1:max(gIX)
    gIX_s(find(gIX==ii))=id(ii);
end
%colorbar
n = round(max(gIX_s)*1.1);
cmap = hsv(max(1,n));
l=0;
h=figure('name','Colorbar');
for zz=1:n
    ll=length(find(gIX_s==zz));
    x = [l l+ll l+ll l];
    y = [0 0 1 1];
    l=l+ll;
    patch(x,y,cmap(zz,:));
end
xlim([1 l]); ylim([0 1]);ax=gca;set(ax,'visible','off');
figure('name','t-SNE');gscatter(mappedA(:,1),mappedA(:,2),gIX_s,cmap,'.',[],'on','t1','t2');axis equal;
%plot trace
% A=[];B=reshape(M_norm',frame.per_cycle,size(act_all_acq_block_mean,2),[]);
% for ii=1:size(B,2)
%     A=cat(3,A,B(:,ii,:));
% end
% A=reshape(A,frame.per_cycle,[])';
[~,gIx_s_ind]=sort(gIX_s);
frame_ind=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
A=reshape(act_all_acq_block_mean(frame_ind,:,:),length(frame_ind)*size(act_all_acq_block_mean,2),[])';
figure('name','imagesc');imagesc(A(gIx_s_ind,:),[min(A(:)) 0.05]);colorbar;hold on;%xlim([frame.cs_start,frame.cs_end]);
%scatter(onset(gIx_s_ind)+(frame.cs_start-1),[length(onset):-1:1],10,'r','filled');
A=act_all_acq_block_mean;%reshape(act_all_acq_block_mean,frame.per_cycle*size(act_all_acq_block_mean,2),[])';
[h,~]=plot_test_2(A,1:size(M,1),gIX_s,frame,cmap,2,[0 0.05],true,true);
figure;
scatter3(graph.nodes(:,1),graph.nodes(:,2),graph.nodes(:,3),12,cmap(gIX_s,:),'filled');axis equal;
%% acq profile of clusters
ind_clu=[];
for ii=1:size(act_all_acq_block_mean,3)
    for jj=1:size(act_all_acq_block_mean,2)
        ind_clu(jj,ii)=gIX_s(find(label_roi== ii & label_block== jj));
    end
end
figure('name','ind_clu'),imagesc(ind_clu');colormap(cmap);colorbar('Ticks',[1:max(gIX_s)]);
%hireachy
D = pdist(ind_clu','euclidean');
clustTree = linkage(D,'average');
cophenet(clustTree,D)
figure('name','hireachy_tree_acq_profile');[h,nodes] = dendrogram(clustTree,0);
h_gca = gca;
h_gca.TickDir = 'out';
h_gca.TickLength = [.002 0];
h_gca.XTickLabel = [];
Ix = cluster(clustTree,'criterion','distance','cutoff',2);
%kmeans
[Ix,~,~,~]=kmeans(ind_clu',4,'Replicates',10);
%sort
[Ix_s,Ix_s_ind]=sort(Ix);figure('name','ind_clu_s'),imagesc(ind_clu(:,Ix_s_ind)');colormap(cmap);colorbar('Ticks',[1:max(gIX_s)]);
cmap=hsv(max(Ix));
[h,~]=plot_test_2(act_all_acq_block_mean(:,:,Ix_s_ind),1:size(M_norm,1),Ix_s,frame,cmap,2,[-0 0.04],true,true);
A=reshape(act_all_acq_block_mean,frame.per_cycle*size(act_all_acq_block_mean,2),[])';
figure;imagesc(A(Ix_s_ind,:),[-0 0.1]);colorbar;
figure;
scatter3(graph.nodes(Ix_s_ind,1),graph.nodes(Ix_s_ind,2),graph.nodes(Ix_s_ind,3),12,cmap(Ix_s,:),'filled');axis equal;

%% GC
mask=getmask_imfreehand(env.vol(:,:,1:24),6);
figure,imshow(mask',[min(mask(:)) max(mask(:))]);
label=unique(mask);node_label_region=zeros(size(graph.nodes,1),1);
node=sub2ind(size(mask),ceil(graph.nodes(:,2)),ceil(graph.nodes(:,1)));
for ii=1:length(label)
    reg=find(mask==label(ii));
    Lia = ismember(node,reg);
    node_label_region(Lia == 1)=label(ii);
end
mean_act2=[];mean_act=[];
for ii=unique(node_label_region)'
    A=graph.node_mean_act2(find(node_label_region==ii),:);
    mean_act2(:,ii+1)=mean(A,1);
    A=graph.node_mean_act(find(node_label_region==ii),:);
    mean_act(:,ii+1)=mean(A,1);
end
[GC_nonzero_act1,GCstrength_act1]=Groupar_mls_concatenated(mean_act,3103,1);
[GC_nonzero_act2,GCstrength_act2]=Groupar_mls_concatenated(mean_act2,3095,1);
for ii=1:7
    for jj=1:7
        if ii==jj
        GCstrength_act1(ii,jj)=inf;
         GCstrength_act2(ii,jj)=inf;
        end
    end
end
figure('name',['GC_strength spon']);
set(gcf,'position',[50,200,1100,200]);
subplot(1,4,1),imagesc(GCstrength_act1);title('pre spon');
subplot(1,4,2),imagesc(GCstrength_act2);title('aft spon');
subplot(1,4,3),imagesc([(GCstrength_act2-GCstrength_act1) > 0]);title('increased');
subplot(1,4,4),imagesc([(GCstrength_act1-GCstrength_act2) > 0]);title('decreased');
figure('name',['GC_nonzero spon']);
set(gcf,'position',[50,200,1100,200]);
subplot(1,4,1),imagesc(GC_nonzero_act1);title('pre spon');
subplot(1,4,2),imagesc(GC_nonzero_act2);title('aft spon');
subplot(1,4,3),imagesc([(GC_nonzero_act2-GC_nonzero_act1) > 0]);title('increased');
subplot(1,4,4),imagesc([(GC_nonzero_act1-GC_nonzero_act2) > 0]);title('decreased');
