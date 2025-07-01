%% Fig.2 whole brain mapping using aglin data
% 1、 neuron recruit：CS,US,MOTION change
% 2、similarity
% 3、functional type
% 4、population act
%% 1、 neuron recruit：CS,US,MOTION change
%CS activation pattern

%% colcalization
%% conditions & compute trial_avg
%state:early/middle/late learnining stage
%cue:CS/CS-US
%behavior:corr/incorr
S=8;%state {'baseline','middle trian','learned'}
D=2;%behavior {'corr','incorr'}
C=3;%cue {'CS','US','Motion'}


%% 2、clustering
% save all fish act
clc;clear all;close all
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
kk=1;activities_dfdf_align_all=[];fishid=[];fishname=[];label={'L','C','NL','FL'};
for batchi=1:4
    for fishi=1:length(Path{batchi})
        path=Path{batchi}{fishi};
        [nn,~]=fileparts(path);[nn,n1]=fileparts(nn);[~,n2]=fileparts(nn);nn=[n2,n1];
        load(fullfile(path, '/activities_dfdf_align.mat'));
        load(fullfile(path, '/env.mat'));
        A=cat(2,activities_dfdf_align,activities_dfdf_align_go,activities_dfdf_align_nogo);
        graph.mean_act_catall=[];
        for ii=1:size(graph.nodes,1)
            graph.mean_act_catall(:,:,ii)=mean( A(:,:,graph.node_label==ii,:),3);
        end
         figure,
         clr=GetColormap('hsv_new',max(graph.node_label));
         for ii=1:max(graph.node_label)
             a=find(graph.node_label==ii);
           scatter3(env.supervoxel(a,1)*0.66,env.supervoxel(a,2)*0.66,env.supervoxel(a,3)*10,5,clr(ii,:),'filled');hold on;
         end
        axis equal;
%         figure,imagesc(reshape(A(:,:,graph.node_label==ii,:),75*24,[])')
%         figure,imagesc(reshape(A(:,:,graph.idx_functional_adj==ii,:),75*24,[])')
        fishid=cat(1,fishid,[string(kk*ones(size(graph.nodes,1),1)),string([1:size(graph.nodes,1)]'),string(repmat(label{batchi},size(graph.nodes,1),1))]);
        fishname=cat(1,fishname,repmat(nn,size(graph.nodes,1),1));
        activities_dfdf_align_all=cat(3,activities_dfdf_align_all,graph.mean_act_catall);
        save(fullfile(path, '/activities_dfdf_align.mat'),'graph','-append','-v7.3');
        kk=kk+1;
    end
end
A_r=reshape(activities_dfdf_align_all,size(activities_dfdf_align_all,1)*size(activities_dfdf_align_all,2),[])';
a=A_r;a(find(isnan(A_r)))=0;
%A_r = zscore(A_r,0,'all','omitnan') ;
[coeff_pca,score_pca,latent,~,explained_pca,~] = pca(a);
explained_var = cumsum(latent)/sum(latent);nPC=min(find(explained_var>0.9));
activities_dfdf_align_all_PC=score_pca(:,1:nPC);
data=activities_dfdf_align_all_PC;
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;
activities_dfdf_align_all_PCA_dist_kD = pdist2(data,data,'euc','Smallest',nPC); 
[activities_dfdf_align_all_PCA_dist1]=pdist(data);
activities_dfdf_align_all_PCA_dist2=squareform(activities_dfdf_align_all_PCA_dist1);

%activities_dfdf_align_all_tsne1 = tsne(A_r');
numDims = 2; pcaDims = nPC; perplexity = 50; theta = .5; alg = 'svd';
activities_dfdf_align_all_tsne2 = fast_tsne(A_r, numDims, pcaDims, perplexity, theta, alg);
data=activities_dfdf_align_all_tsne2;
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;%data = zscore(data); % 使用zscore函数进行标准化
[activities_dfdf_align_all_tsne_dist1]=pdist(data);
activities_dfdf_align_all_tsne_dist2=squareform(activities_dfdf_align_all_tsne_dist1);
activities_dfdf_align_all_tsne_dist_kD = pdist2(data,data,'euc','Smallest',nPC); 

save([savepath '\activities_dfdf_align_all'],'activities_dfdf_align_all','fishid','fishname','A_r',...
    'activities_dfdf_align_all_PC','nPC','coeff_pca','score_pca','explained_pca',...
    'activities_dfdf_align_all_PCA_dist1','activities_dfdf_align_all_PCA_dist2','activities_dfdf_align_all_PCA_dist_kD',...
    'activities_dfdf_align_all_tsne2','activities_dfdf_align_all_tsne_dist1','activities_dfdf_align_all_tsne_dist2','activities_dfdf_align_all_tsne_dist_kD','-v7.3');

%% load and clusering 
clc;clear all;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\activities_dfdf_align_all'],'activities_dfdf_align_all_tsne_dist2');
load([savepath '\activities_dfdf_align_all'],'activities_dfdf_align_all_tsne2');
%A(:,find(isnan(A(1,:,1))),:)=[];
%A_r = zscore(A_r,0,'all','omitnan') ;
figure,plot(mean(A_r,1,'omitnan'))
% **********************dbscan************************
data=activities_dfdf_align_all_tsne2;
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;%data = zscore(data); % 使用zscore函数进行标准化
diary  dbscanlog.txt
diary on
best_score=0;eps=0;
[best_eps, best_min_samples, best_score,eps]=my_dbscan(data,activities_dfdf_align_all_tsne_dist2,nPC,activities_dfdf_align_all_tsne_dist_kD);
diary off
best_eps=0.010;best_min_samples=20;
[idx_dbscan, corepts_dbscan] = dbscan(activities_dfdf_align_all_tsne_dist2,best_eps, best_min_samples,'Distance','precomputed');
max(idx_dbscan)
figure,gscatter(activities_dfdf_align_all_tsne2(:,1),activities_dfdf_align_all_tsne2(:,2),idx_dbscan); legend off;
title([num2str(best_eps),'-',num2str(best_min_samples)]);
save([savepath '\cluster_align_all'],'idx_dbscan', 'corepts_dbscan','best_eps', 'best_min_samples', 'best_score','eps','-v7.3');
% *********************dpc ******************************
[idx_dpc_cl,idx_dpc_halo]=my_cluster_dp(activities_dfdf_align_all_dist2,activities_dfdf_align_all_dist1');
save([savepath '\cluster_align_all'],'mdist','idx_dpc_cl','idx_dpc_halo','-append');
% **********************SNN*****************************
K=[10:5:200];result_snn=[];
for p=1:length(K)
  %If you want to choose centers manually, set AutoPick to 0, otherwise, number of centers
  result_snn(p)=SnnDpc(activities_dfdf_align_all_PC,[1:size(activities_dfdf_align_all_PC,1)],K(p),'AutoPick',3,'Distance',activities_dfdf_align_all_dist1,'Ui',true);
%   fprintf('K=%i,AMI=%f,ARI=%f,FMI=%f\n',K(p),result(p).ami,result(p).ari,result(p).fmi);
end
dashboard_snn=table([result_snn(:).ami]',[result_snn(:).ari]',[result_snn(:).fmi]',K','VariableNames',{'AMI','ARI','FMI','K'});
resultBest_snn=dashboard_snn(dashboard_snn.ARI==max([result_snn(:).ari]),:);
save([savepath '\cluster_align_all'],'result_snn', 'dashboard_snn','resultBest_snn','-v7.3');

%% visualization
clc;clear all;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
load([savepath '\activities_dfdf_align_all'],'activities_dfdf_align_all_tsne_dist2','activities_dfdf_align_all_tsne2','activities_dfdf_align_all','fishid','fishname');

best_eps=0.02;best_min_samples=20;
[idx_dbscan_pca, corepts_dbscan] = dbscan(activities_dfdf_align_all_PCA_dist2,best_eps,best_min_samples,'Distance','precomputed');
max(idx_dbscan_pca)
length(find(idx_dbscan_pca==-1))
figure,gscatter(activities_dfdf_align_all_tsne2(:,1),activities_dfdf_align_all_tsne2(:,2),idx_dbscan_pca); legend off;
title([num2str(best_eps),'-',num2str(best_min_samples)]);
save([savepath '\cluster_align_all_pca'],'idx_dbscan_pca', 'corepts_dbscan','best_eps', 'best_min_samples', 'best_score','eps','-v7.3');


res=[0.66,0.66,10];
load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
A=reshape(activities_dfdf_align_all,size(activities_dfdf_align_all,1),size(activities_dfdf_align_all,2)/3,[],size(activities_dfdf_align_all,3));
activities_dfdf_align=squeeze(A(:,:,1,:));
activities_dfdf_align_go=squeeze(A(:,:,2,:));
activities_dfdf_align_nogo=squeeze(A(:,:,3,:));
ll={'All','Correct','in Correct'};    
sessionx ={'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'};% categorical
label_group={'L','C','NL','FL'};
Fraction_in_region_type={};Loc_in_region={};Loc_in_region_cell={};
Num_in_region={};Brain_region_id={};Label_region={};
idx=idx_dbscan;
node_ind_in_group=cell(max(idx),4);node_name_in_group=cell(max(idx),4);node_supervolxeli_in_group=cell(max(idx),4);
neuron_ind_in_group_cat=cell(max(idx),4);neuron_supervolxeli_in_group_cat=cell(max(idx),4);
neuron_ind_in_group={};neuron_supervolxeli_in_group={};
for idxi=1:max(idx)
    ind=find(idx==idxi);
    %脑区统计
    for batchi=1:4
        node_ind_in_group{idxi,batchi}=fishid(ind(strcmp(fishid(ind,3),label_group{batchi})),:);
        node_name_in_group{idxi,batchi}=string(fishname(ind(strcmp(fishid(ind,3),label_group{batchi})),:));
    end
    for batchi=1:size(node_ind_in_group,2)
        for fishi=1:length(Path{batchi});
            path=Path{batchi}{fishi};nn=[Path{batchi}{fishi}(end-14:end-7),Path{batchi}{fishi}(end-5:end-1)];
            node_ind_in_group_in_fish=node_ind_in_group{idxi,batchi}(strcmp(node_name_in_group{idxi,batchi},nn),:);
            if ~isempty(node_ind_in_group_in_fish)
                load(fullfile(path, '/activities_dfdf_align.mat'),'graph');
                supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                a=double([node_ind_in_group_in_fish(:,2)]);ind_in_group_in_fish_in_neuron=[];bb=string;kk=1;
                for ii=unique(a)'
                    b=find(graph.node_label==ii);
%                     figure,scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
%                     scatter3(supervolxeli(b,1),supervolxeli(b,2),supervolxeli(b,3),10,'r','filled');
                    node_supervolxeli_in_group{idxi,batchi}(kk,:)=mean(supervolxeli(b,1:3),1);
                    bb=cat(2,string(repmat(label_group{batchi},length(b),1)),string(fishi*ones(length(b),1)),b);
                    ind_in_group_in_fish_in_neuron=cat(1,ind_in_group_in_fish_in_neuron,bb);
                    kk=kk+1;
                end
                neuron_ind_in_group{idxi,batchi,fishi}=ind_in_group_in_fish_in_neuron;
                neuron_ind_in_group_cat{idxi,batchi}=cat(1, neuron_ind_in_group_cat{idxi,batchi},ind_in_group_in_fish_in_neuron);
                supervolxeli_in_group_in_fish_in_neuron=supervolxeli(str2double(ind_in_group_in_fish_in_neuron(:,3)),1:3);
%                 figure,scatter3( supervolxeli_in_group_in_fish_in_neuron(:,1), supervolxeli_in_group_in_fish_in_neuron(:,2), supervolxeli_in_group_in_fish_in_neuron(:,3),'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
                neuron_supervolxeli_in_group{idxi,batchi,fishi}=supervolxeli_in_group_in_fish_in_neuron;
                neuron_supervolxeli_in_group_cat{idxi,batchi}=cat(1,neuron_supervolxeli_in_group_cat{idxi,batchi},supervolxeli_in_group_in_fish_in_neuron);
                gIX=ones(1,size(supervolxeli_in_group_in_fish_in_neuron,1))';clrmap = GetColormap('hsv_new',max(gIX));
                for fraction_type=1:2
                    switch fraction_type
                        case 1
                            temp_vol=supervolxeli(:,1:3);
                        case 2
                            temp_vol=temp_supervoxel;
                    end
                    [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervolxeli_in_group_in_fish_in_neuron,reg_mask,reg_name,reg_loc,temp_vol,temp_env,clrmap,nn);
                    Fraction_in_region_type{idxi,batchi,fraction_type}(:,fishi)=fraction_in_region_in_clust(:,1);
                    Loc_in_region{idxi,batchi,fishi,fraction_type}=loc_in_region_in_clust;
                    Loc_in_region_cell{idxi,batchi,fishi,fraction_type}=loc_in_region_cell;
                    Num_in_region{idxi,batchi,fraction_type,3}(:,fishi)= num_in_region_in_clust{1}(1,:);
                    Num_in_region{idxi,batchi,fraction_type,3}(:,fishi)= num_in_region_in_clust{2}(1,:);
                    Num_in_region{idxi,batchi,fraction_type,3}(:,fishi)= num_in_region_in_clust{3}(1,:);
                    Brain_region_id{idxi,batchi,fraction_type,fishi}=brain_region_id(:,1);
                    Label_region{idxi,batchi,fishi,fraction_type}=Label;
                end
            end
        end
    end
end
save([savepath '\cluster_align_all'],'node_ind_in_group', 'node_name_in_group','node_supervolxeli_in_group',...
    'neuron_ind_in_group_cat','neuron_supervolxeli_in_group_cat','neuron_ind_in_group','neuron_supervolxeli_in_group',...
    'Fraction_in_region_type','Loc_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region',...
    '-append');


regionxx=categorical(Label_region{1,1,1,1}(1,:));regionxx = reordercats(regionxx,Label_region{1,1,1,1}(1,:));
mean1= @(x)(mean(x,2,'omitnan'));
error1= @(x)(std(x,[],2,'omitnan')./sqrt(length(x)));error2= @(x)(mean(x,2,'omitnan')-std(x,[],2,'omitnan'));
figure,
clrmap=GetColormap('hsv_new',length(unique(idx)));
gscatter(activities_dfdf_align_all_tsne2(:,1),activities_dfdf_align_all_tsne2(:,2),idx,clrmap);legend off;
figure,
for batchi=[1,4,2,3]
    subplot(2,2,batchi);
    for ii=1:size(node_supervolxeli_in_group,1)
        a=node_supervolxeli_in_group{ii,batchi};
        scatter3(a(:,1),a(:,2),a(:,3),10,clrmap(ii,:),'filled');hold on;
    end
    title(label_group{batchi});axis equal
end
P_all={};
for idxi=1:max(idx)
    ind=find(idx==idxi);
    % ****************plot*******************
    h=figure('position',[1,41,1920,962]);
    tiledlayout(12,6);
    t.TileSpacing = 'compact';t.Padding = 'compact';
    %tsne
    nexttile(1,[2,1]);
    scatter(activities_dfdf_align_all_tsne2(:,1),activities_dfdf_align_all_tsne2(:,2),10,[0.5 0.5 0.5],'filled');hold on;legend off;
    gscatter(activities_dfdf_align_all_tsne2(ind,1),activities_dfdf_align_all_tsne2(ind,2),idx(ind));legend off;
    xlabel('tsne-1');ylabel('tsne-2');
    title(['Cluster. ',num2str(idxi),'- Neuron Num. ',num2str(length(ind))]);
    %heatmap
    for ii=1
        nexttile(2,[2,1]);
        a=squeeze(A(:,:,ii,ind));a=reshape(a,size(a,1)*size(a,2),[]);
        imagesc(a',[0.1 1]);colormap('hot');title(ll{ii});
    end
    for ii=2:3
        nexttile(12+(ii-1),[2,1]);
        a=squeeze(A(:,:,ii,ind));a=reshape(a,size(a,1)*size(a,2),[]);
        imagesc(a',[-0.5 0.5]);colormap('hot');title(ll{ii});
    end
    %mapback : co-local  
    kk=1;
    for batchi=[1,4,2,3]
         nexttile(2+kk,[2,1]);kk=kk+1;
         scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),8,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
         scatter3(node_supervolxeli_in_group{idxi,batchi}(:,1),node_supervolxeli_in_group{idxi,batchi}(:,2),node_supervolxeli_in_group{idxi,batchi}(:,3),12,[1 0 0],'filled');
         view([90 -90]);title([label_group{batchi},'-',num2str(size(node_supervolxeli_in_group{idxi,batchi},1))]);%axis equal
    end
    kk=1;
    for batchi=[1,4,2,3]
        nexttile(14+kk,[2,1]);kk=kk+1;
        scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),8,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
        scatter3(node_supervolxeli_in_group{idxi,batchi}(:,1),node_supervolxeli_in_group{idxi,batchi}(:,2),node_supervolxeli_in_group{idxi,batchi}(:,3),12,[1 0 0],'filled');
        view([110 30]);
    end
    %plot trace 
    %不同session画
    b=0;
     for ii=1:size(activities_dfdf_align,2);
          a=mean(squeeze(A(10:60,ii,:,ind)),3);
          b=max(max(a(:)),b);
     end
     b=b*1.5;
    clrmap=[0,0,0;1,0,0;0,0,1];
    for ii=1:size(activities_dfdf_align,2)/2;
        nexttile(25+12*(ii-1),[2,1]);
        a=squeeze(A(:,ii,:,ind));
        plot_trace_202405(a,[0 b],clrmap,ll);xlim([10,60]);%set(gca,'xlabel','off');
        title(sessionx{ii});
    end
    for ii=size(activities_dfdf_align,2)/2+1:size(activities_dfdf_align,2);
        nexttile(26+12*(ii-size(activities_dfdf_align,2)/2-1),[2,1]);
        a=squeeze(A(:,ii,:,ind));
        plot_trace_202405(a,[0 b],clrmap,ll);xlim([10,60]);%set(gca,'xlabel','off');%
        title(sessionx{ii});
    end
    %all/go/no-go分开画
    clrmap=GetColormap('hsv_new',size(A,2));
    b=0.3;
    for ii=1:3
        nexttile(27+12*(ii-1),[2,1]);%kk=kk+2;
        a=squeeze(A(:,:,ii,ind));plot_trace_202405(a,[0 b],clrmap,sessionx);xlim([10,60]);title(ll{ii});
    end
    %脑区占比
    clrmap=GetColormap('hsv_new',30);
    kk=1;yli=0;
    for batchi=[1,4,2,3]
        a=mean1(Fraction_in_region_type{idxi,batchi,1});
        yli=max(max(a(:)),yli );
    end
    for batchi=[1,4,2,3]
        nexttile(28+12*(kk-1),[2,3]);
        a=Fraction_in_region_type{idxi,batchi,1};
        b=bar(regionxx,mean1(a));hold on;
        for k = 1:size(b,2)
            b(k).CData = clrmap(k,:);
            b(k).FaceColor = 'flat';
            b(k).FaceAlpha = 0.5;
        end
        xtips1 = b.XEndPoints; ytips1 = b.YEndPoints;
        scatter(reshape(repmat(regionxx',1,size(a,2)),1,[]),reshape(a,1,[]),10,'k','filled');hold on;
        hold on;
        er = errorbar(regionxx,mean1(a),error1(a));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        %set(gca,'fontsize',16);
        test_p=string;zz=1;
        for ii=setdiff(1:4,batchi)
            b=Fraction_in_region_type{idxi,ii,1};
            for regioni=1:size(b,1)
                try
                 p=ranksum(a(regioni,:),b(regioni,:));
                 if  p<=0.05;
                    % disp([labels_typei_other{type},'----SESSION ',num2str(jj),'L-C','p<0.05!!!!!!']);
                    a=strcat(label_group{batchi},'-',label_group{ii},':',string(regionxx(regioni)));
                    test_p(zz)=a;zz=zz+1;
                     %text(xtips1(regioni),ytips1(regioni),a,'HorizontalAlignment','center','VerticalAlignment','top','fontsize',14);hold on;
                 end
                catch
                    continue;
                end
            end
            P_all{idxi,batchi,ii}=test_p;
        end
        text(xtips1(1),max(ytips1),test_p,'HorizontalAlignment','center','VerticalAlignment','top','fontsize',8);hold on; 
        kk=kk+1;
        title(label_group{batchi});ylim([0,yli*1.2]);
    end
    a=checkpath(fullfile(savepath,'clusteringall-dbscan'));
    saveas(h,[a,'\',num2str(idxi)],'jpeg');
    savefig(h,[a,'\',num2str(idxi)]);
    close(h);
end
save([savepath '\cluster_align_all'],'P_all', '-append');

%charact of each cluster




% find(isnan(A_r));
% barcodes=strcat(string(1*ones(size(A,3),1)),brain_region_id(:,2));
% %barcodes=double(barcodes);
% genes=repmat([1:size(A_r,2)]',1,2);
% matrix=A_r';
% writematrix(barcodes,'F:\DUlab\FC_analyse\Ca_actvities_analyze\US mapping\R\test\barcodes.csv');
% writematrix(genes,'F:\DUlab\FC_analyse\Ca_actvities_analyze\US mapping\R\test\genes.csv');
% writematrix(matrix,'F:\DUlab\FC_analyse\Ca_actvities_analyze\US mapping\R\test\matrix.csv');
% save('F:\DUlab\FC_analyse\Ca_actvities_analyze\US mapping\R\test\matrix.mat','matrix');