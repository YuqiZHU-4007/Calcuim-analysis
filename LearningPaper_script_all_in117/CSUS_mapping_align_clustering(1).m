%% Fig.2 whole brain mapping using aglin data
% 1¡¢ neuron recruit£ºCS,US,MOTION change
% 2¡¢similarity
% 3¡¢functional type
% 4¡¢population act
%% 1¡¢ neuron recruit£ºCS,US,MOTION change
%% colcalization
%% conditions & compute trial_avg
%state:early/middle/late learnining stage
%cue:CS/CS-US
%behavior:corr/incorr
S=8;%state {'baseline','middle trian','learned'}
D=2;%behavior {'corr','incorr'}
C=3;%cue {'CS','US','Motion'}


%% 2¡¢clustering
clc;clear all;close all
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
kk=1;activities_dfdf_align_all=[];fishid=[];fishname=[];
for batchi=1:4
    for fishi=1:length(Path{batchi})
        path=Path{batchi}{fishi};
        [nn,~]=fileparts(path);[nn,n1]=fileparts(nn);[~,n2]=fileparts(nn);nn=[n2,n1];
        load(fullfile(path, '/activities_dfdf_align.mat'));
        A=cat(2,activities_dfdf_align,activities_dfdf_align_go,activities_dfdf_align_nogo);
        graph.mean_act_catall=[];
        for ii=1:size(graph.nodes,1)
            graph.mean_act_catall(:,:,ii)=mean( A(:,:,graph.idx_functional_adj==ii,:),3);
        end
        fishid=cat(1,fishid,kk*ones(size(graph.nodes,1),1));
        fishname=cat(1,fishname,repmat(nn,size(graph.nodes,1),1));
        activities_dfdf_align_all=cat(3,activities_dfdf_align_all,graph.mean_act_catall);
        save(fullfile(path, '/activities_dfdf_align.mat'),'graph','-append','-v7.3');
        kk=kk+1;
    end
end
A_r=reshape(activities_dfdf_align_all,size(activities_dfdf_align_all,1)*size(activities_dfdf_align_all,2),[])';
%A_r = zscore(A_r,0,'all','omitnan') ;
[coeff_pca,score_pca,latent,~,explained_pca,~] = pca(A_r);
explained_var = cumsum(latent)/sum(latent);nPC=min(find(explained_var>0.9));
activities_dfdf_align_all_PC=score_pca(:,1:nPC);
data=activities_dfdf_align_all_PC;
data=(data-min(data))./(max(data)-min(data));
data(isnan(data))=0;
activities_dfdf_align_all_dist1=squareform(pdist(data));

%activities_dfdf_align_all_tsne1 = tsne(A_r');
numDims = 2; pcaDims = nPC; perplexity = 50; theta = .5; alg = 'svd';
activities_dfdf_align_all_tsne2 = fast_tsne(A_r, numDims, pcaDims, perplexity, theta, alg);
save([savepath '\activities_dfdf_align_all'],'activities_dfdf_align_all','fishid','fishname','A_r','activities_dfdf_align_all_PC','nPC','coeff_pca','score_pca','explained_pca',...
    'activities_dfdf_align_all_dist1','activities_dfdf_align_all_tsne2','-v7.3');

%% 
clc;clear all;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\activities_dfdf_align_all']);
%A(:,find(isnan(A(1,:,1))),:)=[];
%A_r = zscore(A_r,0,'all','omitnan') ;
figure,plot(mean(A_r,1,'omitnan'))
data=activities_dfdf_align_all_PC;
%idx = dbscan(data,1,5);
K=[10:5:200];
for p=1:length(K)
  %If you want to choose centers manually, set AutoPick to 0, otherwise, number of centers
  result(p)=SnnDpc(data,[1:size(data,1)],K(p),'AutoPick',3,'Distance',activities_dfdf_align_all_dist1,'Ui',true);
%   fprintf('K=%i,AMI=%f,ARI=%f,FMI=%f\n',K(p),result(p).ami,result(p).ari,result(p).fmi);
end
dashboard=table([result(:).ami]',[result(:).ari]',[result(:).fmi]',K','VariableNames',{'AMI','ARI','FMI','K'});
resultBest=dashboard(dashboard.ARI==max([result(:).ari]),:);

%% visualization
%T-SNE
figure,gscatter(activities_dfdf_align_all_tsne2(:,1),activities_dfdf_align_all_tsne2(:,2),idx)
%plot trace
A=reshape(activities_dfdf_align_all,size(activities_dfdf_align_all,1),size(activities_dfdf_align_all,2)/3,[],size(activities_dfdf_align_all,3));
activities_dfdf_align=squeeze(A(:,:,1,:));
activities_dfdf_align_go=squeeze(A(:,:,2,:));
activities_dfdf_align_nogo=squeeze(A(:,:,3,:));
ll={'all','Correct','in Correct'};
figure,
for idxi=1:max(idx)
    ind=find(idx==idxi);
    figure,
    tiledlayout(9,4);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    clrmap=GetColormap('hsv_new',size(data,2));
    sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
    for ii=1:3
        a=squeeze(A(:,:,ii,ind));
        nexttile([3,1]),plot_trace_202405(a,[0 1],clrmap,sessionx);
        title(ll{ii})
    end
    clrmap=[0,0,0;1,0,0,0,0,1];
    for ii=1:size(activities_dfdf_align,2);
        a=cat(activities_dfdf_align(:,ii,ind),activities_dfdf_align_go(:,ii,ind),activities_dfdf_align_nogo(:,ii,ind))
        nexttile([1,1]),plot_trace_202405(a,[0 1]);
    end
end
%mapback
%co-local
%ÄÔÇøÍ³¼Æ
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