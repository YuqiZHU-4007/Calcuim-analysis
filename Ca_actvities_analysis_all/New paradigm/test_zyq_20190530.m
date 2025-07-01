
hdf5_dir='F:\DUlab\FC analyse\FishExplorer\test data\subject_6\TimeSeries.h5';
CellResp = h5read(hdf5_dir,'/CellResp');
CellRespZ = h5read(hdf5_dir,'/CellRespZ');
CellRespAvr = h5read(hdf5_dir,'/CellRespAvr');
CellRespAvrZ = h5read(hdf5_dir,'/CellRespAvrZ');
absIX = h5read(hdf5_dir,'/absIX');


%% setup dir
file_masterdir = 'F:\DUlab\FC analyse\';
if ~exist(file_masterdir,'dir')
    file_masterdir = 'C:\Users\xiuye\Dropbox';
end
code_dir = fullfile(file_masterdir,'FishExplorer');

% SET CUSTOM PATH example:
% code_dir = 'C:\Users\xiuye\Dropbox\FishExplorer-master';

addpath(genpath(code_dir));

cd(code_dir)


for kk=[8]
    clust = find(gIX2==kk);
    activities_preCS_dfdf_aftcorrect(:,cIX2(clust))=[];
    env.supervoxel(cIX2(clust),:)=[];
    area.CS_hab_tst(:,cIX2(clust))=[];
    area.CS_acq_block(:,cIX2(clust))=[];
    rawacti.acq_mean_blocks(cIX2(clust),:)=[];
    rawacti.hab(:,cIX2(clust))=[];
    rawacti.tst(:,cIX2(clust))=[];
    rawacti.acq_1st_block(cIX2(clust),:)=[];
    inter=intersect(ind_increased.pagetst_acq_nonadj,cIX2(clust));
    for ii=1:length(inter)
        iinddd=find(ind_increased.pagetst_acq_nonadj==inter(ii));
        if ~isempty(iinddd)
            ind_increased.pagetst_acq_nonadj(iinddd)=[];
        end
    end
end

cut_move=true;
ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.us_start-1)),1);
ind_cut_trial=unique(ind_cut_trial);
aa=[];
for ii=setdiff(trial.acq(2):trial.acq(3),intersect(trial.acq(2):trial.acq(3),ind_cut_trial))
a=activities_preCS_dfdf_aftcorrect((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.us_start-1,:)';
aa=[aa a];
end
stimCS=zeros(1,size(activities_preCS_dfdf_aftcorrect,1));
for ii=1:trial.total
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=23;
end
stimUS=ones(1,size(activities_preCS_dfdf_aftcorrect,1))*3;
for ii=trial.acq(2):trial.acq(3)
    stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=16;
end
find(stimUS==16);
M_0=aa;
%M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
cIX_reg = (1:size(M_0,1))'; 
numK=[];
[gIX,C,numK]=find_best_k_in_range(M_0,4);
cIX=1:length(gIX);
clrmap=pushbutton_popupplot_Callback(activities_preCS_dfdf_aftcorrect',cIX,gIX,env,fs.ca,stimCS,stimUS,0,6,0,0,0,0);
DrawTiledPics_zyq_20190530(cIX,gIX,[1:size(M_0,1)],[env.supervoxel(:,2) env.supervoxel(:,1) env.supervoxel(:,3)],env.vol)
if isempty(numK)
    numK=3;
end
para=struct;
para=setfield(para,'k1',numK);
para=setfield(para,'merge',0.7);
para=setfield(para,'cap',0.7);
para=setfield(para,'reg1',0.7);para=setfield(para,'reg2',0.7);
para=setfield(para,'minSize',10);

ind=ind_increased.pagetst_acq_nonadj(1:17000);
[cIX2,gIX2] = AutoClustering(ind,ones(size(ind)),M_0,cIX_reg,1,para,1,0.7);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
[gIX2, numU] = hierplot_zyq_20190530(cIX2,gIX2,M_0(cIX2,:));
clrmap=pushbutton_popupplot_Callback(activities_preCS_dfdf_aftcorrect',cIX2,gIX2,env,fs.ca,stimCS,stimUS,0,6,0,0,0,0);
DrawTiledPics_zyq_20190530(cIX2,gIX2,[1:size(M_0,1)],[env.supervoxel(:,2) env.supervoxel(:,1) env.supervoxel(:,3)],env.vol)
% %tSNE
% tic;[mappedA, mapping] = compute_mapping(M_0(cIX2,:),'tSNE');toc
% mappedA = tsne(M_0(cIX2,:), [], 2); % no_dims = 2
% figure;
% n = round(numU*1.1);
% cmap = hsv(max(1,n));
% gscatter(mappedA(:,1),mappedA(:,2),gIX2,cmap,'.',[],'off');



figure,
B=[area.CS_hab_tst(1,:); area.CS_acq_block ;area.CS_hab_tst(2,:)];
x=1:size(B,1);
a=[-0.01 0.2];
for i = 1:length(unique(gIX2))
    subplot(ceil(length(unique(gIX2))/2),2,i)
    clust = find(gIX2==i);
%     plot(x,B(:,cIX2(clust)),'color',[0.5 0.5 0.5],...
%         'LineWidth',1,...
%         'markersize',3);
%     hold on
    plot(x,mean(B(:,cIX2(clust)),2),'-o','color',clrmap(i,:),...
        'LineWidth',2,...
        'markersize',3);
    hold on
    set(gca,'FontSize',20,'xtick',x,'TickLength',[0.01 0.01]);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
    text(5,(min(a)+max(a))/2,num2str(length(clust)),'fontsize',25);hold on;
    ylabel('Integrated ¡÷F/F0');xlabel('Blocks');
    ylim(a)
end
hold off

B=rawacti.acq_mean_blocks;
x=1:size(B{1},1);
a=[-0.01 0.03];
for i = 1:length(unique(gIX2))
    figure,
    clust = find(gIX2==i);  
    BB=[];
    for jj=1:length(clust)
        BB(:,:,jj)=B{clust(jj)};
    end
    plot(x,mean(BB,3),'-',...
        'LineWidth',2,...
        'markersize',3);
    hold on
    set(gca,'FontSize',12);box off;xlim([0.5 length(x)+0.5]);%ylim([0 1.5]);
    text(5,(min(a)+max(a))/2,num2str(length(clust)),'fontsize',25);hold on;
    ylabel('¡÷F/F0');xlabel('Time(s)');
    ylim(a)
end
hold off
