set(0,'defaultfigurecolor','w');
%%
% javaaddpath('F:\DUlab\Matlab\snn\lib\javaml\javaml\0.1.7\javaml-0.1.7.jar');
% javaaddpath('F:\DUlab\Matlab\snn\lib\ajt\ajt\2.9\ajt-2.9.jar');
% %javaaddpath('F:\DUlab\Matlab\snn\lib\javaml\javaml\0.1.7\lib\commons-math-1.2.jar');
% javaaddpath('F:\DUlab\Matlab\snn\lib\jama\jama\1.0.2\Jama-1.0.2.jar');
% javaaddpath('F:\DUlab\Matlab\snn\lib\libsvm\libsvm\0.0\libsvm-0.0.jar');
% javaaddpath('F:\DUlab\Matlab\snn\lib\weka\weka\0.0\weka-0.0.jar');
%% load
outpath='G:\data_huc_non_learner\20190605\fish1\';
input=load('G:\data_huc_non_learner\20190605\fish1\activities_aft_process.mat');
load('G:\data_huc_non_learner\20190605\fish1\para.mat');
load('G:\data_huc_non_learner\20190605\fish1\env.mat');
load('G:\data_huc_non_learner\20190605\fish1\region_mask.mat');
X=env.supervoxel;
for ii=unique(X(:,3))'
    id=find(X(:,3)==ii);
    X(id,3)=(ii-1)*8/0.66+1;
end
env.supervoxel=X;
spon=load('E:\A_Data_lightsheet\Data_huc\20190905\fish2\spon1\act.mat');
spon2=load('E:\A_Data_lightsheet\Data_huc\20190905\fish2\spon2\act.mat');
winsize=121;
base=compute_dfdf(spon.activities,winsize);base=base';
base=base((winsize-1)/2+1:end-(winsize-1)/2,:);
base_r=base(1:floor(size(base,1)/trial.total)*trial.total,:);
base_r=reshape(base_r,floor(size(base,1)/trial.total),trial.total,[]);

base2=compute_dfdf(spon2.activities,winsize);base2=base2';
base2=base2((winsize-1)/2+1:end-(winsize-1)/2,:);

%% PLOT PARA
color_cs=[0 0.45 0.74];
fonts_axis=16;
%% behavior
tail_movement_by_trial_fin_cor(y_3sd,re_startpoint,[],[],[]);

figure,
yli=[0 trial.total+2];
y=trial.total-re_startpoint(:,1)+1;
patch1=patch([frameb.cs_start,...
    frameb.cs_end,...
    frameb.cs_end,...
    frameb.cs_start]/fs.behavior,...
    [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
    color_cs,'edgecolor',color_cs,'facecolor',color_cs,'edgealpha',0.2,'facealpha',0.25);hold on
l3=plot([frameb.us_start frameb.us_start]/fs.behavior,trial.total-[trial.acq(2) trial.acq(3)]+1,'r','LineStyle','-','linewidth',1);hold on
scatter(re_startpoint(:,2)/fs.behavior,y,14,'k','s','filled');hold on
set(gca,'fontsize',fonts_axis);%,'ytick',[trial.total:-5:0]
axis equal;xlim([0 time.per_cycle]);ylim([0 trial.total+2]);xlabel('Time(s)');ylabel('Trial number');

figure,n=fs.behavior*0.2;
for t=1:5
    x=[];
    switch t
        case 1
            bin=trial.hab(2):trial.hab(3);
        case 2
            bin=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1;
        case 3
            bin=trial.acq(2)+trial.acq_block_trial:trial.acq(2)+trial.acq_block_trial*4-1;
        case 4
            bin=trial.acq(2)+trial.acq_block_trial*4:trial.acq(3);
        case 5
            bin=trial.test(2):trial.test(3);
    end
    for ii=1:frameb.per_cycle-n
        a=intersect(bin,re_startpoint(find(re_startpoint(:,2)>ii & re_startpoint(:,2)<ii+n),1));
        x(ii)=length(a)/(20/fs.behavior*length(bin));
    end
    plot((1:length(x))/fs.behavior,x,'linewidth',2);hold on;
end
legend({'hab','early','middle','late','test','CS','US'});
yli=[0 3.5];
patch1=patch([frameb.cs_start,...
    frameb.cs_end,...
    frameb.cs_end,...
    frameb.cs_start]/fs.behavior,...
    [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
    color_cs,'edgecolor',color_cs,'facecolor',color_cs,'edgealpha',0.2,'facealpha',0.25);hold on
l3=plot([frameb.us_start frameb.us_start]/fs.behavior,[min(yli(:)) max(yli(:))],'r','LineStyle','-','linewidth',1);hold on
xlim([0 time.per_cycle]);ylim([min(yli(:)) max(yli(:))]);xlabel('Time(s)');ylabel('Shakes/s');set(gca,'fontsize',fonts_axis);
%% task moulate
% Task-modulated cells were identified based on a Wilcoxon rank-sum test of
% average firing rate within the 1 s baseline and three sequential 1 s epochs within the task
% starting at the odor onset, on a per trial basis. A cell was deemed task modulated if any of
% the task epochs differed significantly from the baseline after false discovery rate (FDR)
% correction across all units, with a corrected significance threshold of P < 0.001.
A=input.activities_preCS_dfdf_aftcorrect;
% A_r=zeros(frame.per_cycle,trial.total,size(A,2));
% for ii=1:size(A,2)
%     for jj=1:trial.total
%     A_r(:,jj,ii)=A((jj-1)*frame.per_cycle+1:jj*frame.per_cycle,ii)';
%     end
% end
A_r=reshape(A,frame.per_cycle,trial.total,[]);figure;plot(A_r(:,30,100))
baseline=mean(A_r(1:frame.cs_start-2,:,:),1);
%baseline=mean(base_r,1);
B=zeros(2,trial.total,size(A,2));
for ii=1:2
    switch ii
        case 1
            bin=[1:frame.cs_start-1];
        case 2
            bin=[frame.cs_start:frame.cs_end];
    end
    B(ii,:,:)=mean(A_r(bin,:,:),1);
end
%B_r=reshape(B,ii*trial.total,[]);
p=nan(2,size(A,2));
for ii=1:size(A,2)
    for jj=1:2
        p(jj,ii) = ranksum(baseline(:,:,ii),B(jj,:,ii));
    end
end
%p_r=reshape(p,1,[]);
[FDR,p_adj] = mafdr(p(:));
p_r=false(2,size(A,2));
p_r(find(p_adj<0.001))=true;
task_moulate_ind= find((p_r(1,:) | p_r(2,:))==1);
task_non_moulate_ind=find((p_r(1,:) | p_r(2,:))~=1);
figure,imagesc(A(:,task_non_moulate_ind)',[0 0.5]);colorbar;
figure,
subplot(2,1,1),shadedErrorBar([1:size(base,1)],base(:,task_non_moulate_ind)',{@mean,@(x) std(x)},'lineprops','k','transparent',1,'patchSaturation',0.1); hold on
subplot(2,1,2),shadedErrorBar([1:size(A,1)],A(:,task_non_moulate_ind)',{@mean,@(x) std(x)},'lineprops','k','transparent',1,'patchSaturation',0.1); hold on
figure,
scatter3(env.supervoxel(task_non_moulate_ind,1),env.supervoxel(task_non_moulate_ind,2),env.supervoxel(task_non_moulate_ind,3),12,'k','filled');hold on;
scatter3(env.supervoxel(task_moulate_ind,1),env.supervoxel(task_moulate_ind,2),env.supervoxel(task_moulate_ind,3),12,'r','filled');axis equal;
legend('task non moulate','task moulate');

task_moulate_ind=1:size(A,2);
task_moulate_ind(find(env.supervoxel(task_moulate_ind,3)==1))=[];
frame_ind=1:frame.per_cycle;
task_moulate_act=A_r(frame_ind,:,task_moulate_ind);
%% conditions & compute trial_avg
%state:early/middle/late learnining stage
%cue:CS/CS-US
%behavior:corr/incorr
S=3;%state {'baseline','middle trian','learned'}
D=2;%behavior {'corr','incorr'}
C=2;%cue {'CS','US'}
go_trial_cs=unique(re_startpoint(find(re_startpoint(:,2)>=frameb.cs_start & re_startpoint(:,2)< frameb.us_start ),1))';
ingo_trial_cs=setdiff(1:trial.total,go_trial_cs);
corr_trial_hab=intersect(ingo_trial_cs,[trial.hab(2):trial.hab(3)])';%CS & non mov in hab
corr_trial_tst=intersect(go_trial_cs,[trial.test(2):trial.test(3)])';%CS & mov in test
corr_trial=[corr_trial_hab;corr_trial_tst];

incorr_trial_hab=intersect(go_trial_cs,[trial.hab(2):trial.hab(3)])';%CS & mov in hab
incorr_trial_tst=intersect(ingo_trial_cs,[trial.test(2):trial.test(3)])';%CS & non mov in test
incorr_trial=[incorr_trial_hab;incorr_trial_tst];

corr_trial=[corr_trial;nan(trial.test(1)+trial.hab(1)-length(corr_trial),1)];
incorr_trial=[incorr_trial;nan(trial.test(1)+trial.hab(1)-length(incorr_trial),1)];
conditions_ind.behav=[corr_trial,incorr_trial];

% conditions.state=[trial.acq(2):trial.acq_block_trial+trial.acq(2)-1;...
%                   trial.acq_block_trial*2+trial.acq(2):trial.acq_block_trial*3+trial.acq(2)-1;
%                   trial.acq_block_trial*4+trial.acq(2):trial.acq(3)];
%conditions_ind.state=[[trial.hab(2):trial.hab(3),nan(1,trial.test(1)-trial.hab(1))];trial.test(2):trial.test(3)]';% naive/learned
n=18;conditions_ind.state=[[trial.hab(2):trial.hab(3),nan(1,n-trial.hab(1))];...;
    trial.acq(2)+trial.acq_block_trial:trial.acq(2)+trial.acq_block_trial*4-1;...
    [trial.test(2):trial.test(3),nan(1,n-trial.test(1))]]';% naive/learned
n=trial.acq(1)-trial.hab(1)-trial.test(1);
conditions_ind.cue=[[trial.hab(2):trial.hab(3),trial.test(2):trial.test(3),nan(1,n)];trial.acq(2):trial.acq(3)]';
%trial average across different conditions
% % trail_avg=[];trail_avg.state=[];
% % for ii=1:size(conditions.state,2)
% %     a=conditions.state(:,ii);a(isnan(a))=[];
% %     trail_avg.state(:,ii,:)=squeeze(mean(task_moulate_act(:,a,:),2));
% % end
% % trail_avg.behav=[];
% % for ii=1:size(conditions.behav,2)
% %     a=conditions.behav(:,ii);a(isnan(a))=[];
% %     trail_avg.behav(:,ii,:)=squeeze(mean(task_moulate_act(:,a,:),2));
% % end
trail_avg=zeros(size(task_moulate_act,3),S,D,size(task_moulate_act,1));
for ii=1:S %{'baseline','middle trian','learned'}
    for jj=1:D %{'corr','incorr'}
        %a=intersect(conditions_ind.behav(:,jj),conditions_ind.state(:,ii));a(isnan(a))=[];
        switch jj
            case 1
                if ii==1 x=ingo_trial_cs;else x=go_trial_cs;end
            case 2
                if ii==1 x=go_trial_cs;else x=ingo_trial_cs;end
        end
        a=intersect(conditions_ind.state(:,ii),x);b=mean(task_moulate_act(:,a,:),2);
        b=reshape(b,[],1);b=normalize(b,'zscore');b=reshape(b,size(task_moulate_act,1),1,[]);
        trail_avg(:,ii,jj,:)=[squeeze(b)]';
        %         figure,plot(b(:,100));a=trail_avg(100,ii,jj,:);
        %         figure,plot( a(:));
    end
end
%% PCA
projected_data=[];coeff=[];explained_var=[];
% X=cat(2,trail_avg.state,trail_avg.behav);
% for ii=1:size(X,2)
% [projected_data(ii,:,:,:),~,coeff(:,ii,:),explained_var(ii,:)]=mypca(X(:,ii,:),6,0);
% end
kk=1;
for ii=1:S
    for jj=1:D
        [projected_data(kk,:,:,:),~,coeff(:,kk,:),~]=mypca(reshape(squeeze(trail_avg(:,ii,jj,:))',size(trail_avg,4),1,[]),20,0,frame);
        kk=kk+1;
    end
end
for ii=1:C
    a=conditions_ind.cue(:,ii);a(isnan(a))=[];b=mean(task_moulate_act(:,a,:),2);
    b=reshape(b,[],1);b=normalize(b,'zscore');b=reshape(b,size(task_moulate_act,1),1,[]);
    [projected_data(kk,:,:,:),~,coeff(:,kk,:),~]=mypca(b,20,0,frame);
    kk=kk+1;
end
coeff(:,isnan(projected_data(:,1,1,1)),:)=[];projected_data(isnan(projected_data(:,1,1,1)),:,:,:)=[];
projected_data=squeeze(projected_data);

for jj=1:size(coeff,2)
    figure,
    for ii=1:10
        subplot(1,10,ii);
        barh(coeff(:,jj,ii))
        xlim([-0.005 0.005])
    end
end
%% dPCA
%% CLUSTER BY PC
X=coeff(:,:,1:20);
XX=[];
for ii=1:size(X,2)
    XX=cat(3,XX,X(:,ii,:));
end
XX=squeeze(XX);
%[a,~,numK,~]=find_best_k_in_range(XX,2:20);
[idx_all,~,~,~] = kmeans(XX,6,'Distance','sqeuclidean','Replicates',10,'Display','final','MaxIter',1000);
%% baseline SD VS state
for ii=1:S
    a=task_moulate_act(1:frame.cs_start-1,conditions_ind.state(:,S),:);
    a=reshape(a,conditions_ind.state(:,S)*(frame.cs_start-1),[]);
    state.SD(:,ii)=std(a,[],2);
end
figure;plot(state.SD);

%% regressor
%CS: HAB/TRIANING/TEST
%MOV:SPON/INDUCED/
%STATE:
%US
colorCS=[0.5 0.5 0.9];
prct_const=2;nCells_total = size(A,2);topN = round(prct_const/100 * nCells_total);
reg_thres=0.4;reg_thres2=0.8;%para
stimCS=zeros(1,frame.per_cycle*trial.acq_block_trial);
for ii=1:size(stimCS,2)/frame.per_cycle
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=1;%23
end
stim=zeros(5,trial.acq_block_trial*frame.per_cycle);stim(1,:)=stimCS;
[stregressors,~,~,~] = GetMotorRegressor(stim,1);
regress.CS_acq_block=stregressors(1,1).im;
stimcorr=zeros(size(task_moulate_act,3),trial.acq_block_num+2);
stAvrcorr=zeros(size(task_moulate_act,3),trial.acq_block_num+2);
task_moulate_act_acq_blockmean=zeros(frame.per_cycle,trial.acq_block_num+2,size(task_moulate_act,3));
for ii=1:trial.acq_block_num+2
    if ii>=2 && ii<=6
        frame_ind=1:frame.us_start-1;isus=true;
    else
        frame_ind=1:frame.per_cycle;isus=false;
    end
    switch ii
        case 1
            trial_ind=trial.hab(2):min(trial.hab(3),trial.hab(2)+trial.acq_block_trial-1);
        case 2
            trial_ind=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1;
        case 3
            trial_ind=trial.acq(2)+trial.acq_block_trial*1:trial.acq(2)+trial.acq_block_trial*2-1;
        case 4
            trial_ind=trial.acq(2)+trial.acq_block_trial*2:trial.acq(2)+trial.acq_block_trial*3-1;
        case 5
            trial_ind=trial.acq(2)+trial.acq_block_trial*3:trial.acq(2)+trial.acq_block_trial*4-1;
        case 6
            trial_ind=trial.acq(2)+trial.acq_block_trial*4:trial.acq(3);
        case 7
            trial_ind=trial.test(2):min(trial.test(3),trial.test(2)+trial.acq_block_trial-1);
    end
    a=task_moulate_act(:,trial_ind,:);a=mean(a,2);task_moulate_act_acq_blockmean(:,ii,:)=a;
    M=task_moulate_act(frame_ind,trial_ind,:);M=reshape(M,size(M,1)*size(M,2),[])';
    stim_output=[regress.CS_acq_block];stim_output=reshape(stim_output,[],length(trial_ind));stim_output=stim_output(frame_ind,:);stim_output=reshape(stim_output,1,[]);
    [stimcorr(:,ii),~] = MotorSourceCorrelation(M,stim_output,[]);
    [~,IX] = sort(stimcorr(:,ii),'descend');
    thr = stimcorr(IX(topN),ii);
    ind_output = (find(stimcorr(:,ii)>max(reg_thres,thr)))';
    stim_output=mean(M(ind_output,:),1);
    [stim_tAvr,~,~,~,~] = GetTrialAvrLongTrace_zyq_20190730(stim_output,length(frame_ind));
    stAvrcorr(:,ii) = corr(stim_tAvr',M');
    [~,IX] = sort(stAvrcorr(:,ii),'descend');
    x0 = stAvrcorr(IX(topN),ii);
    thr=max(x0,reg_thres2);IX_passX = find(stAvrcorr(:,ii)>=thr);
    figure('name',['regressor vs comt temp ' num2str(ii)]),
    histogram(stimcorr(:,ii),'BinEdges',-1:0.01:1,'FaceColor','b');hold on
    histogram(stAvrcorr(:,ii),'BinEdges',-1:0.01:1,'FaceColor','r');hold on
    figure('name',num2str(ii)),
    x=stim_output';a=[-0.05 0.4];
    patch1=patch([[frame.cs_start:length(frame_ind):size(x,1)]'...
        [min(frame.cs_end,frame_ind(end)):length(frame_ind):size(x,1)]'...
        [min(frame.cs_end,frame_ind(end)):length(frame_ind):size(x,1)]'...
        [frame.cs_start:length(frame_ind):size(x,1)]']',...
        repmat([min(a) min(a) max(a) max(a)],length([frame.cs_start:length(frame_ind):size(x,1)]),1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    p1=plot(reshape(stim_output,1,[]),'k','linewidth',1);hold on;
    p2=plot(stim_tAvr,'b','linewidth',1);hold on;
    %     if isus
    %         l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
    %     end
    p3=plot(mean(M(IX_passX,:),1),'r','linewidth',1.5);hold on;
    if ii==7 legend([patch1 p1,p2,p3],{'CS','Regressor','Comp. Templete','>thr Avg. trace'},'location','northeastoutside');end
    ylim(a);xlim([1 length(x)]);
    box off;
    %line([thr thr],[0 1500],'linewidth',1.5,'color','r','linestyle','--');hold on
end
%cluster
%figure('name','hireachy_tree_regress profile');
stAvrcorr_n=normalize(stAvrcorr,1);
[gIX_CSregressor,~,~,~] = kmeans(stAvrcorr_n,5,'Distance','correlation','Replicates',10,'Display','final','MaxIter',1000);
cmap= GetColormap('hsv_new',max(gIX_CSregressor));
[h,~]=plot_test_2(task_moulate_act,1:length(task_moulate_ind),gIX_CSregressor,frame,cmap,2,[0 0.05],true,true);
figure;scatter3(env.supervoxel(task_moulate_ind,2),env.supervoxel(task_moulate_ind,1),env.supervoxel(task_moulate_ind,3),12,cmap(gIX_CSregressor,:),'filled');axis equal;
[~,Ix_s_ind]=sort(gIX_CSregressor);
figure('name','ind_clu_s'),imagesc(stAvrcorr_n(Ix_s_ind,:),[-1 1]);colorbar;%colormap(cmap);%colorbar('Ticks',[1:max(gIX_CSregressor)]);
%% getmask
e=load('E:\A_Data_lightsheet\Data_huc\20190905\fish2\training\env.mat');
% %correct
% mask=zeros(size(env.vol(:,:,1:24)));
% for ii=1:24
% mask(:,:,ii)=getmask_imfreehand(env.vol(:,:,ii),6);
% end
% mask(:,2047:end,:)=[];
% %show_spv_GUI(mask)
% show_spv_GUI(reg_mask);
% reg_mask_c=reg_mask;kk=max(reg_mask_c(:));
% for rr=[4,11]
%     for jj=1:6
%         for ii=1:24
%
%             a=squeeze([mask(:,:,ii)]);a=flip(a,1);
%             id=find(a'==jj);
%             id_r=find(squeeze([reg_mask(:,:,ii)])==4);
%             idd=intersect(id,id_r);
%             %         a=zeros(2046,1348);a(id)=1;a(id_r)=2;
%             %         figure,imshow(a,[0 5]);
%             a=reg_mask_c(:,:,ii);a(idd)=kk+1;
%             reg_mask_c(:,:,ii)=a;
%         end
%         kk=kk+1;
%     end
% end
% show_spv_GUI(reg_mask_c);
%[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,id_in_region_cell]=get_region_fraction(reg_mask,reg_name,reg_loc',gIX_CSregressor,e.env.supervoxel(task_moulate_ind,1:3),cmap);
%% GC
mean_act2=[];mean_act=[];mean_act_training=[];GC_nonzero_train={};GCstrength_train={};GCstrength_spon={};
id=[1,8,2,9,3,10,4,11,5,12,6,13,7];
I=[1:length(id_in_region_cell);1:length(id_in_region_cell)]';
for ii=1:length(id_in_region_cell)
    for jj=1:trial.acq_block_num
        ind=trial.hab(3)+trial.acq_block_trial*(jj-1)+1:trial.hab(3)+trial.acq_block_trial*jj;
        mean_act_training(:,jj,ii)=reshape(mean(task_moulate_act(1:frame.cs_start-1,ind,id_in_region_cell{ii}),3),[],1);
    end
    mean_act(:,ii)=mean(base(:,task_moulate_ind(id_in_region_cell{ii})),2);
    mean_act2(:,ii)=mean(base2(:,task_moulate_ind(id_in_region_cell{ii})),2);
end
for ii=1:trial.acq_block_num
    [GC_nonzero_train{ii},GCstrength_train{ii}]=Groupar_mls_concatenated(squeeze(mean_act_training(:,ii,:)),frame.cs_start-1);
    for jj=1:length(id_in_region_cell)
        GCstrength_train{ii}(jj,jj)=inf;
    end
end
[GC_nonzero_spon{1},GCstrength_spon{1}]=Groupar_mls_concatenated(mean_act,size(mean_act,1),1);
[GC_nonzero_spon{2},GCstrength_spon{2}]=Groupar_mls_concatenated(mean_act2,size(mean_act2,1),1);
for ii=1:2
    for jj=1:length(id_in_region_cell)
        GCstrength_spon{ii}(jj,jj)=inf;
    end
end
figure('name',['GC_strength spon']);set(gcf,'position',[50,200,1600,400]);
for ii=1:2
    subplot(1,2,ii),imagesc(GCstrength_spon{ii}(id,id));
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end
figure('name',['GC_strength training']);set(gcf,'position',[50,200,1600,400]);
for ii=1:trial.acq_block_num
    subplot(1,trial.acq_block_num,ii),imagesc(GCstrength_train{ii}(id,id));title(['Acq. ' num2str(ii)]);
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end
%increase
figure('name',['GC_strength spon increase']);set(gcf,'position',[50,200,1600,400]);
for ii=1:2
    subplot(1,2,ii),imagesc([(GCstrength_spon{ii}(id,id)-GCstrength_spon{1}(id,id)) > 0]);
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end
figure('name',['GC_strength training increase']);set(gcf,'position',[50,200,1600,400]);
for ii=1:trial.acq_block_num
    subplot(1,trial.acq_block_num,ii),imagesc([GCstrength_train{ii}(id,id)-GCstrength_train{1}(id,id)] > 0);title(['Acq. ' num2str(ii)]);
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end
%decrease
figure('name',['GC_strength spon increase']);set(gcf,'position',[50,200,1600,400]);
for ii=1:2
    subplot(1,2,ii),imagesc([(GCstrength_spon{ii}(id,id)-GCstrength_spon{1}(id,id)) < 0]);
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end
figure('name',['GC_strength training increase']);set(gcf,'position',[50,200,1600,400]);
for ii=1:trial.acq_block_num
    subplot(1,trial.acq_block_num,ii),imagesc([GCstrength_train{ii}(id,id)-GCstrength_train{1}(id,id)] < 0);title(['Acq. ' num2str(ii)]);
    set(gca,'xticklabels',{reg_name{id}},'xtick',1:length(id_in_region_cell),'ytick',1:length(id_in_region_cell),'XTickLabelRotation',45,'yticklabels',{reg_name{id}},'yTickLabelRotation',45,'fontsize',12);
end

%% save
save([outpath 'diff_component.mat'],'task_moulate_ind','conditions_ind','trail_avg','projected_data','coeff','idx_all','stimcorr','stAvrcorr');
%% visualization
cmap= GetColormap('hsv_new',max(idx_all));
%PC
figure,a=3;n=6;size(projected_data,3);
clrmap= GetColormap('hsv_new',size(projected_data,1));hsv(max(1,size(projected_data,1)));
for ii=1:n;%
    yli=[min(min(projected_data(:,:,ii))) max(max(projected_data(:,:,ii)))];
    subplot(ceil(n/a),a,ii);%set(gca,'colororder',GetColormap('hsv_new',n))
    for jj=1:size(projected_data,1)
    plot(squeeze(projected_data(jj,:,ii))','color',clrmap(jj,:),'linewidth',1.5);hold on;
    end
    patch1=patch([frame.cs_start,...
        frame.cs_end,...
        frame.cs_end,...
        frame.cs_start],...
        [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
        color_cs,'edgecolor',color_cs,'facecolor',color_cs,'edgealpha',0.2,'facealpha',0.25);hold on
    l3=plot([frame.us_start frame.us_start],[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',3);hold on
    title(['PC ' num2str(ii)]); xlim([1 frame.per_cycle]);xlabel('frame');
    %set(gca,'fontsize',fonts_axis);
end
legend({'Baseline','Training&shake','Training&nonshake','Test&shake','Test&nonshake','CS','US'});xlabel('frame');
%cluster
clrmap=[1 0.2 0.2; 0.2 0.2 1 ];yli=[-0.5 0.5];
for id=unique(idx_all)'
    figure;kk=1;
    for ii=1:S
        subplot(S,1,kk);
        for jj=1:D
            a=squeeze(trail_avg((find(idx_all==id)),ii,jj,:));
            %shadedErrorBar(1:size(a,2),a,{@(x) mean(x),@(x) std(x,[],1)/sqrt(length(x))},'lineprops',{clrmap(jj,:,:)},'transparent',1,'patchSaturation',0.1);hold on
            plot(mean(a),'linewidth',2,'color',clrmap(jj,:,:));hold on;
        end
        kk=kk+1;
        patch1=patch([frame.cs_start,...
            frame.cs_end,...
            frame.cs_end,...
            frame.cs_start],...
            [min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
            color_cs,'edgecolor',color_cs,'facecolor',color_cs,'edgealpha',0.2,'facealpha',0.25);hold on
        l3=plot([frame.us_start frame.us_start],[min(yli(:)) max(yli(:))],'r','LineStyle','--','linewidth',2);hold on
        xlim([1 frame.per_cycle]);
        switch ii
            case 1
                title(['Baseline ' num2str(id)]);
            case 2
                title(['Learning ' num2str(id)]);
            case 3
                title(['Learned ' num2str(id)]);
        end
    end
    legend('corr','incorr','CS','US');
end
% typical case of each cluster
for id=unique(idx_all)'
    a=find(idx_all==id);
    switch id
        case 1
            a1=squeeze(mean(trail_avg(a,2,1,frame.cs_start:frame.us_start-1),4));
            a2=squeeze(mean(trail_avg(a,2,2,frame.cs_start:frame.us_start-1),4));
            [~,ii]=max(a1-a2);
        case 2
            a1=squeeze(mean(mean(trail_avg(a,2,1:2,frame.cs_start:frame.us_start-1),4),3));
            a2=squeeze(mean(mean(trail_avg(a,3,1:2,frame.cs_start:frame.us_start-1),4),3));
            [~,ii]=max(a2-a1);
        case 3
            a1=squeeze(mean(mean(trail_avg(a,2,1:2,frame.cs_start:frame.us_start-1),4),3));
            a2=squeeze(mean(mean(trail_avg(a,3,1:2,frame.cs_start:frame.us_start-1),4),3));
            a3=squeeze(mean(mean(trail_avg(a,1,1,frame.cs_start:frame.us_start-1),4),3));
            [~,ii]=min(sum(abs(diff([a3,a2,a1],1,2)),2));
        case 4
            ii=1000;
    end
    a=a(ii);
    figure;
    plot(squeeze(trail_avg(a,1,1,:)),'k','linewidth',2);hold on
    plot(squeeze(trail_avg(a,2,1,:)),'r','linewidth',2);hold on;
    plot(squeeze(trail_avg(a,2,2,:)),'r--','linewidth',2);hold on;
    plot(squeeze(trail_avg(a,3,1,:)),'b','linewidth',2);hold on;
    plot(squeeze(trail_avg(a,3,2,:)),'b--','linewidth',2);hold on;
    legend({'Baseline','Learning & Go','Learning & NoGo','Learned & Go','Learned & NoGo'});
end
%mapback
%gscatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),idx_all,cmap,'.',[],'on','t1','t2');axis equal;
%scatter3(env.supervoxel(task_moulate_ind,1),env.supervoxel(task_moulate_ind,2),env.supervoxel(task_moulate_ind,3),12,cmap(idx_all,:),'filled');hold on;
kk=1;
id_plot=1:length(idx_all);id_plot(find(idx_all== 7 | idx_all== 9 | idx_all== 10 ))=[];
idx_all_plot=idx_all;idx_all_plot(find(idx_all==2))=5;
cmap= GetColormap('hsv_new',7);max(idx_all_plot(id_plot))+1
for ii=1:env.depth
    if mod(ii,6)==1
        figure('position',[50,50,1800,550]),kk=1;
    end
    if ii==env.depth
        a='on';
    else
        a='off';
    end
    subplot(1,6,kk);
    id=find(env.supervoxel(task_moulate_ind(id_plot),3)==((ii-1)*8/0.66+1));
    id_n=find(env.supervoxel(task_non_moulate_ind,3)==((ii-1)*8/0.66+1));
    gscatter(env.supervoxel(task_non_moulate_ind(id_n),2),env.supervoxel(task_non_moulate_ind(id_n),1),ones(1,length(id_n)),[0.5 0.5 0.5],'.',[],'off');hold on;
    gscatter(env.supervoxel(task_moulate_ind(id_plot(id)),2),env.supervoxel(task_moulate_ind(id_plot(id)),1),idx_all_plot(id_plot(id)),cmap,'.',[],a);hold on;
    %     scatter(env.supervoxel(task_non_moulate_ind(id_n),1),env.supervoxel(task_non_moulate_ind(id_n),2),12,[0.5 0.5 0.5],'filled');hold on;
    %     scatter(env.supervoxel(task_moulate_ind(id),1),env.supervoxel(task_moulate_ind(id),2),12,cmap(idx_all(id),:),'filled');
    xlim([0 1400]);axis equal;
    set(gca,'visible','off');
    kk=kk+1;
end
[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,id_in_region_cell]=get_region_fraction(reg_mask,reg_name,reg_loc',idx_all_plot(id_plot),e.env.supervoxel(task_moulate_ind(id_plot),1:3),cmap);
%check & sort by motor
e=load('E:\A_Data_lightsheet\Data_huc\20190905\fish2\training\env.mat');
for id=[1,2,3,4,5,6,8];
    a=task_moulate_ind(find(idx_all==id));%a=a(1);
    %showspv= mapback(ones(length(a),1), e.env.supervoxel(a,:),[e.env.height e.env.width e.env.depth],1:length(a));
    %show_spv_GUI(showspv);
    %figure,plot_rawtrace_trials(A(:,a),[],fs,frame,trial,startpoint,1);
    b=zeros(16,length(startpoint),length(a));
    for ii=1:length(startpoint)
        t=max(round(startpoint(ii)/fs.behavior/fs.ca)-5,1):min(round(startpoint(ii)/fs.behavior/fs.ca)+10,size(A,1));
        b(1:length(t),ii,:)=A(t,a);
    end
    bb=squeeze(mean(b,2));
    figure,
    line([6 6],[-0.02 0.1]);
    plot(bb,'color',[0.5 0.5 0.5]);
    %shadedErrorBar([1:16],bb',{@mean,@(x) std(x)},'lineprops','k','transparent',1,'patchSaturation',0.1); hold on
end
