%% Single cell analysis of state,valence,behavior
%zyq,202505
clc;clear all;close all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
for batchi=1:4
    for fishi=1:length(Path{batchi})
        %load single fish data
        p=char(Path{batchi}(fishi));
        nn=[p(end-14:end-7),p(end-5:end-1)];
        load(fullfile(p,'activities_dfdf_align.mat'),'align_win');
        load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(p,'para.mat'));
        A=zscore(activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:),[],'all');
        A_r=reshape(A,frame.per_cycle,trial.total,[]);
        go_trial={};nogo_trial={};trial_ind_session={};
        CR=re_startpoint(find(re_startpoint(:,2)<=frameb.cs_end-0.5*fs.behavior & re_startpoint(:,2)>=frameb.cs_start),:);
        for ii=1:size(align_win,2)
            actm1=align_win(:,ii);actm1(isnan(actm1))=[];trial_ind_session{ii}=actm1;
            go_trial{ii}=unique(CR(find(CR(:,1)<=trial_ind_session{ii}(end) & CR(:,1)>=trial_ind_session{ii}(1)),1));
            nogo_trial{ii}=setdiff(trial_ind_session{ii},go_trial{ii});
        end
        % set period para
        trial_ind=struct;period_ind=struct;Difference=struct;
        trial_ind.state{1}{1}=trial.hab(2):trial.hab(3);trial_ind.state{1}{2}=trial.test(2):trial.test(3);
        period_ind.state{1}=frame.us_start+4/fs.ca:frame.per_cycle;
        trial_ind.cue{1}{1}=trial.hab(2):trial.hab(3);trial_ind.cue{1}{2}=trial.test(2):trial.test(3);
        trial_ind.cue{2}{1}=intersect(trial.hab(2):trial.hab(3),nogo_trial{1});trial_ind.cue{2}{2}=intersect(trial.test(2):trial.test(3),nogo_trial{end});
        period_ind.cue{1}=frame.cs_start:frame.us_start-1;
        trial_ind.behav{1}{1}=intersect(trial.test(2):trial.test(3),go_trial{end});trial_ind.behav{1}{2}=intersect(trial.hab(2):trial.hab(3),nogo_trial{end});
        period_ind.behav{1}=frame.cs_start:frame.us_start-1;
        % average activity acorss period
        actm1=nan(size(A_r,3),length(trial_ind.state{1}{2}));actm1(:,:)=squeeze(nanmean(A_r(period_ind.state{1},trial_ind.state{1}{2},:),1))';
        actm2=nan(size(A_r,3),length(trial_ind.state{1}{1}));actm2(:,:)=squeeze(nanmean(A_r(period_ind.state{1},trial_ind.state{1}{1},:),1))';
        Difference.state{1,1}=actm1;Difference.state{1,2}=actm2;
        actm1=nan(size(A_r,3),length(trial_ind.cue{1}{2}));actm1(:,:)=squeeze(nanmean(A_r(period_ind.cue{1},trial_ind.cue{1}{2},:),1))';
        actm2=nan(size(A_r,3),length(trial_ind.cue{1}{1}));actm2(:,:)=squeeze(nanmean(A_r(period_ind.cue{1},trial_ind.cue{1}{1},:),1))';
        Difference.cue{1,1}=actm1;Difference.cue{1,2}=actm2;
        actm1=nan(size(A_r,3),length(trial_ind.cue{2}{2}));actm1(:,:)=squeeze(nanmean(A_r(period_ind.cue{1},trial_ind.cue{2}{2},:),1))';
        actm2=nan(size(A_r,3),length(trial_ind.cue{2}{1}));actm2(:,:)=squeeze(nanmean(A_r(period_ind.cue{1},trial_ind.cue{2}{1},:),1))';
        Difference.cue{2,1}=actm1;Difference.cue{2,2}=actm2;
        actm1=nan(size(A_r,3),length(trial_ind.behav{1}{2}));actm1(:,:)=squeeze(nanmean(A_r(period_ind.behav{1},trial_ind.behav{1}{2},:),1))';
        actm2=nan(size(A_r,3),length(trial_ind.behav{1}{1}));actm2(:,:)=squeeze(nanmean(A_r(period_ind.behav{1},trial_ind.behav{1}{1},:),1))';
        Difference.behav{1,1}=actm1;Difference.behav{1,2}=actm2;
        %% behav: POSTGO-PRENOGO / 加alignwin后期-post
        % statistic test
        for neuroni=1:size(A_r,3)
            actm1=Difference.state{1,1}(neuroni,:);actm2=Difference.state{1,2}(neuroni,:);
            if  (~(sum(isnan(actm1))==length(actm1))  & ~(sum(isnan(actm2))==length(actm2)))
                [Difference.state{1,3}(neuroni,:),~]=ranksum(actm1,actm2);
            else
                Difference.state{1,3}(neuroni,:)=nan;
            end
            actm1=Difference.cue{1,1}(neuroni,:);actm2=Difference.cue{1,2}(neuroni,:);
            if  (~(sum(isnan(actm1))==length(actm1))  & ~(sum(isnan(actm2))==length(actm2)))
                [Difference.cue{1,3}(neuroni,:),~]=ranksum(actm1,actm2);
            else
                Difference.cue{1,3}(neuroni,:)=nan;
            end
            actm1= Difference.cue{2,1}(neuroni,:);actm2=Difference.cue{2,2}(neuroni,:);
            if  (~(sum(isnan(actm1))==length(actm1))  & ~(sum(isnan(actm2))==length(actm2)))
                [Difference.cue{2,3}(neuroni,:),~]=ranksum(actm1,actm2);
            else
                Difference.cue{2,3}(neuroni,:)=nan;
            end
            actm1=Difference.behav{1,1}(neuroni,:);actm2=Difference.behav{1,2}(neuroni,:);
            if (~(sum(isnan(actm1))==length(actm1))  & ~(sum(isnan(actm2))==length(actm2)))
                [Difference.behav{1,3}(neuroni,:),~]=ranksum(actm1,actm2);
            else
                Difference.behav{1,3}(neuroni,:)=nan;
            end
        end
        alpha = 0.05;
        Difference.state{1,4}=Difference.state{1,3}<alpha;
        Difference.cue{1,4}=Difference.cue{1,3}<alpha;Difference.cue{2,4}=Difference.cue{2,3}<alpha;
        Difference.behav{1,4}=Difference.behav{1,3}<alpha;
        % Bonferroni校正
        Difference.state{1,5}=Difference.state{1,3}*size(A_r,3);Difference.state{1,6}=Difference.state{1,5}<alpha;
        Difference.cue{1,5}=Difference.cue{1,3}*size(A_r,3);Difference.cue{1,6}=Difference.cue{1,5}<alpha;
        Difference.cue{2,5}=Difference.cue{2,3}*size(A_r,3);Difference.cue{2,6}=Difference.cue{2,5}<alpha;
        Difference.behav{1,5}=Difference.behav{1,3}*size(A_r,3);Difference.behav{1,6}=Difference.behav{1,5}<alpha;
        save([p,'3component_Singlecell.mat'],"Difference","trial_ind","period_ind",'-v7.3')
    end
end
%% visualization
% average trace/ heatmap / mapback/ brain-region/ functional type across conditioning
clc;clear all;close all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
fishi=9;batchi=1;p=char(Path{batchi}(fishi));
res=[0.66,0.66,10];radium=floor(20);
load([p,'3component_Singlecell.mat']);
load(fullfile(p,'para.mat'));
load(fullfile(p,'env.mat'));
load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
load(fullfile(p,'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
supervoxel(:,1)=env.supervoxel(:,1)*res(1);supervoxel(:,2)=env.supervoxel(:,2)*res(2);supervoxel(:,3)=(env.supervoxel(:,3)-1)*res(3);

A=zscore(activities_preCS_dfdf_aftcorrect(1:frame.per_cycle*trial.total,:),[],'all');
A_r=reshape(A,frame.per_cycle,trial.total,[]);
regionxx=categorical(Label);regionxx = reordercats(regionxx,Label);
stimCS=zeros(1,frame.per_cycle*trial.total);
for ii=1:size(stimCS,2)/frame.per_cycle
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=1;%23
end
stimUS=ones(1,frame.per_cycle*trial.total);
for ii=trial.acq(2):trial.acq(3)
    stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=3;%16
end
%figure,plot(stimCS,'r');hold on;plot(stimUS,'k')

Brain_region_id={};Index_in_region={};Fraction_in_region_type={};Loc_in_region_cell={};Num_in_region={};
for typei=1:3
    switch typei
        case 1
            triali=trial_ind.state;diffi=Difference.state;periodi=period_ind.state;leg='state';
            seg={'Pre Cond.','Post Cond.'};
        case 2
            triali=trial_ind.cue;diffi=Difference.cue;periodi=period_ind.cue;leg='cue';
            seg={'Pre Cond.','Post Cond.'};
        case 3
            triali=trial_ind.behav;diffi=Difference.behav;periodi=period_ind.behav;leg='behav';
            seg={'Go trials','Nogo trials'};
    end
    yli=[-0.2 0.8];
    ind=find(diffi{1,4}==1);locate=supervoxel(ind,1:3);
    act1=A_r(:, triali{1}{1},ind);act2=A_r(:, triali{1}{2},ind);actm1=squeeze(nanmean(act1,2));actm2=squeeze(nanmean(act2,2));
    act_all=A_r(:,:,ind);
    act_all_r=reshape(act_all,frame.per_cycle*trial.total,[]);
    if length(ind)<6000
        act_all_r_dr = tsne( act_all_r','Algorithm','exact','Distance','euclidean','Perplexity',50);
    else
        act_all_r_dr = tsne( act_all_r','Algorithm','exact','Distance','euclidean','Perplexity',ceil(size(A,2)/100));
    end
    %plot
    for ii=1 %average trace
        color=[0.5 0.5 0.5];
        figure,
        patch1=patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[min(yli(:)) min(yli(:)) max(yli(:)) max(yli(:))],...
            color,'edgecolor',color,'facecolor',color,'edgealpha',0,'facealpha',0.3);hold on
        shadedErrorBar([1:size(actm1,1)],actm2',{@mean,@(x) std(x)*(1*1/sqrt(length(ind)))},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.5);hold on
        shadedErrorBar([1:size(actm1,1)],actm1',{@mean,@(x) std(x)*(1*1/sqrt(length(ind)))},'lineprops',{[0.5,0.5,0.5]},'transparent',1,'patchSaturation',0.5);hold on
        legend(seg);set(gca,'fontsize',14,'FontName','Arial');ylim(yli);title(leg);
    end
    for ii=1 %heatmap
        a=nanmean( actm2(periodi{1,1},:),1);[~,inds]=sort(a,'descend');
        figure('position',[842,114,560,816]),
        subplot(1,2,1),imagesc(actm1(:,inds)',yli);colormap('hot');hold on;title(seg{1});
        line([frame.cs_start,frame.cs_end;frame.cs_start,frame.cs_end],[1,1;length(ind),length(ind)],'color','w','linewidth',2,'linestyle','--');
        ylabel(['Cell Num.' num2str(length(ind))]); set(gca,'fontsize',14,'FontName','Arial');
        subplot(1,2,2),imagesc(actm2(:,inds)',yli);colormap('hot');hold on;
        line([frame.cs_start,frame.cs_end;frame.cs_start,frame.cs_end],[1,1;length(ind),length(ind)],'color','w','linewidth',2,'linestyle','--');
        title(seg{2});set(gca,'fontsize',14,'FontName','Arial');
    end
    for ii=1 %mapback
        figure,
        scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),5,[0.5,0.5,0.5],'filled');hold on;
        scatter3(locate(:,1),locate(:,2),locate(:,3),10,'filled','r');title(leg);
        axis equal;xlim([0 2048*0.66]);ylim([0 2048*0.66]);set(gca,'fontsize',14,'FontName','Arial');view([250,50])
    end
    for ii=1 %density mapback
        figure,
        [In,On]=count_spatial_location_density(locate,supervoxel(:,:), radium);view([105,-45])
    end
    for ii=1 %brain-region
        Brain_region_id{typei}=brain_region_id(ind,:);
        for regioni=1:length(Label)
            if length(index_in_region_in_clust{regioni,1})<=10
                ind_inregion=[];
            else
                ind_inregion=intersect(ind,index_in_region_in_clust{regioni,1});
            end
            Index_in_region{typei}{regioni}=ind_inregion;
            Fraction_in_region_type{typei}(regioni)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
            Loc_in_region_cell{typei}{regioni}=supervoxel(ind_inregion,:);
            Num_in_region{typei}(regioni,1)= length(ind_inregion)./length(ind);
            Num_in_region{typei}(regioni,2)= length(ind_inregion);
            Num_in_region{typei}(regioni,3)= length(ind);
        end
        figure('position',[680,78,379,900]),
        myStackedBarwithErrorbarh(regionxx,Fraction_in_region_type{typei},zeros(1,length(Label)),[],[],regionxx);  hold on
    end
    for ii=1 %functional type across conditioning
        number_cluster=10;%%%%%%%换聚类
        [gIX,~,~ ]= kmeans(act_all_r',number_cluster,'Distance','sqeuclidean','Replicates',50,'MaxIter',10);
        clrmap = hsv(max(gIX)+3);clrmap = clrmap(1:end-3,:);%GetColormap('cool',max(gIX));
        cIX=1:size(act_all_r,2);
        figure, gscatter(act_all_r_dr(:,1),act_all_r_dr(:,2),gIX,clrmap);
        xlabel('Tsne-1');ylabel('Tsne-2');set(gca,'fontsize',16);legend off
        [h]=pushbutton_popupplot_Callback(act_all_r',cIX,gIX,clrmap,[],fs.ca,stimCS,stimUS,1,'all',0,0,0,0);
        [h,ratio]=plot_test_2(act_all_r,cIX,gIX,frame,clrmap,4,[-2 10],true,true);
    end
    for ii=1 %functional type mapback
        figure('Position',[282,93,1103,772]),
        scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),5,[0.5,0.5,0.5],'filled');hold on;
        for zz=unique(gIX)'
            id=find(gIX==zz);
            scatter3(locate(id,1),locate(id,2),locate(id,3),10,repmat(clrmap(zz,:),length(id),1),'filled');title(leg);
        end
        axis equal;xlim([0 2048*0.66]);ylim([0 2048*0.66]);set(gca,'fontsize',14,'FontName','Arial');view([250,50])
    end
end
%%%overlap
for typei=1
    switch typei
        case 1
            triali=trial_ind.state;diffi=Difference.state;periodi=period_ind.state;leg='state';
            seg={'Pre Cond.','Post Cond.'};
        case 2
            triali=trial_ind.cue;diffi=Difference.cue;periodi=period_ind.cue;leg='cue';
            seg={'Pre Cond.','Post Cond.'};
        case 3
            triali=trial_ind.behav;diffi=Difference.behav;periodi=period_ind.behav;leg='behav';
            seg={'Go trials','Nogo trials'};
    end
    yli=[-0.2 0.8];
    ind1=find(diffi{1,4}==1);locate=supervoxel(ind,1:3);
end
for typei=2
    switch typei
        case 1
            triali=trial_ind.state;diffi=Difference.state;periodi=period_ind.state;leg='state';
            seg={'Pre Cond.','Post Cond.'};
        case 2
            triali=trial_ind.cue;diffi=Difference.cue;periodi=period_ind.cue;leg='cue';
            seg={'Pre Cond.','Post Cond.'};
        case 3
            triali=trial_ind.behav;diffi=Difference.behav;periodi=period_ind.behav;leg='behav';
            seg={'Go trials','Nogo trials'};
    end
    yli=[-0.2 0.8];
    ind2=find(diffi{1,4}==1);
end
ind3=intersect(ind1,ind2);locate3=supervoxel(ind3,1:3);
figure,
scatter3(supervoxel(:,1),supervoxel(:,2),supervoxel(:,3),5,[0.5,0.5,0.5],'filled');hold on;
scatter3(locate3(:,1),locate3(:,2),locate3(:,3),10,'filled','r');title(leg);

%% 分7 type加brain——id



