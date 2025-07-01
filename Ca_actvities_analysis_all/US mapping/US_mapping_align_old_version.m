function US_mapping_align(Path)
set(0,'defaultfigurecolor','w')
addpath(genpath('F:\DUlab\FC_analyse\FishExplorer'));
if isempty(Path)
    p=uigetdir('H:\1.Test US\2.Tail free——Data from 117\')
else
    p=Path
end
if ~(exist([p,'CR_ind_summary_align.mat'],'file')==2)
    if ~[exist(fullfile(p,'activities_dfdf_align.mat'),'file') && exist(fullfile(p,'para.mat'),'file') && exist(fullfile(p,'env.mat'),'file')]
        warning([p '_file not exist']);
    else
        load(fullfile(p,'activities_dfdf_align.mat'),'activities_dfdf_align');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'env.mat'));
        %trial.total=trial.total-4;trial.test(3)=trial.test(3)-4;
         A_r=activities_dfdf_align;
         A=reshape(activities_dfdf_align,size(A_r,1)*size(A_r,2),[]);

        res=[0.66,0.66,10];
        load('X:\calcium data 20230224\脑区分割\segmentation_file_0525_DSregion_mask.mat');
        temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');
        warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
        temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
        nn=[Path(end-14:end-7),Path(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        
        trial.total=size(A_r,2);
        trial.acq_block_num=size(A_r,2)-2;
        trial.acq_block_trial=1;
        trial.hab=[1,1,1];
        trial.acq=[size(A_r,2)-2,2,size(A_r,2)-1];
        trial.test=[1,size(A_r,2),size(A_r,2)];
        trial_hab=trial.hab(2):trial.hab(3);
        trial_test=trial.test(2):trial.test(3);
        supervoxel_raw=env.supervoxel;
        for ii=unique(supervoxel_raw(:,3))'
            id=find(supervoxel_raw(:,3)==ii);
            supervoxel_raw(id,3)=(ii-1)*10/0.66+1;
        end
        supervoxel_raw(:,1)=env.width-supervoxel_raw(:,1);
        stimUS=zeros(1,frame.per_cycle*trial.total);ind_US=[];kk=1;
        trial.acq_block_interval=0;frame.us_dur=1;
        for ss=1:trial.acq_block_num
            for tt=1:trial.acq_block_trial
                ind=(trial.hab(3)+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frame.per_cycle+frame.us_start;
                stimUS(ind)=2;%23
                ind_US(kk,1)=ind;kk=kk+1;
            end
        end
        stimCS=zeros(1,frame.per_cycle*trial.total);ind_CS=[];kk=1;
        for ss=1:trial.total
            ind=(ss-1)*frame.per_cycle+[frame.cs_start:frame.cs_end];
            stimCS(ind)=1;%23
            ind_CS(kk,1)=ind(1);ind_CS(kk,2)=ind(end);kk=kk+1;
        end
        %figure, plot(stimCS);hold on;plot(stimUS);xlim([0 2700])
        colorCS=[0.9 0.5 0.5];
        nCells_total = size(A,2);prct_const = 2 ;topN = round(prct_const/100 * nCells_total);
        %% 有无CR in hab & tst
        pre_cs_time=[6 1];
        base_win=frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca;
        us_win=frame.cs_start:frame.us_start-1;
        CR_ind_up=[];CR_ind_down=[];CS_response=cell(2,2);
        for ii=1:2
            switch ii
                case 1
                    ind=trial_hab;
                case 2
                    ind=trial_test;
            end
            a=[];
            if isempty(a) || length(a)==length(ind)
                ind_base=ind;
            else
                ind_base=setdiff(ind,a);
            end
            aa=squeeze(mean(A_r(:,ind_base,:),2));
            m=mean(aa(base_win,:),1);sd=std(aa(base_win,:),[],1);
            %         aa=squeeze(reshape(A_r(base_win,ind_base,:),length(base_win)*length(ind_base),[]));
            %         m=mean(aa,1,'omitnan');sd=std(aa,[],1,'omitnan');
            if ii==1
                baseline_up=m+ 2*sd; baseline_down=m+ 2*sd;baseline_UR=m+2*sd;
            end
            aa=A_r(:,ind,:);
            for iscutmov=1 %1:no_cut;2:cut
                b=nan(length(us_win),length(ind),size(A,2));
                for jj=1:length(ind)
                        b(:,jj,:)=aa(us_win,jj,:);
                end
                m_CS=squeeze(mean(mean(b,1,'omitnan'),2,'omitnan'))';
                CS_response{ii,iscutmov}(:,:,:)=b;
                CR_ind_up(:,ii,iscutmov)=m_CS > baseline_up;
                CR_ind_down(:,ii,iscutmov)=m_CS < - baseline_down;
            end
        end
        %% 有无CR in acq
        us_win=frame.cs_start:frame.us_start-1;
        CR_ind_up_acq=[];CR_ind_down_acq=[];CS_response_acq=cell(trial.acq_block_num,2);
        for ii=1:trial.acq_block_num
            ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-1)*trial.acq_block_trial;
            aa=A_r(:,ind,:);
            for iscutmov=1
                b=nan(length(us_win),length(ind),size(A,2));
                for jj=1:length(ind)
                        b(:,jj,:)=aa(us_win,jj,:);
                end
                m_CS=squeeze(mean(mean(b,1,'omitnan'),2,'omitnan'))';
                CS_response_acq{ii,iscutmov}=b;
                CR_ind_up_acq(:,ii,iscutmov)=m_CS > baseline_up;
                CR_ind_down_acq(:,ii,iscutmov)=m_CS < - baseline_down;
            end
        end
        %% 有无UR
        us_win=frame.us_start:frame.us_start+2/fs.ca;
        UR_ind=[];US_response=cell(trial.acq_block_num,1);
        for ii=1:trial.acq_block_num
            ind=(trial.acq(2):trial.acq(2)+trial.acq_block_trial-1)+(ii-1)*trial.acq_block_trial;
            aa=squeeze(mean(A_r(:,ind,:),2));
            %figure,plot(aa(:,3000));
            %m=mean(aa(base_win,:),1);sd=std(aa(base_win,:),[],1);
            m_US=mean(aa(us_win,:),1);sd_US=std(aa(us_win,:),1);
            UR_ind(:,ii)=m_US> baseline_UR;
            US_response{ii}=A_r(us_win,ind,:);
        end
        ind=find(sum(UR_ind(:,1:max(ceil(size(UR_ind,2)/2),1)),2)>=2);
        h1=figure;
        scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,'r','filled');hold on
        axis equal;legend('US-responsive');
        %saveas(h1,fullfile(p,'US-responsive'),'png');
        view(90,90);
        figure,plot([1:size(A,1)].*fs.ca,mean(A(:,ind),2),'linewidth',2,'color','k');hold on;
        y=[-0.5 1.5];
        patch1=patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
            [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
            [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
            [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']'.*fs.ca,...
            repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
            colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
        line([frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1));frame.us_start+trial.hab(1)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab(1)+trial.acq(1))].*fs.ca,...
            [min(y) max(y)],'color','r','linestyle','--','linewidth',2);hold on;
        xlim([1 size(A,1)].*fs.ca);ylim(y);set(gca,'fontsize',16);box off;
        
        %% hab期间 CS引起视觉反应
        trial_ind=trial_hab;
        bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
        y=reshape(bb,[],1);
        for ii=1
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=[stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            thr=mean(stimcorr)+2*std(stimcorr);
            figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
            CS_related_hab_ind=find(stimcorr>thr);
            CS_related_hab=stimcorr;
        end
        %% Acq期间 CS引起视觉反应
        CS_related_acq=[];CS_related_acq_ind={};
        for ii=1:trial.acq_block_num
            trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-1)*trial.acq_block_trial;
            bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
            y=reshape(bb,[],1);
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=[stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
            CS_related_acq_ind{ii}=find(stimcorr>thr);
            CS_related_acq(:,ii)=stimcorr';
        end
        %% test期间 CS引起视觉反应
        trial_ind=trial_test;
        bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
        y=reshape(bb,[],1);
        for ii=1
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=[stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            %thr=0.2068;mean(stimcorr)+1*std(stimcorr)
            figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
            CS_related_tst_ind=find(stimcorr>thr);
            CS_related_tst=stimcorr;
        end
        set(gca,'fontsize',16)
        %% hab期间 CS引起向下的视觉反应
        trial_ind=trial_hab;
        bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
        y=reshape(bb,[],1);
        for ii=1
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=-[stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            thr2=mean(stimcorr)+2*std(stimcorr);
            figure,hist(stimcorr);hold on;line([thr,thr],[0 2000],'color','r','linewidth',3);
            CS_related_hab_down_ind=find(stimcorr>thr);
            CS_related_hab_down=stimcorr;
        end
        %% Acq期间 CS引起向下的视觉反应
        CS_related_acq_down_ind={};CS_related_acq_down=[];
        for ii=1:trial.acq_block_num
            trial_ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-1)*trial.acq_block_trial;
            bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
            y=reshape(bb,[],1);
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=[-stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            figure,hist(stimcorr);hold on;line([thr2,thr2],[0 2000],'color','r','linewidth',3);
            CS_related_acq_down_ind{ii}=find(stimcorr>thr2);
            CS_related_acq_down(:,ii)=stimcorr';
        end
        %% test 期间 CS引起向下的视觉反应
        trial_ind=trial_test;
        bb=reshape(stimCS,frame.per_cycle,trial.total);bb=bb(:,trial_ind);
        nCells_total = size(A,2);prct_const = 2 ;topN = round(prct_const/100 * nCells_total);
        y=reshape(bb,[],1);
        for ii=1
            stim=zeros(5,length(trial_ind)*frame.per_cycle);stim(1,:)=y;
            [stregressors,~,~,~] = GetMotorRegressor(stim,2);
            stAvrcorr_all=[];
            stim_output=-[stregressors(1,1).im];
            M=A_r(:,trial_ind,:);M=reshape(M,frame.per_cycle*length(trial_ind),[]);
            [stimcorr,~] = MotorSourceCorrelation(M',stim_output,[]);
            figure,plot(stim(1,:));hold on;
            plot(stim_output,'r');legend('raw','regressor')
            figure,hist(stimcorr);hold on;line([thr2,thr2],[0 2000],'color','r','linewidth',3);
            CS_related_tst_down_ind=find(stimcorr>thr2);
            CS_related_tst_down=stimcorr;
        end
        set(gca,'fontsize',16)
        
        close all
        
        %% save
        corr_CS_up=[];corr_CS_down=[];
        corr_CS_up(:,1)=CS_related_hab;
        corr_CS_up(:,2:1+trial.acq_block_num)=CS_related_acq;
        corr_CS_up(:,trial.acq_block_num+2)=CS_related_tst;
        corr_CS_down(:,1)=CS_related_hab_down;
        corr_CS_down(:,2:1+trial.acq_block_num)=CS_related_acq_down;
        corr_CS_down(:,trial.acq_block_num+2)=CS_related_tst_down;
        
        %% part of plot from mapping_index_plot.m
        labels_typei={'CS-avtivation','CS-inhibition','US-activation','CS-up regulate','CS-down regulate','CS-stable regulate','US-up regulate','US-down regulate','US-stable regulate','CS-new emerged neuron'};
        NUM_DOWN_CR=[]; NUM_UP_CR=[]; NUM_UP_UR=[];NUM_NEW_EMERGED=[];
        num_of_ind_all_emergeCSUS=nan(length(labels_typei),trial.acq_block_num+2,trial.acq_block_num+2,2,3);
        if exist('p_CS_up_bef_aft')==1
            if sum(isnan(p_CS_up_bef_aft(:,1)))<size(A,2);
                ranksumtest=0;
            else
                p_CS_up_bef_aft=nan(size(A,2),2);p_CS_down_bef_aft=nan(size(A,2),2);
                p_CS_up_bef_cond=nan(size(A,2),trial.acq_block_num,2);p_CS_down_bef_cond=nan(size(A,2),trial.acq_block_num,2);
                p_US_up_bef_cond=nan(size(A,2),trial.acq_block_num);p_US_down_bef_cond=nan(size(A,2),trial.acq_block_num);
                ranksumtest=1;
            end
        else
            p_CS_up_bef_aft=nan(size(A,2),2);p_CS_down_bef_aft=nan(size(A,2),2);
            p_CS_up_bef_cond=nan(size(A,2),trial.acq_block_num,2);p_CS_down_bef_cond=nan(size(A,2),trial.acq_block_num,2);
            p_US_up_bef_cond=nan(size(A,2),trial.acq_block_num);p_US_down_bef_cond=nan(size(A,2),trial.acq_block_num);
            ranksumtest=1;
        end
        q_CS_up_bef_aft=p_CS_up_bef_aft;q_CS_up_bef_cond=p_CS_up_bef_cond;q_US_up_bef_cond=p_US_up_bef_cond;q_CS_down_bef_cond=p_CS_down_bef_cond;q_US_down_bef_cond=p_US_down_bef_cond;q_CS_down_bef_aft=p_CS_down_bef_aft;
        fdr_CS_up_bef_aft=p_CS_up_bef_aft;fdr_CS_up_bef_cond=p_CS_up_bef_cond;fdr_US_up_bef_cond=p_US_up_bef_cond;fdr_CS_down_bef_cond=p_CS_down_bef_cond;fdr_US_down_bef_cond=p_US_down_bef_cond;fdr_CS_down_bef_aft=p_CS_down_bef_aft;
        CS_up_ind=cell(trial.acq_block_num+2,2); CS_down_ind=cell(trial.acq_block_num+2,2); CS_stable_ind=cell(trial.acq_block_num+2,2);
        CS_up=nan(trial.acq_block_num+2,2);CS_down=nan(trial.acq_block_num+2,2);CS_stable=nan(trial.acq_block_num+2,2);
        US_up_ind=cell(trial.acq_block_num+2,2); US_down_ind=cell(trial.acq_block_num+2,2); US_stable_ind=cell(trial.acq_block_num+2,2);
        US_up=nan(trial.acq_block_num+2,2);US_down=nan(trial.acq_block_num+2,2);US_stable=nan(trial.acq_block_num+2,2);
        ind_all_CSUS_RESPONSIVE=cell(length(labels_typei),trial.acq_block_num+2,2);
        ind_all_emergeCSUS=cell(length(labels_typei),trial.acq_block_num+2,trial.acq_block_num+2,2,3);
        Fraction_in_region_type=cell(length(labels_typei),2,2);Loc_in_region=cell(length(labels_typei),trial.acq_block_num+2,2,2);
        Loc_in_region_cell=cell(length(labels_typei),trial.acq_block_num+2,2,2);Num_in_region=cell(length(labels_typei),2,2,3);
        Brain_region_id=cell(length(labels_typei),trial.acq_block_num+2,2,2);Label_region=cell(length(labels_typei),2,2);
        Fraction_in_region_type_emerged=cell(length(labels_typei),trial.acq_block_num+2,2,3,2);Loc_in_region_emerged=cell(length(labels_typei),trial.acq_block_num+2,trial.acq_block_num+2,2,3,2);
        Loc_in_region_cell_emerged=cell(length(labels_typei),trial.acq_block_num+2,trial.acq_block_num+2,2,3,2);Num_in_region_emerged=cell(length(labels_typei),trial.acq_block_num+2,2,3,2,3);
        Brain_region_id_emerged=cell(length(labels_typei),trial.acq_block_num+2,trial.acq_block_num+2,2,3,2);Label_region_emerged=cell(length(labels_typei),trial.acq_block_num+2,2,3,2);
        ind_type_sort={};ind_type_sort_regionlabel={};
        sessionx = {'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Cond.7','Cond.8','Post Cond'};
        for iscutmov=1
            if isfolder(fullfile(p,['iscutmov',num2str(iscutmov)]))
                try
                    rmdir(fullfile(p,['iscutmov',num2str(iscutmov)]),'s');
                catch
                    warning(fullfile(p,['iscutmov',num2str(iscutmov)]))
                end
            end
            checkpath(fullfile(p,['iscutmov',num2str(iscutmov)]));
            v=[90,90];
            radium=floor(15/0.66);
            % CR in hab
            ind=find(CR_ind_up(:,1, iscutmov)==1);
            %ind=intersect(ind,CS_related_hab_ind);
            ind2=find(CR_ind_down(:,1,iscutmov)==1);
            %ind2=intersect(ind2,CS_related_hab_down_ind);
            ind3=union(ind,ind2);
            ind3=setdiff(1:size(A,2),ind3);
            h2=figure;t2=tiledlayout(h2,2,2);t2.TileSpacing = 'compact';t2.Padding = 'compact';
            ax21=nexttile(t2);a=[0 0.03];
            patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
                [min(a) min(a) max(a) max(a)],...
                colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
            trial_ind=[trial.hab(2):trial.hab(3)];
            a=mean(mean(A_r(:,trial_ind,ind3),3),2);plot(a,'color',[0.5 0.5 0.5],'linewidth',2)
            b=mean(mean(A_r(:,trial_ind,ind),3),2);plot(b,'color','r','linewidth',2)
            c=mean(mean(A_r(:,trial_ind,ind2),3),2);plot(c,'color','b','linewidth',2)
            ax22=nexttile(t2);b=bar([length(ind3),length(ind2),length(ind)]./size(A,2));
            xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;labels1 = num2str(b(1).YData','%0.4f');
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
            
            h=figure('position',[313,74,1099,719]);t=tiledlayout(h,3,2);t.TileSpacing = 'compact';t.Padding = 'compact';
            ax1=nexttile(t,[2 1]);
            scatter3(supervoxel_raw(ind3,1),supervoxel_raw(ind3,2),supervoxel_raw(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
            scatter3(supervoxel_raw(ind2,1),supervoxel_raw(ind2,2),supervoxel_raw(ind2,3),10,'b','filled');hold on
            scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,'r','filled');hold on
            grid off;set(gca,'visible','off');
            axis equal;legend('CS unresponsive in hab','CS inhibition in hab','CS activation in hab');axis equal;grid off;
            view(v(1),v(2))
            ax1=nexttile(t,5,[1 1]);
            scatter3(supervoxel_raw(ind3,1),supervoxel_raw(ind3,2),supervoxel_raw(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
            scatter3(supervoxel_raw(ind2,1),supervoxel_raw(ind2,2),supervoxel_raw(ind2,3),10,'b','filled');hold on
            scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,'r','filled');hold on
            grid off;set(gca,'visible','off');
            axis equal;axis equal;grid off;
            view(0,0)
            
            % CR in test
            ind=find(CR_ind_up(:,2,iscutmov)==1);
            %ind=intersect(ind,CS_related_tst_ind);%ind=setdiff(ind,CSmotor_related_tst_ind );
            ind2=find(CR_ind_down(:,2,iscutmov)==1);
            %ind2=intersect(ind2, CS_related_tst_down_ind);
            ind3=union(ind,ind2);
            ind3=setdiff(1:size(A,2),ind3);
            
            ax23=nexttile(t2);
            patch1=patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start],...
                [min(a) min(a) max(a) max(a)],...
                colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
            trial_ind=[trial.test(2):trial.test(3)];
            a=mean(mean(A_r(:,trial_ind,ind3),3),2);plot(a,'color',[0.5 0.5 0.5],'linewidth',2)
            b=mean(mean(A_r(:,trial_ind,ind),3),2);plot(b,'color','r','linewidth',2)
            c=mean(mean(A_r(:,trial_ind,ind2),3),2);plot(c,'color','b','linewidth',2)
            ax24=nexttile(t2);b=bar([length(ind3),length(ind2),length(ind)]./size(A,2));
            xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;labels1 = num2str(b(1).YData','%0.4f');
            text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
            ax2=nexttile(t,[2 1]);
            scatter3(supervoxel_raw(ind3,1),supervoxel_raw(ind3,2),supervoxel_raw(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
            scatter3(supervoxel_raw(ind2,1),supervoxel_raw(ind2,2),supervoxel_raw(ind2,3),10,'b','filled');hold on
            scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,'r','filled');hold on
            axis equal;legend('CS unresponsive in tst','CS inhibition in tst','CS activation in tst');
            title(p(end-13:end));
            view(v(1),v(2))
            ax2=nexttile(t,6,[1 1]);
            scatter3(supervoxel_raw(ind3,1),supervoxel_raw(ind3,2),supervoxel_raw(ind3,3),2,[0.5 0.5 0.5],'filled');hold on
            scatter3(supervoxel_raw(ind2,1),supervoxel_raw(ind2,2),supervoxel_raw(ind2,3),10,'b','filled');hold on
            scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,'r','filled');hold on
            axis equal;
            title(p(end-13:end));
            view(0,0)
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Location of CS responsive neuron']),'jpeg');
            saveas(h2,fullfile(p,['iscutmov',num2str(iscutmov)],['Responsive trace of CS']),'jpeg');
            %up CR in different session
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num+2
                switch ii
                    case 1
                        ind=find(CR_ind_up(:,1,iscutmov)==1);%ind=intersect(ind,CS_related_hab_ind);
                        leg='Pre';
                    case trial.acq_block_num+2
                        ind=find(CR_ind_up(:,2,iscutmov)==1);%ind=intersect(ind,CS_related_tst_ind);
                        leg='Post';
                    otherwise
                        ind=find(CR_ind_up_acq(:,ii-1,iscutmov)==1);%ind=intersect(ind,CS_related_acq_ind{ii-1});
                        leg=['Acq.' num2str(ii)];
                end
                NUM_UP_CR(ii,iscutmov)=length(ind);
                nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
                set(gca,'fontsize',16)
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of CS activation neuron']),'jpeg');
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num+2
                switch ii
                    case 1
                        ind=find(CR_ind_up(:,1,iscutmov)==1);cr=squeeze(mean(mean(CS_response{1,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_hab_ind);
                        leg='Pre';
                    case trial.acq_block_num+2
                        ind=find(CR_ind_up(:,2,iscutmov)==1);cr=squeeze(mean(mean(CS_response{2,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_tst_ind);
                        leg='Post';
                    otherwise
                        ind=find(CR_ind_up_acq(:,ii-1,iscutmov)==1);cr=squeeze(mean(mean(CS_response_acq{ii-1,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_acq_ind{ii-1});
                        leg=['Acq.' num2str(ii)];
                end
                nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(leg);
                set(gca,'fontsize',16)
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['CR AMP of CS activation neuron']),'jpeg');
            %down CR in different session
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num+2
                switch ii
                    case 1
                        ind=find(CR_ind_down(:,1,iscutmov)==1);%ind=intersect(ind,CS_related_hab_ind);
                        leg='Pre';
                    case trial.acq_block_num+2
                        ind=find(CR_ind_down(:,2,iscutmov)==1);%ind=intersect(ind,CS_related_tst_ind);
                        leg='Post';
                    otherwise
                        ind=find(CR_ind_down_acq(:,ii-1,iscutmov)==1);%ind=intersect(ind,CS_related_acq_ind{ii-1});
                        leg=['Acq.' num2str(ii)];
                end
                NUM_DOWN_CR(ii,iscutmov)=length(ind);
                nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
                set(gca,'fontsize',16)
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of CS inhibition neuron']),'jpeg');
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num+2
                switch ii
                    case 1
                        ind=find(CR_ind_down(:,1,iscutmov)==1);cr=squeeze(mean(mean(CS_response{1,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_hab_ind);
                        leg='Pre';
                    case trial.acq_block_num+2
                        ind=find(CR_ind_down(:,2,iscutmov)==1);cr=squeeze(mean(mean(CS_response{2,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_tst_ind);
                        leg='Post';
                    otherwise
                        ind=find(CR_ind_down_acq(:,ii-1,iscutmov)==1);cr=squeeze(mean(mean(CS_response_acq{ii-1,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));%ind=intersect(ind,CS_related_acq_ind{ii-1});
                        leg=['Acq.' num2str(ii)];
                end
                nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(leg);
                set(gca,'fontsize',16)
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['CR AMP of CS inhibition neuron']),'jpeg');
            %UR in different session
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num
                ind=find(UR_ind(:,ii)==1);    NUM_UP_UR(ii)=length(ind);
                nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));
                set(gca,'fontsize',16);title(['Acq.' num2str(ii)]);
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of US responsive neuron']),'jpeg');
            h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.acq_block_num+2)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
            for ii=1:trial.acq_block_num
                ind=find(UR_ind(:,ii)==1);
                nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),squeeze(mean(mean(US_response{ii,1}(:,:,ind),1,'omitnan'),2,'omitnan')));view(v(1),v(2));
                set(gca,'fontsize',16);title(['Acq.' num2str(ii)]);
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['UR AMP of US responsive neuron']),'jpeg');
            %bar of CR AND UR
            x=[0 NUM_UP_UR 0];
            h=figure;b=bar([NUM_UP_CR(:,iscutmov)';x;NUM_DOWN_CR(:,iscutmov)']'*100./size(A,2));
            xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;labels = num2str(b(1).YData','%0.2f');
            xtips3 = b(3).XEndPoints;ytips3 = b(3).YEndPoints;labels3 = num2str(b(3).YData','%0.2f');
            xtips2 = b(2).XEndPoints;ytips2 = b(2).YEndPoints;labels2 = num2str(b(2).YData','%0.2f');
            text(xtips1,ytips1,labels,'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom')
            ylabel('Neurons(%)');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);legend({'CS activation';'US';'CS inhibition'});box off;
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Fraction of CS and US responsive neuron']),'jpeg');
            
            %% newly emerged  distribution
            ind=find(CR_ind_up(:,1,iscutmov)==1);
            %ind=intersect(ind,CS_related_hab_ind);
            ind2=find(CR_ind_down(:,1,iscutmov)==1);
            %ind2=intersect(ind2,CS_related_hab_down_ind);
            ind3=union(ind,ind2);
            ind3=setdiff(1:size(CR_ind_up(:,1,iscutmov),1),ind3);
            ind2=find(CR_ind_up(:,2,iscutmov)==1);
            %ind2=intersect(ind2,CS_related_tst_ind);
            CS_new_emerged=intersect(ind2,ind3);
            
            h=figure('position',[35,531,1883,447]);
            ind=find(CR_ind_up(:,1,iscutmov)==1);
            %ind=intersect(ind,CS_related_hab_ind);
            subplot(1,3,1),count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);title('pre');view(v(1),v(2))
            ind=find(CR_ind_up(:,2,iscutmov)==1);
            %ind=intersect(ind,CS_related_tst_ind);
            subplot(1,3,2),count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);title('post');view(v(1),v(2))
            subplot(1,3,3),count_spatial_location_density(supervoxel_raw(CS_new_emerged,:),supervoxel_raw(:,:), radium);title('newly emerged');view(v(1),v(2))
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of newly emerged CS responsive neuron']),'jpeg');
            
            h=figure('position',[35,531,1883,447]);
            ind=find(CR_ind_up(:,1,iscutmov)==1);
            %ind=intersect(ind,CS_related_hab_ind);
            subplot(1,3,1),mapback_colorcode_responses(supervoxel_raw(ind,:),squeeze(mean(mean(CS_response{1,1}(:,:,ind),1,'omitnan'),2,'omitnan')));title('pre');view(v(1),v(2))
            ind=find(CR_ind_up(:,2,iscutmov)==1);
            %ind=intersect(ind,CS_related_tst_ind);
            subplot(1,3,2),mapback_colorcode_responses(supervoxel_raw(ind,:),squeeze(mean(mean(CS_response{2,1}(:,:,ind),1,'omitnan'),2,'omitnan')));title('post');view(v(1),v(2))
            subplot(1,3,3),mapback_colorcode_responses(supervoxel_raw(CS_new_emerged,:),squeeze(mean(mean(CS_response{2,1}(:,:,CS_new_emerged),1,'omitnan'),2,'omitnan')));title('newly emerged');view(v(1),v(2))
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['CR AMP of newly emerged CS responsive neuron']),'jpeg');
            close(h);
            %figure,scatplot(supervoxel(CS_new_emerged,1),supervoxel(CS_new_emerged,2),'circles', radium,100,[],[],[]);
            %heatmap of actvity
            for ii=1:2
                switch ii
                    case 1
                        trial_ind=[trial.hab(2):trial.hab(3)];leg='Pre';
                    case 2
                        trial_ind=[trial.test(2):trial.test(3)];leg='Post';
                end
                a=reshape(A_r(:,trial_ind,CS_new_emerged),[],length(CS_new_emerged));
                a=smoothdata(a,1,'movmedian',3);
                [idx_functional,~,~,~] = kmeans(a',4,'Distance','correlation','Replicates',10);
                [idxAs,idxAa_ind]=sort(idx_functional);
                h=figure;t=tiledlayout(5,1);%t.TileSpacing = 'compact';%t.Padding = 'compact';
                ax1=nexttile([1,1]);plot(mean(a(:,idxAa_ind),2,'omitnan'),'k','linewidth',2);xlim([1 frame.per_cycle*length(trial_ind)]);
                ylim([-0.008 2]);set(gca,'fontsize',16);hold on;
                y=[0.01 2];
                ax2=nexttile([4,1]);imagesc(a(:,idxAa_ind)',y);colormap('hot');hold on;
                line(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*length(trial_ind)],2,1),repmat([0 length(CS_new_emerged)],length(trial_ind),1)','color','w','linestyle','--','linewidth',1.5);
                linkaxes([ax1,ax2],'x');title(t,leg);xticklabels(ax1,{});
                set(gca,'fontsize',16);colorbar('Ticks',y,'Ticklabels',y*100);
                saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Trace of newly emerged CS responsive neuron-',leg]),'jpeg');
                close(h);
            end
            
            %% UP/DOWN CS/US-responsive neuron
            %down/up-regulate
            if ranksumtest==1;
                x=CS_response{1,iscutmov};y=CS_response{2,iscutmov};z=CS_response_acq;
                for ii=1:size(A,2)
                    q_CS_up_bef_aft(ii,iscutmov)=mean(squeeze(x(:,:,ii)),1,'omitnan') < mean(squeeze(y(:,:,ii)),1,'omitnan');
                    q_CS_down_bef_aft(ii,iscutmov)=mean(squeeze(x(:,:,ii)),1,'omitnan')> mean(squeeze(y(:,:,ii)),1,'omitnan');
                    for jj=1:trial.acq_block_num
                        q_CS_up_bef_cond(ii,jj,iscutmov)=mean(squeeze(x(:,:,ii)),1,'omitnan')<mean(squeeze(z{jj,iscutmov}(:,:,ii)),1,'omitnan');
                        q_CS_down_bef_cond(ii,jj,iscutmov)=mean(squeeze(x(:,:,ii)),1,'omitnan')>mean(squeeze(z{jj,iscutmov}(:,:,ii)),1,'omitnan');
                    end
                end
                %[FDR_base,p_CS_responsive_bef_aft_adj] = mafdr( p_CS_responsive_bef_aft);
                x=US_response{1};y=US_response;
                for ii=1:size(A,2)
                    for jj=1:trial.acq_block_num
                        q_US_up_bef_cond(ii,jj)=mean(squeeze(x(:,:,ii)),1,'omitnan')<mean(squeeze(y{jj}(:,:,ii)),1,'omitnan');
                        q_US_down_bef_cond(ii,jj)=mean(squeeze(x(:,:,ii)),1,'omitnan')>mean(squeeze(y{jj}(:,:,ii)),1,'omitnan');
                    end
                end
            end
%             for ii=1
%                 [fdr_CS_up_bef_aft,q_CS_up_bef_aft] = mafdr(squeeze(p_CS_up_bef_aft(:,iscutmov)),'Method','polynomial');
%                 [fdr_CS_down_bef_aft,q_CS_down_bef_aft,~,~] = mafdr(squeeze(p_CS_down_bef_aft(:,iscutmov)),'BHFDR',true);
%             end
%             for jj=1:trial.acq_block_num
%                 [fdr_CS_up_bef_cond,q_CS_up_bef_cond,~,~] = mafdr(p_CS_up_bef_cond,'BHFDR',true);
%                 [fdr_CS_down_bef_cond,q_CS_down_bef_cond,~,~] = mafdr(p_CS_down_bef_cond,'BHFDR',true);
%                 [fdr_US_up_bef_cond,q_US_up_bef_cond,~,~] = mafdr(p_US_up_bef_cond,'BHFDR',true);
%                 [fdr_US_down_bef_cond,q_US_down_bef_cond,~,~] = mafdr(p_US_down_bef_cond,'BHFDR',true);
%             end
      
            %plot pie
            explode = [1 1 0 ];labels_typei = {'Up','Down','Stable'};
            h=figure('position',[14,187,1833,716]);
            CS_all=length(find(sum(CR_ind_up(:,:,iscutmov),2)>=1 | sum(CR_ind_up_acq(:,:,iscutmov),2)>=1));
            US_all=length(find(sum(UR_ind,2)>=1));
            for ii=1:trial.acq_block_num+2
                switch ii
                    case trial.acq_block_num+2
                        CS_up_ind{ii,iscutmov}=find(squeeze(q_CS_up_bef_aft(:,iscutmov))<=0.05 & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        CS_down_ind{ii,iscutmov}=find(squeeze(q_CS_down_bef_aft(:,iscutmov))<=0.05 & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        CS_stable_ind{ii,iscutmov}=find((squeeze(q_CS_up_bef_aft(:,iscutmov))>0.05 & squeeze(q_CS_down_bef_aft(:,iscutmov))>0.05) & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        US_up_ind{ii,iscutmov}=[];
                        US_down_ind{ii,iscutmov}=[];
                        US_stable_ind{ii,iscutmov}=[];
                    case 1
                        CS_up_ind{ii,iscutmov}=find(squeeze(q_CS_up_bef_aft(:,iscutmov))<=0.05 & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        CS_down_ind{ii,iscutmov}=find(squeeze(q_CS_down_bef_aft(:,iscutmov))<=0.05 & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        CS_stable_ind{ii,iscutmov}=find((squeeze(q_CS_up_bef_aft(:,iscutmov))>0.05 & squeeze(q_CS_down_bef_aft(:,iscutmov))>0.05) & sum(squeeze(CR_ind_up(:,:,iscutmov)),2)>=1);
                        US_up_ind{ii,iscutmov}=[];
                        US_down_ind{ii,iscutmov}=[];
                        US_stable_ind{ii,iscutmov}=[];
                    otherwise
                        CS_up_ind{ii,iscutmov}=find(squeeze(q_CS_up_bef_cond(:,ii-1,iscutmov))<=0.05 & squeeze(CR_ind_up_acq(:,ii-1,iscutmov)+ CR_ind_up(:,1,iscutmov))>=1);
                        CS_down_ind{ii,iscutmov}=find(squeeze(q_CS_down_bef_cond(:,ii-1,iscutmov))<=0.05 &  squeeze(CR_ind_up_acq(:,ii-1,iscutmov)+ CR_ind_up(:,1,iscutmov))>=1);
                        CS_stable_ind{ii,iscutmov}=find((squeeze(q_CS_up_bef_cond(:,ii-1,iscutmov))>0.05 & squeeze(q_CS_down_bef_cond(:,ii-1,iscutmov))>0.05) & squeeze(CR_ind_up_acq(:,ii-1,iscutmov)+ CR_ind_up(:,1,iscutmov))>=1);
                        US_up_ind{ii,iscutmov}=find(squeeze(q_US_up_bef_cond(:,ii-1))<=0.05 & squeeze(UR_ind(:,ii-1)+UR_ind(:,1))>=1);
                        US_down_ind{ii,iscutmov}=find(squeeze(q_US_down_bef_cond(:,ii-1))<=0.05 & squeeze(UR_ind(:,ii-1)+UR_ind(:,1))>=1);
                        US_stable_ind{ii,iscutmov}=find((squeeze(q_US_up_bef_cond(:,ii-1))>0.05 & squeeze(q_US_down_bef_cond(:,ii-1))>0.05) & squeeze(UR_ind(:,ii-1)+UR_ind(:,1))>=1);
                end
                CS_up(ii,iscutmov)=length(CS_up_ind{ii,iscutmov})/CS_all;
                CS_down(ii,iscutmov)=length(CS_down_ind{ii,iscutmov})/CS_all;
                CS_stable(ii,iscutmov)=length(CS_stable_ind{ii,iscutmov})/CS_all;
                US_up(ii,iscutmov)=length(US_up_ind{ii,iscutmov})/US_all;
                US_down(ii,iscutmov)=length(US_down_ind{ii,iscutmov})/US_all;
                US_stable(ii,iscutmov)=length(US_stable_ind{ii,iscutmov})/US_all;
                subplot(2,trial.acq_block_num+2,ii),
                pie([CS_up(ii,iscutmov) CS_down(ii,iscutmov) CS_stable(ii,iscutmov)],explode);title(['CS Session',num2str(ii)]);
                set(gca,'fontsize',14,'FontWeight','bold','linewidth',2)
                %subplot(2,trial.acq_block_num+1,ii+trial.acq_block_num),pie([US_up(ii,iscutmov) US_down(ii,iscutmov) US_stable(ii,iscutmov)],explode);legend(labels,'position',[0.008,0.026,0.12,0.28]);title(['US Session',num2str(ii)]);set(gca,'fontsize',14,'FontWeight','bold','linewidth',2)
                if ii==trial.acq_block_num+2
                    subplot(2,trial.acq_block_num+2,ii+trial.acq_block_num+2),
                    title('Responsive');b=bar(categorical({'CS','US','All'}),[CS_all US_all size(A,2)]);
                    xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;label = string(b(1).YData);
                    text(xtips1,ytips1,label,'HorizontalAlignment','center',...
                        'VerticalAlignment','bottom');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off
                else
                    subplot(2,trial.acq_block_num+2,ii+trial.acq_block_num+2),pie([US_up(ii,iscutmov) US_down(ii,iscutmov) US_stable(ii,iscutmov)],explode);legend(labels_typei,'position',[0.008,0.026,0.12,0.28]);
                    title(['US Session',num2str(ii)]);set(gca,'fontsize',14,'FontWeight','bold','linewidth',2)
                end
            end
            saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Fraction of Up-down regulate CS-US responsive neuron']),'jpeg');
            close(h);
            %% 脑区统计
            %cs-responsive
            %us-responsive
            %cs/us up/down/stable regulate
            labels_typei={'CS-avtivation','CS-inhibition','US-activation','CS-up regulate','CS-down regulate','CS-stable regulate','US-up regulate','US-down regulate','US-stable regulate','CS-new emerged neuron'};
            for ii=1:trial.acq_block_num+2
                switch ii
                    case 1
                        ind_all_CSUS_RESPONSIVE{1,ii,iscutmov}=find(CR_ind_up(:,1,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{2,ii,iscutmov}=find(CR_ind_down(:,1,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{3,ii,iscutmov}=[];
                    case trial.acq_block_num+2
                        ind_all_CSUS_RESPONSIVE{1,ii,iscutmov}=find(CR_ind_up(:,2,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{2,ii,iscutmov}=find(CR_ind_down(:,2,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{3,ii,iscutmov}=[];
                    otherwise
                        ind_all_CSUS_RESPONSIVE{1,ii,iscutmov}=find(CR_ind_up_acq(:,ii-1,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{2,ii,iscutmov}=find(CR_ind_down_acq(:,ii-1,iscutmov)==1);
                        ind_all_CSUS_RESPONSIVE{3,ii,iscutmov}=find(UR_ind(:,ii-1)==1);
                end
                ind_all_CSUS_RESPONSIVE{4,ii,iscutmov}=CS_up_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{5,ii,iscutmov}=CS_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{6,ii,iscutmov}=CS_stable_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{7,ii,iscutmov}=US_up_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{8,ii,iscutmov}=US_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{9,ii,iscutmov}=US_stable_ind{ii,iscutmov};
                ind=find(CR_ind_up(:,1,iscutmov)==1);
                ind2=find(CR_ind_down(:,1,iscutmov)==1);
                ind3=union(ind,ind2);
                ind3=setdiff(1:size(CR_ind_up(:,1,iscutmov),1),ind3);
                ind_all_CSUS_RESPONSIVE{10,ii,iscutmov}=intersect(ind_all_CSUS_RESPONSIVE{1,ii,iscutmov},ind3);
                NUM_NEW_EMERGED(ii,iscutmov)=length(ind_all_CSUS_RESPONSIVE{10,ii,iscutmov});
            end
            for ii=1:1:trial.acq_block_num+2
                for jj=1:length(labels_typei)
                    for zz=1:trial.acq_block_num+2
                        %                     if jj>=3 && jj<=10
                        %                         a=ind_all_CSUS_RESPONSIVE{jj,2,iscutmov};
                        %                     else
                        %                         a=ind_all_CSUS_RESPONSIVE{jj,1,iscutmov};
                        %                     end
                        a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        b=ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov};
                        c=ind_all_CSUS_RESPONSIVE{jj,end,iscutmov};
                        index_in_session_precond=setdiff(b,a);
                        index_in_session_postcond=setdiff(b,c);
                        index_in_session_prepostcond=intersect(index_in_session_precond,index_in_session_postcond);
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}=index_in_session_precond;
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,2}=index_in_session_postcond;
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,3}=index_in_session_prepostcond;
                        num_of_ind_all_emergeCSUS(jj,ii,zz,iscutmov,1)=length(index_in_session_precond);
                        num_of_ind_all_emergeCSUS(jj,ii,zz,iscutmov,2)=length(index_in_session_postcond);
                        num_of_ind_all_emergeCSUS(jj,ii,zz,iscutmov,3)=length(index_in_session_prepostcond);
                    end
                end
            end
            for jj=1:length(labels_typei)
                for ii=1:trial.acq_block_num+2
                    ind=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                    if ~isempty(ind)
                        gIX=ones(1,length(ind))';supervoxel_in_fish=supervolxeli(ind,1:3);clrmap = GetColormap('hsv_new',max(gIX));
                        for fraction_type=1:2
                            switch fraction_type
                                case 1
                                    temp_vol=supervolxeli(:,1:3);
                                case 2
                                    temp_vol=temp_supervoxel;
                            end
                            [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervoxel_in_fish,reg_mask,reg_name,reg_loc,temp_vol,temp_env,clrmap,nn);
                            Fraction_in_region_type{jj,iscutmov,fraction_type}(:,ii)=fraction_in_region_in_clust(:,1);
                            Loc_in_region{jj,ii,iscutmov,fraction_type}=loc_in_region_in_clust;
                            Loc_in_region_cell{jj,ii,iscutmov,fraction_type}=loc_in_region_cell;
                            Num_in_region{jj,iscutmov,fraction_type,1}(:,ii)= num_in_region_in_clust{1}(1,:);
                            Num_in_region{jj,iscutmov,fraction_type,2}(:,ii)= num_in_region_in_clust{2}(1,:);
                            Num_in_region{jj,iscutmov,fraction_type,3}(:,ii)= num_in_region_in_clust{3}(1,:);
                            Brain_region_id{jj,ii,iscutmov,fraction_type}=brain_region_id(:,1);
                            Label_region{jj,iscutmov,fraction_type}(:,ii)=Label;
                        end
                    end
                end
            end
            for zz=1:trial.acq_block_num+2
                for jj=1:length(labels_typei)
                    for ii=1:trial.acq_block_num+2
                        for typei=1:size(ind_all_emergeCSUS,5)
                            ind=ind_all_emergeCSUS{jj,zz,ii,iscutmov,typei};
                            if ~isempty(ind)
                                gIX=ones(1,length(ind))';supervoxel_in_fish=supervolxeli(ind,1:3);clrmap = GetColormap('hsv_new',max(gIX));
                                for fraction_type=1:2
                                    switch fraction_type
                                        case 1
                                            temp_vol=supervolxeli(:,1:3);
                                        case 2
                                            temp_vol=temp_supervoxel;
                                    end
                                    [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervoxel_in_fish,reg_mask,reg_name,reg_loc,temp_vol,temp_env,clrmap,nn);
                                    Fraction_in_region_type_emerged{jj,zz,iscutmov,typei,fraction_type}(:,ii)=fraction_in_region_in_clust(:,1);
                                    Loc_in_region_emerged{jj,zz,ii,iscutmov,typei,fraction_type}=loc_in_region_in_clust;
                                    Loc_in_region_cell_emerged{jj,zz,ii,iscutmov,typei,fraction_type}=loc_in_region_cell;
                                    Num_in_region_emerged{jj,zz,iscutmov,typei,fraction_type,1}(:,ii)= num_in_region_in_clust{1}(1,:);
                                    Num_in_region_emerged{jj,zz,iscutmov,typei,fraction_type,2}(:,ii)= num_in_region_in_clust{2}(1,:);
                                    Num_in_region_emerged{jj,zz,iscutmov,typei,fraction_type,3}(:,ii)= num_in_region_in_clust{3}(1,:);
                                    Brain_region_id_emerged{jj,zz,ii,iscutmov,typei,fraction_type}=brain_region_id(:,1);
                                    Label_region_emerged{jj,zz,iscutmov,typei,fraction_type}(:,ii)=Label;
                                end
                            end
                        end
                    end
                end
            end
            
            %% mapback
            CS_response_all={CS_response{1,iscutmov},CS_response_acq{:,iscutmov},CS_response{2,iscutmov}};
            im=double(env.vol(:,:,1:end-1));rescalegd(im2double(env.vol), [1/10000 1/10000]);
            ismapback=false;
            if ismapback==true
                for jj=1:length(labels_typei)
                    h=figure('position',[10,10,3500,824]);t=tiledlayout(2,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.acq_block_num+2
                        ind=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        if ~isempty(ind)
                            nexttile(ii),
                            count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(sessionx(ii));
                            %cb = colorbar; cb.Layout.Tile = 'east';
                            nexttile(ii+trial.acq_block_num+2),
                            cr=squeeze(mean(mean(CS_response_all{:,ii}(:,:,ind),1,'omitnan'),2,'omitnan'));
                            mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(sessionx(ii));
                            %set(gca,'fontsize',16)
                        end
                    end
                    saveas(h,fullfile(checkpath(fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of ind_all_CSUS_RESPONSIVE'])),labels_typei{jj}),'jpeg');
                    close(h);
                end
                for zz=1:trial.acq_block_num+2
                    for jj=1:length(labels_typei)
                        for typei=1:size(ind_all_emergeCSUS,5)
                            h=figure('position',[7,154,3500,824]);t=tiledlayout(2,trial.acq_block_num+2);t.TileSpacing = 'compact';t.Padding = 'compact';
                            for ii=1:trial.acq_block_num+2
                                ind=ind_all_emergeCSUS{jj,zz,ii,iscutmov,typei};
                                if ~isempty(ind)
                                    nexttile(ii),
                                    count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(sessionx(ii));
                                    nexttile(ii+trial.acq_block_num+2),
                                    cr=squeeze(mean(mean(CS_response_all{:,ii}(:,:,ind),1,'omitnan'),2,'omitnan'));
                                    mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(sessionx(ii));
                                end
                                %set(gca,'fontsize',16)
                            end
                            saveas(h,fullfile(checkpath(fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of ind_all_emergeCSUS'],['Compare to session-',char(sessionx(zz))])),[labels_typei{jj},'-type',num2str(typei)]),'jpeg');
                            close(h);
                        end
                    end
                end
            end
            if false
                for zz=1:trial.acq_block_num+2
                    for jj=1:length(labels_typei)
                        if exist(fullfile(p,labels_typei{jj}),'dir')==7
                            try(rmdir(fullfile(p,labels_typei{jj}),'s'));
                            catch
                                warning(['fullfile(p,labels_typei{jj})',' already exist']);
                            end
                        end
                        for ii=[1:trial.acq_block_num+2]
                            ii
                            for typei=1:size(ind_all_emergeCSUS,5)
                                showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
                                id=ind_all_emergeCSUS{jj,zz,ii,iscutmov,typei};
                                clr=[1 0 0];
                                zi=unique(env.supervoxel(id,3));
                                for zz=zi'
                                    pti =  id(find(env.supervoxel(id,3)==zz));
                                    slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5))+3],'color',clr(kk,:));
                                    showspv(:,:,:,zz)=slice;
                                end
                                seqwrite(showspv,checkpath(char(fullfile(p,['iscutmov',num2str(iscutmov)],'Mapback of ind_all_emergeCSUS',['Compare to session-',char(sessionx(zz))],labels_typei{jj},['typei',num2str(typei)],char(sessionx(ii))))));
                            end
                        end
                    end
                end
                for jj=1:length(labels_typei)
                    if exist(fullfile(p,labels_typei{jj}),'dir')==7
                        try(rmdir(fullfile(p,labels_typei{jj}),'s'));
                        catch
                            warning(['fullfile(p,labels_typei{jj})',' already exist']);
                        end
                    end
                    for ii=[1:trial.acq_block_num+2]
                        ii
                        showspv=[];showspv(:,:,1,:)=im;showspv(:,:,2,:)=im;showspv(:,:,3,:)=im;showspv=showspv/255/2;
                        id=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        clr=[1 0 0];kk=1;
                        zi=unique(env.supervoxel(id,3));
                        for zz=zi'
                            pti =  id(find(env.supervoxel(id,3)==zz));
                            slice = insertShape(showspv(:,:,:,zz), 'filledcircle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5))+3],'color',clr(kk,:));
                            showspv(:,:,:,zz)=slice;
                        end
                        seqwrite(showspv, checkpath(char(fullfile(p,['iscutmov',num2str(iscutmov)],'Mapback of ind_all_CSUS_RESPONSIVE',labels_typei{jj},char(sessionx(ii))))));kk=kk+1;
                    end
                end
            end
            %% plot heatmap
            Colorregion=GetColormap('hsv_new',30);n=4;
            for typei=1:length((labels_typei))
                h=figure('position',[7,154,1910,824]);t=tiledlayout(1,(trial.acq_block_num+2)*(n+1));t.TileSpacing =  'none';t.Padding =  'none';
                for sessioni=1:trial.acq_block_num+2
                    ind_type=ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov};
                    region_id= double(Brain_region_id{typei,sessioni,iscutmov});
                    if ~isempty(region_id)
                        switch sessioni
                            case 1
                                ind_trial=trial.hab(2):trial.hab(3);
                            case trial.acq_block_num+2
                                ind_trial=trial.test(2):trial.test(3);
                            otherwise
                                ind_trial=(trial.acq(2):trial.acq(2)+trial.acq_block_trial-1)+(sessioni-2)*trial.acq_block_trial;
                        end
                        cr_trace=squeeze(mean(A_r(:,ind_trial,:),2));
                        ind_type_in_region_sort=[];
                        ind_type_in_region_sort_regionlable=[];
                        for regioni=1:length(Label_region{1,1}(:,1))%unique(region_id)'
                            ind_type_in_region=ind_type(find(region_id==regioni));
                            %cr_trace_in_region=cr_trace(:,ind_type_in_region);
                            cr=squeeze(mean(CS_response_all{:,sessioni}(:,:,ind_type_in_region),2,'omitnan'));
                            m=mean(cr,1);[B,I]=sort(m,'descend');
                            ind_type_in_region_sort=[ind_type_in_region_sort;ind_type_in_region(I)];
                            ind_type_in_region_sort_regionlable=[ind_type_in_region_sort_regionlable;repmat(regioni,length(ind_type_in_region),1)];
                            %nexttile(h);imagesc(cr_trace_in_region(:,I))
                        end
                        ind_type_sort{typei,sessioni,iscutmov}=ind_type_in_region_sort;
                        ind_type_sort_regionlabel{typei,sessioni,iscutmov}=ind_type_in_region_sort_regionlable;
                        a=length(ind_type_sort{typei,sessioni,iscutmov});
                    else
                        a=0;
                    end
                    if a~=0
                        nexttile([1,1]); axis ij;l=0;
                        for zz=unique(ind_type_sort_regionlabel{typei,sessioni,iscutmov})'
                            ll=length(find(ind_type_sort_regionlabel{typei,sessioni,iscutmov}==zz));
                            y = [l l+ll l+ll l];
                            x = [0 0 1 1];
                            l=l+ll;
                            patch(x,y,Colorregion(zz,:),'EdgeColor','none');hold on;
                        end
                        l=0;
                        for zz=unique(ind_type_sort_regionlabel{typei,sessioni,iscutmov})'
                            ll=length(find(ind_type_sort_regionlabel{typei,sessioni,iscutmov}==zz));
                            y = [l l+ll l+ll l];
                            x = [0 0 1 1];
                            text(0,(l+l+ll)/2,Label_region{typei,iscutmov}(zz,sessioni),'Rotation',0,'fontweight','bold','fontsize',8);hold on;
                            l=l+ll;
                        end
                        ylim([1 a]); xlim([0 1]);set(gca,'ytick',[],'xtick',[]);
                        ylabel(['Cell number:',num2str(length(ind_type_sort{typei,sessioni,iscutmov}))],'fontsize',12,'fontweight','bold');
                        nexttile([1,n]);
                        imagesc(cr_trace(:,ind_type_sort{typei,sessioni,iscutmov})',[0 0.15]);colormap('hot');colorbar('Location','northoutside');hold on
                        line([frame.cs_start frame.cs_start],[0 a],'Color','w','LineStyle','--','linewidth',2);hold on
                        line([frame.cs_end frame.cs_end],[0 a],'Color','w','LineStyle','--');hold on;
                        set(gca,'ytick',[]);%set(gca,'yticklabelrotation',90,'ytick',[0:ceil(a/400)*100:a,a]);
                        xlabel('Time (s)');
                    else
                        nexttile([1,1]);nexttile([1,n]);
                    end
                end
                saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],[labels_typei{typei},' heatmap(sort in region)']),'jpeg');
                close(h);
            end
            close all;
        end
    end
    %% SAVE
    save(fullfile(p,'CR_ind_summary_align.mat'),'A_r','trial_hab','trial_test',...
        'CR_ind_up','CR_ind_down','CR_ind_down_acq','CR_ind_up_acq',...
        'UR_ind',...
        'corr_CS_up','corr_CS_down',...
        'CS_response','CS_response_acq','US_response',...
        'CS_related_hab_ind','CS_related_hab_down_ind','CS_related_tst_ind','CS_related_tst_down_ind','CS_related_acq_ind','CS_related_acq_down_ind',...
        'NUM_UP_UR','NUM_UP_CR','NUM_DOWN_CR','NUM_NEW_EMERGED','labels_typei',...
        'num_of_ind_all_emergeCSUS',...
        'p_CS_up_bef_aft','p_CS_down_bef_aft','p_US_up_bef_cond','p_US_down_bef_cond','p_CS_up_bef_cond','p_CS_down_bef_cond',...
        'q_CS_up_bef_aft','q_CS_down_bef_aft','q_US_up_bef_cond','q_US_down_bef_cond','q_CS_up_bef_cond','q_CS_down_bef_cond',...
        'fdr_CS_up_bef_aft','fdr_CS_down_bef_aft','fdr_US_up_bef_cond','fdr_US_down_bef_cond','fdr_CS_up_bef_cond','fdr_CS_down_bef_cond',...
        'CS_up_ind','CS_down_ind','CS_stable_ind','CS_up','CS_down','CS_stable',...
        'US_up_ind','US_down_ind','US_stable_ind','US_up','US_down','US_stable',...
        'ind_all_CSUS_RESPONSIVE','Fraction_in_region_type','Loc_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region',...
        'ind_type_sort','ind_type_sort_regionlabel',...
        'ind_all_emergeCSUS','Fraction_in_region_type_emerged','Loc_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged',...
        '-v7.3');
end
% NUM_UP_CR:number of cs activation neuron in different sessions
% CS_up_ind:index of cs activation neuron with up amplitude in different sessions
% CS_up:number of cs activation neuron with up amplitude in different sessions
% ind_type_sort & ind_type_sort_regionlabel:不同type的neuron按脑区sort后的index以及region label
% ind_all_CSUS_RESPONSIVE :
%1：CS activation；2：CS inhibition；3：US activation；4：CS-up
%regulate；5：CS-down regulate；6：CS-stable regulate；7：US-up
%regulate；8：US-down regulate；9：US-stable regulate；10：CS-emerged（pre cond无，sessionx有）
% Fraction_in_region_type：不同type的neuron的index及脑区统计
% ind_all_emergeCSUS ：ind_all_CSUS_RESPONSIVE 不同类型的emerged neuron（type1：不存在在pre cond，存在在sessionx；type2：不存在在post cond，存在在sessionx；type1：不存在在pre 和post cond，存在在sessionx）
% 'num_of_ind_all_emergeCSUS':number of ind_all_CSUS_RESPONSIVE
%fraction type: 1:除这条鱼每个脑区分割出的神经元;2:除模板分割出的神经元
end
