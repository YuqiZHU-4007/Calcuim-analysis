function US_mapping_singletrial(Path)
set(0,'defaultfigurecolor','w')
addpath(genpath('F:\DUlab\FC_analyse\FishExplorer'));
if isempty(Path)
    p=uigetdir('H:\1.Test US\2.Tail free——Data from 117\')
else
    p=Path
end
if 1%~(exist([fullfile(p,'singletrial'),'\CR_ind_summary_singletrial.mat'],'file')==2)
    %     if (exist([fullfile(p,'singletrial'),'\CR_ind_summary_singletrial.mat'],'file')==2)
    %         delete([fullfile(p,'singletrial'),'\CR_ind_summary_singletrial.mat']);
    %     end
    if ~[exist(fullfile(p,'activities_aft_process.mat'),'file') && exist(fullfile(p,'para.mat'),'file') && exist(fullfile(p,'env.mat'),'file')]
        warning([p '_file not exist']);
    else
        load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'env.mat'));
        p=checkpath(fullfile(p,'singletrial'));
        A=activities_preCS_dfdf_aftcorrect;%A=zscore(A,0,'all');
        A_r=reshape(A,frame.per_cycle,trial.total,[]);
        
        res=[0.66,0.66,10];
        v=[90,90];
        radium=floor(15/0.66);
        load('X:\calcium data 20230224\脑区分割\segmentation_file_0525_DSregion_mask.mat');
        temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');
        warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
        temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
        nn=[Path(end-14:end-7),Path(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        
        trial_hab=trial.hab(2):trial.hab(3);
        trial_test=trial.test(2):trial.test(3);
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
        cut_criterion=[3,2,3];% CS activation;CS inhibition;US activation
        %% 有无CR on in hab & tst & acq
        pre_cs_time=[6 0];
        base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
        cs_win=frame.cs_start:frame.us_start-2;
        CRon_activation_ind=nan(size(A,2),trial.total,2);CRon_inhibition_ind=nan(size(A,2),trial.total,2);CSon_activation_response=cell(trial.total,2);
        trial_ind_base=trial_hab;
        aa=A_r(base_win,trial_ind_base,:);aa=squeeze(mean(aa,2,'omitnan'));
        m=mean(aa,1,'omitnan');sd=std(aa,[],1,'omitnan');
        baseline_up=m+ cut_criterion(1)*sd; baseline_down=m- cut_criterion(2)*sd;baseline_UR=m+cut_criterion(3)*sd;
        for ii=1:trial.total
            ind=ii;
            aa=A_r(:,ind,:);
            for iscutmov=1 %1:no_cut;2:cut
                b=nan(length(cs_win),length(ind),size(A,2));
                for jj=1:length(ind)
                    b(:,jj,:)=aa(cs_win,jj,:);
                end
                m_CS=squeeze(mean(mean(b,1,'omitnan'),2,'omitnan'))';
                CSon_activation_response{ii,iscutmov}(:,:,:)=b;
                CRon_activation_ind(:,ii,iscutmov)=m_CS > baseline_up;
                CRon_inhibition_ind(:,ii,iscutmov)=m_CS < baseline_down;
            end
        end
        %% 有无CR off
        cs_win=frame.us_start:frame.us_start+2/fs.ca;
        CRoff_activation_ind=nan(size(A,2),trial.total,2);CRoff_inhibition_ind=nan(size(A,2),trial.total,2);CSoff_activation_response=cell(trial.total,2);
        for ii=1:trial.total
            ind=ii;
            aa=A_r(:,ind,:);
            for iscutmov=1 %1:no_cut;2:cut
                b=nan(length(cs_win),length(ind),size(A,2));
                for jj=1:length(ind)
                    b(:,jj,:)=aa(cs_win,jj,:);
                end
                m_CS=squeeze(mean(mean(b,1,'omitnan'),2,'omitnan'))';
                CSoff_activation_response{ii,iscutmov}(:,:,:)=b;
                CRoff_activation_ind(:,ii,iscutmov)=m_CS > baseline_up;
                CRoff_inhibition_ind(:,ii,iscutmov)=m_CS < baseline_down;
            end
        end
        %% 有无UR
        us_win=frame.us_start+1:frame.us_start+2/fs.ca;
        UR_activation_ind=nan(size(A,2),trial.total);US_activation_response=cell(trial.total,1);
        for ii=trial.acq(2):trial.acq(3)
            ind=ii;
            aa=squeeze(mean(A_r(:,ind,:),2));
            m_US=mean(aa(us_win,:),1);
            UR_activation_ind(:,ii)=m_US> baseline_UR;
            US_activation_response{ii}=A_r(us_win,ind,:);
        end
        %% 有无UR cut trail
        us_win=frame.us_start:frame.us_start+2/fs.ca;
        cs_win=frame.cs_start:frame.us_start-2;
        UR_activation_ind_cuttrail=nan(size(A,2),trial.total);US_activation_response_cuttrail=cell(trial.total,1);
        for ii=trial.acq(2):trial.acq(3)
            ind=ii;
            aa=squeeze(mean(A_r(:,ind,:),2));
            m_US=mean(aa(us_win,:),1)-mean(aa(cs_win,:),1);
            UR_activation_ind_cuttrail(:,ii)=m_US> baseline_UR;
            US_activation_response_cuttrail{ii}=A_r(us_win,ind,:);
        end
        close all
        %% onset
        onset_trial=[];
        onset_win=frame.cs_start-4:frame.per_cycle;
        for ii=1:trial.total
            [onset_trial(:,ii),~,~,~]=findonset( squeeze(A_r(:,ii,:)),[floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1)],onset_win);
        end
        
        %% part of plot from mapping_index_plot.m
        labels_typei={'CSon-avtivation','CSon-inhibition','CSoff-avtivation','CSoff-inhibition',...
            'US-activation','US-activation-cuttrail',...
            'CSon-up regulate','CSon-down regulate','CSon-stable regulate',...
            'CSoff-up regulate','CSoff-down regulate','CSoff-stable regulate',...
            'US-up regulate','US-down regulate','US-stable regulate',...,
            'US-up regulate-cuttrail','US-down regulate-cuttrail','US-stable regulate-cuttrail',...,
            'CSon-up regulate-all','CSon-down regulate-all','CSon-stable regulate-all',...
            'CSoff-up regulate-all','CSoff-down regulate-all','CSoff-stable regulate-all',...
            'US-up regulate-all','US-down regulate-all','US-stable regulate-all',...,
            'US-up regulate-cuttrail-all','US-down regulate-cuttrail-all','US-stable regulate-cuttrail-all',...
            'CSUS-shift','CSUS-shift-cuttrail','CSUS-shift-all','CSUS-shift-cuttrail-all'};
        if isfile(fullfile(checkpath(p),'CR_ind_summary_singletrial_statistic.mat'))
            ranksumtest=0;
        else
            p_CSon_up=nan(size(A_r,3),trial.total,2);p_CSon_down=nan(size(A_r,3),trial.total,2);
            p_CSoff_up=nan(size(A_r,3),trial.total,2);p_CSoff_down=nan(size(A_r,3),trial.total,2);
            p_US_up=nan(size(A_r,3),trial.total);p_US_down=nan(size(A_r,3),trial.total);
            p_US_up_cuttrail=nan(size(A_r,3),trial.total);p_US_down_cuttrail=nan(size(A_r,3),trial.total);
            q_CSon_up=p_CSon_up;q_CSon_down=p_CSon_down;q_US_up=p_US_up;q_US_down=p_US_down;
            q_CSoff_up=p_CSoff_up;q_CSoff_down=p_CSoff_down;q_US_up_cuttrail=p_US_up_cuttrail;q_US_down_cuttrail=p_US_down_cuttrail;
            fdr_CSon_up=p_CSon_up;fdr_CSon_down=p_CSon_down;fdr_US_up=p_US_up;fdr_US_down=p_US_down;
            fdr_CSoff_up=p_CSoff_up;fdr_CSoff_down=p_CSoff_down;fdr_US_up_cuttrail=p_US_up_cuttrail;fdr_US_down_cuttrail=p_US_down_cuttrail;
            ranksumtest=1;
        end
        CSon_up_ind=cell(trial.total,2); CSon_down_ind=cell(trial.total,2); CSon_stable_ind=cell(trial.total,2);
        CSon_up_ind_all=cell(trial.total,2); CSon_down_ind_all=cell(trial.total,2); CSon_stable_ind_all=cell(trial.total,2);
        CSon_up=nan(trial.total,2);CSon_down=nan(trial.total,2);CSon_stable=nan(trial.total,2);
        CSoff_up_ind=cell(trial.total,2); CSoff_down_ind=cell(trial.total,2); CSoff_stable_ind=cell(trial.total,2);
        CSoff_up_ind_all=cell(trial.total,2); CSoff_down_ind_all=cell(trial.total,2); CSoff_stable_ind_all=cell(trial.total,2);
        CSoff_up=nan(trial.total,2);CSoff_down=nan(trial.total,2);CSoff_stable=nan(trial.total,2);
        US_up_ind=cell(trial.total,2); US_down_ind=cell(trial.total,2); US_stable_ind=cell(trial.total,2);
        US_up_ind_all=cell(trial.total,2); US_down_ind_all=cell(trial.total,2); US_stable_ind_all=cell(trial.total,2);
        US_up=nan(trial.total,2);US_down=nan(trial.total,2);US_stable=nan(trial.total,2);
        US_up_ind_cuttrail=cell(trial.total,2); US_down_ind_cuttrail=cell(trial.total,2); US_stable_ind_cuttrail=cell(trial.total,2);
        US_up_ind_all_cuttrail=cell(trial.total,2); US_down_ind_all_cuttrail=cell(trial.total,2); US_stable_ind_all_cuttrail=cell(trial.total,2);
        US_up_cuttrail=nan(trial.total,2);US_down_cuttrail=nan(trial.total,2);US_stable_cuttrail=nan(trial.total,2);
        CS_US_shift_ind=cell(trial.total,2); CS_US_shift_ind_cuttrail=cell(trial.total,2);
        CS_US_shift_ind_all=cell(trial.total,2); CS_US_shift_ind_all_cuttrail=cell(trial.total,2);
        
        ind_all_CSUS_RESPONSIVE=cell(length(labels_typei),trial.total,2);
        ind_all_emergeCSUS_inthissession=cell(length(labels_typei),trial.total,2,3);
        ind_all_emergeCSUS=cell(length(labels_typei),trial.total,trial.total,2,1);
       
        Response_all_CSUS_RESPONSIVE=nan(size(A_r,1),trial.total,length(labels_typei));
        Response_all_emergeCSUS_inthissession=nan(size(A_r,1),trial.total,length(labels_typei),3);
        Response_all_emergeCSUS=nan(size(A_r,1),trial.total,length(labels_typei),1);
        CS_all_tst=cell(length(labels_typei),trial.total,2); CS_all_hab=cell(length(labels_typei),trial.total,2);
        
        Fraction_in_region_type=cell(length(labels_typei),2,2);Index_in_region=cell(length(labels_typei),trial.total,2,2);
        Loc_in_region_cell=cell(length(labels_typei),trial.total,2,2);Num_in_region=cell(length(labels_typei),2,2,3);
        Brain_region_id=cell(length(labels_typei),trial.total,2,2);Label_region=cell(length(labels_typei),2,2);
        
        Fraction_in_region_type_emerged=cell(length(labels_typei),2,2,1);Index_in_region_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);
        Loc_in_region_cell_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);Num_in_region_emerged=cell(length(labels_typei),2,2,3,1);
        Brain_region_id_emerged=cell(length(labels_typei),trial.total,trial.total,2,2,1);Label_region_emerged=cell(length(labels_typei),2,2,1);
        
        Fraction_in_region_type_emerged_inthissession=cell(length(labels_typei),2,2,3);Index_in_region_emerged_inthissession=cell(length(labels_typei),trial.total,2,3);
        Loc_in_region_cell_emerged_inthissession=cell(length(labels_typei),trial.total,2,2,3);Num_in_region_emerged_inthissession=cell(length(labels_typei),2,2,3,3);
        Brain_region_id_emerged_inthissession=cell(length(labels_typei),trial.total,2,2,3);Label_region_emerged_inthissession=cell(length(labels_typei),2,2,3);


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
            %% UP/DOWN CS/US-responsive neuron
            %down/up-regulate
            if ranksumtest==1;
                x={CSon_activation_response{:,iscutmov}};a=[];
                for ii=trial.hab(2):trial.hab(3)
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A,2)
                    for jj=1:trial.total
                        p_CSon_up(ii,jj,iscutmov)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                        p_CSon_down(ii,jj,iscutmov)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','right');
                    end
                end
                x={CSoff_activation_response{:,iscutmov}};a=[];
                for ii=trial.hab(2):trial.hab(3)
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A,2)
                    for jj=1:trial.total
                        p_CSoff_up(ii,jj,iscutmov)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                        p_CSoff_down(ii,jj,iscutmov)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','right');
                    end
                end
                x=US_activation_response;
                a=[];
                for ii=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A,2)
                    for jj=trial.acq(2):trial.acq(3)
                        p_US_up(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                        p_US_down(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                    end
                end
                x=US_activation_response_cuttrail;
                a=[];
                for ii=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A,2)
                    for jj=trial.acq(2):trial.acq(3)
                        p_US_up_cuttrail(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                        p_US_down_cuttrail(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','left');
                    end
                end
                for jj=1:trial.total
                    [fdr_CSon_up(:,jj,iscutmov),q_CSon_up(:,jj,iscutmov)] = FDR_base(p_CSon_up(:,jj,iscutmov),'BH');
                    [fdr_CSon_down(:,jj,iscutmov),q_CSon_down(:,jj,iscutmov)] = FDR_base(p_CSon_down(:,jj,iscutmov),'BH');
                    [fdr_CSoff_up(:,jj,iscutmov),q_CSoff_up(:,jj,iscutmov)] = FDR_base(p_CSoff_up(:,jj,iscutmov),'BH');
                    [fdr_CSoff_down(:,jj,iscutmov),q_CSoff_down(:,jj,iscutmov)] = FDR_base(p_CSoff_down(:,jj,iscutmov),'BH');
                    [fdr_US_up(:,jj,iscutmov),q_US_up(:,jj,iscutmov)] = FDR_base(p_US_up(:,jj,iscutmov),'BH');
                    [fdr_US_down(:,jj,iscutmov),q_US_down(:,jj,iscutmov)] = FDR_base(p_US_down(:,jj,iscutmov),'BH');
                    [fdr_US_up_cuttrail(:,jj,iscutmov),q_US_up_cuttrail(:,jj,iscutmov)] = FDR_base(p_US_up_cuttrail(:,jj,iscutmov),'BH');
                    [fdr_US_down_cuttrail(:,jj,iscutmov),q_US_down_cuttrail(:,jj,iscutmov)] = FDR_base(p_US_down_cuttrail(:,jj,iscutmov),'BH');
                end
                save(fullfile(checkpath(p),'CR_ind_summary_singletrial_statistic.mat'),'p_CSon_up','p_CSon_down','q_CSon_up','q_CSon_down','fdr_CSon_up','fdr_CSon_down',...
                    'p_CSoff_up','p_CSoff_down','q_CSoff_up','q_CSoff_down','fdr_CSoff_up','fdr_CSoff_down',...
                    'p_US_up','p_US_down','q_US_up','q_US_down','fdr_US_up','fdr_US_down',...
                    'p_US_up_cuttrail','p_US_down_cuttrail','q_US_up_cuttrail','q_US_down_cuttrail','fdr_US_up_cuttrail','fdr_US_down_cuttrail','-v7.3');
            else
                load(fullfile(checkpath(p),'CR_ind_summary_singletrial_statistic.mat'),'p_CSon_up','p_CSon_down','p_CSoff_up','p_CSoff_down', 'p_US_up','p_US_down', 'p_US_up_cuttrail','p_US_down_cuttrail');
            end
            for ii=1:trial.total
                CSon_up_ind{ii,iscutmov}=find(squeeze(p_CSon_up(:,ii,iscutmov))<=0.05 & squeeze(sum(CRon_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSon_down_ind{ii,iscutmov}=find(squeeze(p_CSon_down(:,ii,iscutmov))<=0.05 & squeeze(sum(CRon_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSon_stable_ind{ii,iscutmov}=find((squeeze(p_CSon_up(:,ii,iscutmov))>0.05 & squeeze(p_CSon_down(:,ii,iscutmov))>0.05) & squeeze(sum(CRon_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSon_up_ind_all{ii,iscutmov}=find(squeeze(p_CSon_up(:,ii,iscutmov))<=0.05 );
                CSon_down_ind_all{ii,iscutmov}=find(squeeze(p_CSon_down(:,ii,iscutmov))<=0.05 );
                CSon_stable_ind_all{ii,iscutmov}=find((squeeze(p_CSon_up(:,ii,iscutmov))>0.05 & squeeze(p_CSon_down(:,ii,iscutmov))>0.05));
                
                CSoff_up_ind{ii,iscutmov}=find(squeeze(p_CSoff_up(:,ii,iscutmov))<=0.05 & squeeze(sum(CRoff_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSoff_down_ind{ii,iscutmov}=find(squeeze(p_CSoff_down(:,ii,iscutmov))<=0.05 & squeeze(sum(CRoff_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSoff_stable_ind{ii,iscutmov}=find((squeeze(p_CSoff_up(:,ii,iscutmov))>0.05 & squeeze(p_CSoff_down(:,ii,iscutmov))>0.05) & squeeze(sum(CRoff_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1);
                CSoff_up_ind_all{ii,iscutmov}=find(squeeze(p_CSoff_up(:,ii,iscutmov))<=0.05 );
                CSoff_down_ind_all{ii,iscutmov}=find(squeeze(p_CSoff_down(:,ii,iscutmov))<=0.05 );
                CSoff_stable_ind_all{ii,iscutmov}=find((squeeze(p_CSoff_up(:,ii,iscutmov))>0.05 & squeeze(p_CSoff_down(:,ii,iscutmov))>0.05));
                
                US_up_ind{ii,iscutmov}=find(squeeze(p_US_up(:,ii))<=0.05 & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1);
                US_down_ind{ii,iscutmov}=find(squeeze(p_US_down(:,ii))<=0.05 & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1);
                US_stable_ind{ii,iscutmov}=find((squeeze(p_US_up(:,ii))>0.05 & squeeze(p_US_down(:,ii))>0.05) & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1);
                US_up_ind_all{ii,iscutmov}=find(squeeze(p_US_up(:,ii))<=0.05);
                US_down_ind_all{ii,iscutmov}=find(squeeze(p_US_down(:,ii))<=0.05);
                US_stable_ind_all{ii,iscutmov}=find((squeeze(p_US_up(:,ii))>0.05 & squeeze(p_US_down(:,ii))>0.05));
                
                US_up_ind_cuttrail{ii,iscutmov}=find(squeeze(p_US_up_cuttrail(:,ii))<=0.05 & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1);
                US_down_ind_cuttrail{ii,iscutmov}=find(squeeze(p_US_down_cuttrail(:,ii))<=0.05 & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1);
                US_stable_ind_cuttrail{ii,iscutmov}=find((squeeze(p_US_up_cuttrail(:,ii))>0.05 & squeeze(p_US_down_cuttrail(:,ii))>0.05) & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1);
                US_up_ind_all_cuttrail{ii,iscutmov}=find(squeeze(p_US_up_cuttrail(:,ii))<=0.05);
                US_down_ind_all_cuttrail{ii,iscutmov}=find(squeeze(p_US_down_cuttrail(:,ii))<=0.05);
                US_stable_ind_all_cuttrail{ii,iscutmov}=find((squeeze(p_US_up_cuttrail(:,ii))>0.05 & squeeze(p_US_down_cuttrail(:,ii))>0.05));
                
                a=[squeeze(p_CSon_up(:,ii,iscutmov))<=0.05 & squeeze(sum(CRon_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1];
                b=[squeeze(p_US_down(:,ii))<=0.05 & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1];
                CS_US_shift_ind{ii,iscutmov}=find((a & b) == 1);
                a=[squeeze(p_CSon_up(:,ii,iscutmov))<=0.05];
                b=[squeeze(p_US_down(:,ii))<=0.05 ];
                CS_US_shift_ind_all{ii,iscutmov}=find((a & b) == 1);
                
                a=[squeeze(p_CSon_up(:,ii,iscutmov))<=0.05 & squeeze(sum(CRon_activation_ind(:,[trial_hab,ii],iscutmov),2))>=1];
                b=[squeeze(p_US_down_cuttrail(:,ii))<=0.05 & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1];
                CS_US_shift_ind_cuttrail{ii,iscutmov}=find((a & b) == 1);
                a=[squeeze(p_CSon_up(:,ii,iscutmov))<=0.05];
                b=[squeeze(p_US_down_cuttrail(:,ii))<=0.05];
                CS_US_shift_ind_all_cuttrail{ii,iscutmov}=find((a & b) == 1);
            end
            if ranksumtest==1;
                clear('p_CSon_up','p_CSon_down','q_CSon_up','q_CSon_down','fdr_CSon_up','fdr_CSon_down',...
                    'p_CSoff_up','p_CSoff_down','q_CSoff_up','q_CSoff_down','fdr_CSoff_up','fdr_CSoff_down',...
                    'p_US_up','p_US_down','q_US_up','q_US_down','fdr_US_up','fdr_US_down',...
                    'p_US_up_cuttrail','p_US_down_cuttrail','q_US_up_cuttrail','q_US_down_cuttrail','fdr_US_up_cuttrail','fdr_US_down_cuttrail');
            else
                clear('p_CSon_up','p_CSon_down','p_CSoff_up','p_CSoff_down', 'p_US_up','p_US_down', 'p_US_up_cuttrail','p_US_down_cuttrail');
                
            end
            
            
            %% 脑区统计
            %cs-responsive
            %us-responsive
            %cs/us up/down/stable regulate
            for jj=1:length(labels_typei)
                for ii=trial_hab
                    CS_all_hab{jj,ii,iscutmov}=union(CS_all_hab{jj,ii,iscutmov},ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov});
                end
                for ii=trial_test
                    CS_all_tst{jj,ii,iscutmov}=union(CS_all_tst{jj,ii,iscutmov},ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov});
                end
            end
            for ii=1:trial.total
                ind_all_CSUS_RESPONSIVE{1,ii,iscutmov}=find(CRon_activation_ind(:,ii,iscutmov)==1);
                ind_all_CSUS_RESPONSIVE{2,ii,iscutmov}=find(CRon_inhibition_ind(:,ii,iscutmov)==1);
                ind_all_CSUS_RESPONSIVE{3,ii,iscutmov}=find(CRoff_activation_ind(:,ii,iscutmov)==1);
                ind_all_CSUS_RESPONSIVE{4,ii,iscutmov}=find(CRoff_inhibition_ind(:,ii,iscutmov)==1);
                ind_all_CSUS_RESPONSIVE{5,ii,iscutmov}=find(UR_activation_ind(:,ii)==1);
                ind_all_CSUS_RESPONSIVE{6,ii,iscutmov}=find(UR_activation_ind_cuttrail(:,ii)==1);
                ind_all_CSUS_RESPONSIVE{7,ii,iscutmov}=CSon_up_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{8,ii,iscutmov}=CSon_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{9,ii,iscutmov}=CSon_stable_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{10,ii,iscutmov}=CSoff_up_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{11,ii,iscutmov}=CSoff_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{12,ii,iscutmov}=CSoff_stable_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{13,ii,iscutmov}=US_up_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{14,ii,iscutmov}=US_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{15,ii,iscutmov}=US_stable_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{16,ii,iscutmov}=US_up_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{17,ii,iscutmov}=US_down_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{18,ii,iscutmov}=US_stable_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{19,ii,iscutmov}=CSon_up_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{20,ii,iscutmov}=CSon_down_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{21,ii,iscutmov}=CSon_stable_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{22,ii,iscutmov}=CSoff_up_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{23,ii,iscutmov}=CSoff_down_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{24,ii,iscutmov}=CSoff_stable_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{25,ii,iscutmov}=US_up_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{26,ii,iscutmov}=US_down_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{27,ii,iscutmov}=US_stable_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{28,ii,iscutmov}=US_up_ind_all_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{29,ii,iscutmov}=US_down_ind_all_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{30,ii,iscutmov}=US_stable_ind_all_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{31,ii,iscutmov}=CS_US_shift_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{32,ii,iscutmov}=CS_US_shift_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{33,ii,iscutmov}=CS_US_shift_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{34,ii,iscutmov}=CS_US_shift_ind_cuttrail{ii,iscutmov};
            end
            for jj=1:length(labels_typei)
                for ii=1:1:trial.total
                    for zz=1:trial.total
                        a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        b=ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov};
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}=setdiff(a,b);
                    end
                end
            end
            for jj=1:length(labels_typei)
                for ii=1:1:trial.total
                    b=[];
                    for zz=1:ii-1
                        b=union(b,ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov});
                    end
                    a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                    ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,1}=setdiff(a,b);
                    ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,2}=setdiff(a,CS_all_hab{jj,ii,iscutmov});
                    ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,3}=setdiff(a,CS_all_tst{jj,ii,iscutmov});
                end
            end
            for jj=1:length(labels_typei)
                for ii=1:1:trial.total
                    for zz=1:trial.total
                        a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        b=ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov};
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}=setdiff(a,b);
                    end
                end
            end
            for jj=1:length(labels_typei)
                for ii=1:trial.total
                    Response_all_CSUS_RESPONSIVE(:,ii,jj)=squeeze(mean(A_r(:,ii,ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov}),3,'omitnan'));
                    Response_all_emergeCSUS_inthissession(:,ii,jj,1)=squeeze(mean(A_r(:,ii,ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,1}),3,'omitnan'));
                    Response_all_emergeCSUS_inthissession(:,ii,jj,2)=squeeze(mean(A_r(:,ii,ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,2}),3,'omitnan'));
                    Response_all_emergeCSUS_inthissession(:,ii,jj,3)=squeeze(mean(A_r(:,ii,ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,3}),3,'omitnan'));
                    for zz=1:trial.total
                        Response_all_emergeCSUS(:,ii,zz,jj,1)=squeeze(mean(A_r(:,ii,ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}),3,'omitnan'));
                    end
                end
            end
            save(fullfile(p,'CR_ind_summary_singletrial.mat'),'A_r','trial_hab','trial_test',...
                'cut_criterion','cs_win','us_win','base_win','pre_cs_time','CRon_activation_ind','CRon_inhibition_ind','CSon_activation_response',...
                'CRoff_activation_ind','CRoff_inhibition_ind','CSoff_activation_response',...
                'UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail',...
                'onset_trial',...
                'CS_all_hab','CS_all_tst',...
                'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS','ind_all_emergeCSUS_inthissession',...
                'CSon_up_ind','CSon_down_ind','CSon_stable_ind','CSon_up_ind_all','CSon_down_ind_all','CSon_stable_ind_all',...
                'CSoff_up_ind','CSoff_down_ind','CSoff_stable_ind','CSoff_up_ind_all','CSoff_down_ind_all','CSoff_stable_ind_all',...
                'US_up_ind','US_down_ind','US_stable_ind','US_up_ind_all','US_down_ind_all','US_stable_ind_all',...
                'US_up_ind_cuttrail','US_down_ind_cuttrail','US_stable_ind_cuttrail','US_up_ind_all_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
                'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
                'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS','-v7.3');
            clear('A_r','trial_hab','trial_test',...
                'cut_criterion','cs_win','us_win','base_win','pre_cs_time','CRon_activation_ind','CRon_inhibition_ind','CSon_activation_response',...
                'CRoff_activation_ind','CRoff_inhibition_ind','CSoff_activation_response',...
                'UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail',...
                'onset_trial',...
                'CS_all_hab','CS_all_tst',...
                'CSon_up_ind','CSon_down_ind','CSon_stable_ind','CSon_up_ind_all','CSon_down_ind_all','CSon_stable_ind_all',...
                'CSoff_up_ind','CSoff_down_ind','CSoff_stable_ind','CSoff_up_ind_all','CSoff_down_ind_all','CSoff_stable_ind_all',...
                'US_up_ind','US_down_ind','US_stable_ind','US_up_ind_all','US_down_ind_all','US_stable_ind_all',...
                'US_up_ind_cuttrail','US_down_ind_cuttrail','US_stable_ind_cuttrail','US_up_ind_all_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
                'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
                'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS');
            
            %fraction
            fraction_type=1;
            for jj=1:length(labels_typei)
                for ii=1:trial.total
                    ind_type=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                    Brain_region_id{jj,ii,iscutmov,fraction_type}=brain_region_id(ind_type,:);
                    for regioni=1:length(Label)
                        if length(index_in_region_in_clust{regioni,1})<=10
                            ind_inregion=[];
                        else
                            ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                        end
                        Index_in_region{jj,ii,iscutmov,fraction_type}{regioni,:}=ind_inregion;
                        Fraction_in_region_type{jj,iscutmov,fraction_type}(regioni,ii)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                        Loc_in_region_cell{jj,ii,iscutmov,fraction_type}{regioni,:}=supervolxeli(ind_inregion,:);
                        Num_in_region{jj,iscutmov,fraction_type,1}(regioni,ii)= length(ind_inregion)./length(ind_type);
                        Num_in_region{jj,iscutmov,fraction_type,2}(regioni,ii)= length(ind_inregion);
                        Num_in_region{jj,iscutmov,fraction_type,3}(regioni,ii)= length(ind_type);
                        Label_region{jj,iscutmov,fraction_type}(regioni,:)=Label{regioni};
                    end
                end
            end
            for jj=1:length(labels_typei);%[1:18,31:34]
                for typei=1:size(ind_all_emergeCSUS_inthissession,4)
                    for ii=1:trial.total
                        ind_type=ind_all_emergeCSUS_inthissession{jj,ii,iscutmov,typei};
                        Brain_region_id_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}=brain_region_id(ind_type,:);
                        for regioni=1:length(Label)
                            if length(index_in_region_in_clust{regioni,1})<=10
                                ind_inregion=[];
                            else
                                ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                            end
                            Index_in_region_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}{regioni,:}=ind_inregion;
                            Fraction_in_region_type_emerged_inthissession{jj,iscutmov,fraction_type,typei}(regioni,ii)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                            Loc_in_region_cell_emerged_inthissession{jj,ii,iscutmov,fraction_type,typei}{regioni,:}=supervolxeli(ind_inregion,:);
                            Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,1,typei}(regioni,ii)= length(ind_inregion)./length(ind_type);
                            Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,2,typei}(regioni,ii)= length(ind_inregion);
                            Num_in_region_emerged_inthissession{jj,iscutmov,fraction_type,3,typei}(regioni,ii)= length(ind_type);
                            Label_region_emerged_inthissession{jj,iscutmov,fraction_type,typei}(regioni,:)=Label{regioni};
                        end
                    end
                end
            end
            for jj=1:length(labels_typei);%[1:18,31:34]
                for  typei=1:size(ind_all_emergeCSUS,5)
                    for ii=1:trial.total
                        for zz=1:trial.total
                            ind_type=ind_all_emergeCSUS{jj,ii,zz,iscutmov,typei};
                            Brain_region_id_emerged{jj,ii,zz,iscutmov,fraction_type,typei}=brain_region_id(ind_type,:);
                            for regioni=1:length(Label)
                                if length(index_in_region_in_clust{regioni,1})<=10
                                    ind_inregion=[];
                                else
                                    ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                                end
                                Index_in_region_emerged{jj,ii,zz,iscutmov,fraction_type,typei}{regioni,:}=ind_inregion;
                                Fraction_in_region_type_emerged{jj,iscutmov,fraction_type,typei}(regioni,ii,zz)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                                Loc_in_region_cell_emerged{jj,ii,zz,iscutmov,fraction_type,typei}{regioni,:}=supervolxeli(ind_inregion,:);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,1,typei}(regioni,ii,zz)= length(ind_inregion)./length(ind_type);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,2,typei}(regioni,ii,zz)= length(ind_inregion);
                                Num_in_region_emerged{jj,iscutmov,fraction_type,3,typei}(regioni,ii,zz)= length(ind_type);
                                Label_region_emerged{jj,iscutmov,fraction_type,typei}(regioni,:)=Label{regioni};
                            end
                        end
                    end
                end
            end
            save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'ind_all_emergeCSUS','Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region',...
                'Fraction_in_region_type_emerged_inthissession','Index_in_region_emerged_inthissession','Loc_in_region_cell_emerged_inthissession','Num_in_region_emerged_inthissession','Brain_region_id_emerged_inthissession','Label_region_emerged_inthissession',...
                'Fraction_in_region_type_emerged','Index_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged',...
                '-append','-v7.3');
            %% mapback
            for ii=[]
                supervoxel_raw=env.supervoxel;
                for ii=unique(supervoxel_raw(:,3))'
                    id=find(supervoxel_raw(:,3)==ii);
                    supervoxel_raw(id,3)=(ii-1)*10/0.66+1;
                end
                supervoxel_raw(:,1)=env.width-supervoxel_raw(:,1);

                CS_response_all={CSon_activation_response{:,iscutmov}};
                im=double(env.vol(:,:,1:end-1));rescalegd(im2double(env.vol), [1/10000 1/10000]);
                ismapback=false;
                if ismapback==true
                    for jj=1:length(labels_typei)
                        h=figure('position',[10,10,3500,824]);t=tiledlayout(2,trial.total);t.TileSpacing = 'compact';t.Padding = 'compact';
                        for ii=1:trial.total
                            ind=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                            if ~isempty(ind)
                                nexttile(ii),
                                count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(sessionx(ii));
                                %cb = colorbar; cb.Layout.Tile = 'east';
                                nexttile(ii+trial.total),
                                cr=squeeze(mean(mean(CS_response_all{:,ii}(:,:,ind),1,'omitnan'),2,'omitnan'));
                                mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(sessionx(ii));
                                %set(gca,'fontsize',16)
                            end
                        end
                        saveas(h,fullfile(checkpath(fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of ind_all_CSUS_RESPONSIVE'])),labels_typei{jj}),'jpeg');
                        close(h);
                    end
                    for zz=1:trial.total
                        for jj=1:length(labels_typei)
                            for typei=1:size(ind_all_emergeCSUS,5)
                                h=figure('position',[7,154,3500,824]);t=tiledlayout(2,trial.total);t.TileSpacing = 'compact';t.Padding = 'compact';
                                for ii=1:trial.total
                                    ind=ind_all_emergeCSUS{jj,zz,ii,iscutmov,typei};
                                    if ~isempty(ind)
                                        nexttile(ii),
                                        count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(sessionx(ii));
                                        nexttile(ii+trial.total),
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
                %% plot heatmap
                for ii=[]
                    Colorregion=GetColormap('hsv_new',30);n=4;
                    for typei=1:length((labels_typei))
                        h=figure('position',[7,154,1910,824]);t=tiledlayout(4,(trial.total)*(n+1)/4);t.TileSpacing =  'none';t.Padding =  'none';
                        for sessioni=1:trial.total
                            ind_type=ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov};
                            region_id= double(Brain_region_id{typei,sessioni,iscutmov});
                            if ~isempty(region_id)
                                ind_trial=sessioni;
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
                                imagesc(cr_trace(:,ind_type_sort{typei,sessioni,iscutmov})',[0 0.8]);colormap('hot');colorbar('Location','northoutside');hold on
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
                end
                close all;
                %% plot spatial location
                for ii=[]
                    %up CR in different session
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(CRon_activation_ind(:,ii,iscutmov)==1);leg=[num2str(ii)];
                        NUM_UP_CR(ii,iscutmov)=length(ind);
                        nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
                        set(gca,'fontsize',16)
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of CS activation neuron']),'jpeg');
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(CRon_activation_ind(:,ii,iscutmov)==1);
                        cr=squeeze(mean(mean(CSon_activation_response{ii,iscutmov}(:,:,ind),1,'omitnan'),2,'omitnan'));leg=[num2str(ii)];
                        nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(leg);
                        set(gca,'fontsize',16)
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['CR AMP of CS activation neuron']),'jpeg');
                    %down CR in different session
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(CRon_inhibition_ind(:,ii,iscutmov)==1);leg=[num2str(ii)];
                        NUM_DOWN_CR(ii,iscutmov)=length(ind);
                        nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
                        set(gca,'fontsize',16)
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of CS inhibition neuron']),'jpeg');
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(CRon_inhibition_ind(:,1,iscutmov)==1);leg=[num2str(ii)];
                        nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),cr);view(v(1),v(2));title(leg);
                        set(gca,'fontsize',16)
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['CR AMP of CS inhibition neuron']),'jpeg');
                    %UR in different session
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(UR_activation_ind(:,ii)==1);
                        NUM_UP_UR(ii)=length(ind);
                        nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));
                        set(gca,'fontsize',16);title(['Acq.' num2str(ii)]);
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Spatial_location_density of US responsive neuron']),'jpeg');
                    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((trial.total)/5),5);t.TileSpacing = 'compact';t.Padding = 'compact';
                    for ii=1:trial.total
                        ind=find(UR_activation_ind(:,ii)==1);
                        nexttile,mapback_colorcode_responses(supervoxel_raw(ind,:),squeeze(mean(mean(US_activation_response{ii,1}(:,:,ind),1,'omitnan'),2,'omitnan')));view(v(1),v(2));
                        set(gca,'fontsize',16);title(['Acq.' num2str(ii)]);
                    end
                    saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['UR AMP of US responsive neuron']),'jpeg');
                    %bar of CR AND UR
                    %             x=[NUM_UP_UR];
                    %             h=figure;b=bar([NUM_UP_CR(:,iscutmov)';x;NUM_DOWN_CR(:,iscutmov)']'*100./size(A,2));
                    %             xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;labels = num2str(b(1).YData','%0.2f');
                    %             xtips3 = b(3).XEndPoints;ytips3 = b(3).YEndPoints;labels3 = num2str(b(3).YData','%0.2f');
                    %             xtips2 = b(2).XEndPoints;ytips2 = b(2).YEndPoints;labels2 = num2str(b(2).YData','%0.2f');
                    %             text(xtips1,ytips1,labels,'HorizontalAlignment','center','VerticalAlignment','bottom');
                    %             text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
                    %             text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom')
                    %             ylabel('Neurons(%)');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);legend({'CS activation';'US';'CS inhibition'});box off;
                    %             saveas(h,fullfile(p,['iscutmov',num2str(iscutmov)],['Fraction of CS and US responsive neuron']),'jpeg');
                end
            end
        end
    end
    %     %% SAVE
    %     if  ~isfile(fullfile(checkpath(p),'CR_ind_summary_singletrial.mat'))
    %         a=[];save(fullfile(p,'CR_ind_summary_singletrial.mat'),'a');
    %     else
    % %         save(fullfile(p,'CR_ind_summary_singletrial.mat'),'A_r','trial_hab','trial_test',...
    % %             'cut_criterion','cs_win','us_win','base_win','pre_cs_time','CRon_activation_ind','CRon_inhibition_ind','CSon_activation_response',...
    % %             'CRoff_activation_ind','CRoff_inhibition_ind','CSoff_activation_response',...
    % %             'UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail',...
    % %             'onset_trial',...
    % %             'CS_all_hab','CS_all_tst',...
    % %             'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS_inthissession',...
    % %             'CSon_up_ind','CSon_down_ind','CSon_stable_ind','CSon_up_ind_all','CSon_down_ind_all','CSon_stable_ind_all',...
    % %             'CSoff_up_ind','CSoff_down_ind','CSoff_stable_ind','CSoff_up_ind_all','CSoff_down_ind_all','CSoff_stable_ind_all',...
    % %             'US_up_ind','US_down_ind','US_stable_ind','US_up_ind_all','US_down_ind_all','US_stable_ind_all',...
    % %             'US_up_ind_cuttrail','US_down_ind_cuttrail','US_stable_ind_cuttrail','US_up_ind_all_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
    % %             'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
    % %             'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS','-v7.3');
    % %          clear('A_r','trial_hab','trial_test',...
    % %             'cut_criterion','cs_win','us_win','base_win','pre_cs_time','CRon_activation_ind','CRon_inhibition_ind','CSon_activation_response',...
    % %             'CRoff_activation_ind','CRoff_inhibition_ind','CSoff_activation_response',...
    % %             'UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail',...
    % %             'onset_trial',...
    % %             'CS_all_hab','CS_all_tst',...
    % %             'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS_inthissession',...
    % %             'CSon_up_ind','CSon_down_ind','CSon_stable_ind','CSon_up_ind_all','CSon_down_ind_all','CSon_stable_ind_all',...
    % %             'CSoff_up_ind','CSoff_down_ind','CSoff_stable_ind','CSoff_up_ind_all','CSoff_down_ind_all','CSoff_stable_ind_all',...
    % %             'US_up_ind','US_down_ind','US_stable_ind','US_up_ind_all','US_down_ind_all','US_stable_ind_all',...
    % %             'US_up_ind_cuttrail','US_down_ind_cuttrail','US_stable_ind_cuttrail','US_up_ind_all_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
    % %             'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
    % %             'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS');
    % %
    %         save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region',...
    %             'Fraction_in_region_type_emerged_inthissession','Index_in_region_emerged_inthissession','Loc_in_region_cell_emerged_inthissession','Num_in_region_emerged_inthissession','Brain_region_id_emerged_inthissession','Label_region_emerged_inthissession',...
    %             'Fraction_in_region_type_emerged','Index_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged',...
    %             '-append','-v7.3');
    %     end
    
    
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
% 'CSon-avtivation': 对CS onset有显著大的响应的neuron（CS时间窗的平均幅值 > m +3sd）
% 'CSon-inhibition':对CS onset有显著小的响应的neuron（CS时间窗的平均幅值 < m +2sd）
% 'CSoff-avtivation',对CS offset有显著大的响应的neuron（CS offset后2s时间窗的平均幅值 > m +3sd）
% 'CSoff-inhibition',对CS offset有显著小的响应的neuron（CS offset后2s时间窗的平均幅值 < m +2sd）
% 'US-activation',对US onset（CS offset）有显著大的响应的neuron（USonset后2s时间窗的平均幅值 > m +3sd）
% 'US-activation-cuttrail',对US onset（CS offset）有显著大的响应的neuron【减去尾迹】（USonset后2s时间窗的平均幅值-CS时间窗的平均幅值 > m +3sd）
% 'CSon-up regulate',对CS onset有响应在当前session相比于pre cond显著增大的neuron，且在当前session或pre cond中任一trial被判断为CSon-avtivation（CS期间幅值与pre cond期间平均幅值做left tail ranksum test，p < 0.05）
% 'CSon-down regulate',对CS onset有响应在当前session相比于pre cond显著减小的neuron，且在当前session或pre cond中任一trial被判断为CSon-avtivation（CS期间幅值与pre cond期间平均幅值做right tail ranksum test，p < 0.05）
% 'CSon-stable regulate':既不属于'CSon-up regulate',又不属于'CSon-down regulate'，且在当前session或pre cond中任一trial被判断为CSon-avtivation
% 'CSoff-up regulate',对CS offset有响应在当前session相比于pre cond显著增大的neuron，且在当前session或pre cond中任一trial被判断为CSoff-avtivation（CS期间幅值与pre cond期间平均幅值做left tail ranksum test，p < 0.05）
% 'CSoff-down regulate',对CS offset有响应在当前session相比于pre cond显著减小的neuron，且在当前session或pre cond中任一trial被判断为CSoff-avtivation（CS期间幅值与pre cond期间平均幅值做right tail ranksum test，p < 0.05）
% 'CSoff-stable regulate'：既不属于'CSoff-up regulate',又不属于'CSoff-down regulate'，且在当前session或pre cond中任一trial被判断为CSoff-avtivation
% 'CSon-up regulate-all','CSon-down regulate-all','CSon-stable regulate-all','CSoff-up regulate-all','CSoff-down regulate-all','CSoff-stable regulate-all'：同上，不同的是不需要在当前session或pre cond中任一trial被判断为CSon-avtivation或CSoff-avtivation
% 'US-up regulate':对US onset有响应在当前session相比于Acq.1显著增大的neuron，且在当前session或Acq.1中任一trial被判断为US-avtivation（USonset后2s时间窗幅值与Acq.1期间平均幅值做left tail ranksum test，p < 0.05）
% 'US-down regulate':对US onset有响应在当前session相比于Acq.1显著减小的neuron，且在当前session或Acq.1中任一trial被判断为US-avtivation（USonset后2s时间窗幅值与Acq.1期间平均幅值做right tail ranksum test，p < 0.05）
% 'US-stable regulate'：既不属于'US-up regulate'，又不属于'US-down regulate'，且在当前session或Acq.1中任一trial被判断为US-avtivation
% 'US-up regulate-cuttrail','US-down regulate-cuttrail','US-stable regulate-cuttrail'：同上，不同的是用US期间减去尾迹（USonset后2s时间窗的平均幅值-CS时间窗的平均幅值）的幅值做显著检验
% 'US-up regulate-all','US-down regulate-all','US-stable regulate-all'，'US-up regulate-cuttrail-all','US-down regulate-cuttrail-all','US-stable regulate-cuttrail-all',：同上，不同的是不需要在当前session或Acq.1中任一trial被判断为US-avtivation
% 'CSUS-shift'：CS响应在当前session相比于pre cond显著增大，US响应在当前session相比于Acq.1显著减小，且在当前session或pre cond中任一trial被判断为CSon-avtivation,以及在当前session或Acq.1中任一trial被判断为US-avtivation
% 'CSUS-shift-cuttrail'：同上，不同的是用US期间减去尾迹（USonset后2s时间窗的平均幅值-CS时间窗的平均幅值）的幅值做显著检验
% 'CSUS-shift-all'：CS响应在当前session相比于pre cond显著增大，US响应在当前session相比于Acq.1显著减小
% 'CSUS-shift-cuttrail-all'：同上，不同的是用US期间减去尾迹（USonset后2s时间窗的平均幅值-CS时间窗的平均幅值）的幅值做显著检验
end
