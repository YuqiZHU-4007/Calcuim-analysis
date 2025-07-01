function US_mapping_singletrial_onlyUS(Path)
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
        %load(fullfile(p,'activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
        load(fullfile(p,'env.mat'));
        p=checkpath(fullfile(p,'singletrial'));
        %A=activities_preCS_dfdf_aftcorrect;%A=zscore(A,0,'all');
        %         A_r=reshape(A,frame.per_cycle,trial.total,[]);
        
        res=[0.66,0.66,10];
        load('X:\calcium data 20230224\脑区分割\segmentation_file_0525_DSregion_mask.mat');
        temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');
        warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
        temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
        nn=[Path(end-14:end-7),Path(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        
        load(fullfile(p,'CR_ind_summary_singletrial_statistic.mat'),'p_US_up','p_US_up_cuttrail','p_CSon_up');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'A_r','trial_hab','trial_test','CS_all_hab','CS_all_tst','CRon_activation_ind');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'), 'UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'labels_typei','ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS','ind_all_emergeCSUS_inthissession');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS');
        %load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Fraction_in_region_type');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Index_in_region');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Loc_in_region_cell');
        load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Num_in_region');
        %load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Brain_region_id');
       % load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Label_region');
        p_US_down=nan(size(A_r,3),trial.total);
        p_US_down_cuttrail=nan(size(A_r,3),trial.total);
        q_US_down=p_US_down;
        q_US_down_cuttrail=p_US_down_cuttrail;
        fdr_US_down=p_US_down;
        fdr_US_down_cuttrail=p_US_down_cuttrail;
        ranksumtest=1;
         Brain_region_id={};
        US_down_ind=cell(trial.total,2); US_stable_ind=cell(trial.total,2);
        US_down_ind_all=cell(trial.total,2); US_stable_ind_all=cell(trial.total,2);
        US_down_ind_cuttrail=cell(trial.total,2); US_stable_ind_cuttrail=cell(trial.total,2);
        US_down_ind_all_cuttrail=cell(trial.total,2); US_stable_ind_all_cuttrail=cell(trial.total,2);
        CS_US_shift_ind=cell(trial.total,2); CS_US_shift_ind_cuttrail=cell(trial.total,2);
        CS_US_shift_ind_all=cell(trial.total,2); CS_US_shift_ind_all_cuttrail=cell(trial.total,2);
        
        for iscutmov=1
            %% UP/DOWN CS/US-responsive neuron
            %down/up-regulate
            if ranksumtest==1;
                x=US_activation_response;
                a=[];
                for ii=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A_r,3)
                    for jj=trial.acq(2):trial.acq(3)
                        p_US_down(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','right');%%%%%%%%%%%%%%%%
                    end
                end
                x=US_activation_response_cuttrail;
                a=[];
                for ii=trial.acq(2):trial.acq(2)+trial.acq_block_trial-1
                    a(:,ii,:)=squeeze(x{ii}(:,:,:));
                end
                base_CR=squeeze(mean(a,2,'omitnan'));
                for ii=1:size(A_r,3)
                    for jj=trial.acq(2):trial.acq(3)
                        p_US_down_cuttrail(ii,jj)=ranksum(base_CR(:,ii),squeeze(x{jj}(:,:,ii)),'tail','right');%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    end
                end
                for jj=1:trial.total
                    [fdr_US_down(:,jj,iscutmov),q_US_down(:,jj,iscutmov)] = FDR_base(p_US_down(:,jj,iscutmov),'BH');
                    [fdr_US_down_cuttrail(:,jj,iscutmov),q_US_down_cuttrail(:,jj,iscutmov)] = FDR_base(p_US_down_cuttrail(:,jj,iscutmov),'BH');
                end
                save(fullfile(checkpath(p),'CR_ind_summary_singletrial_statistic.mat'),'p_US_down','q_US_down','fdr_US_down','p_US_down_cuttrail','q_US_down_cuttrail','fdr_US_down_cuttrail','-append');
            end
            for ii=1:trial.total
                US_down_ind{ii,iscutmov}=find(squeeze(p_US_down(:,ii))<=0.05 & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1);
                US_stable_ind{ii,iscutmov}=find((squeeze(p_US_up(:,ii))>0.05 & squeeze(p_US_down(:,ii))>0.05) & squeeze(sum(UR_activation_ind(:,[ii,trial.acq(2)]),2))>=1);
                US_down_ind_all{ii,iscutmov}=find(squeeze(p_US_down(:,ii))<=0.05);
                US_stable_ind_all{ii,iscutmov}=find((squeeze(p_US_up(:,ii))>0.05 & squeeze(p_US_down(:,ii))>0.05));
                
                US_down_ind_cuttrail{ii,iscutmov}=find(squeeze(p_US_down_cuttrail(:,ii))<=0.05 & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1);
                US_stable_ind_cuttrail{ii,iscutmov}=find((squeeze(p_US_up_cuttrail(:,ii))>0.05 & squeeze(p_US_down_cuttrail(:,ii))>0.05) & squeeze(sum(UR_activation_ind_cuttrail(:,[ii,trial.acq(2)]),2))>=1);
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
            clear('p_US_up','p_US_up_cuttrail','p_CSon_up','CRon_activation_ind''p_US_down','q_US_down','fdr_US_down','p_US_down_cuttrail','q_US_down_cuttrail','fdr_US_down_cuttrail');
            
            %% 脑区统计
            for ii=1:trial.total
                ind_all_CSUS_RESPONSIVE{14,ii,iscutmov}=US_down_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{15,ii,iscutmov}=US_stable_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{17,ii,iscutmov}=US_down_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{18,ii,iscutmov}=US_stable_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{26,ii,iscutmov}=US_down_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{27,ii,iscutmov}=US_stable_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{29,ii,iscutmov}=US_down_ind_all_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{30,ii,iscutmov}=US_stable_ind_all_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{31,ii,iscutmov}=CS_US_shift_ind{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{32,ii,iscutmov}=CS_US_shift_ind_cuttrail{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{33,ii,iscutmov}=CS_US_shift_ind_all{ii,iscutmov};
                ind_all_CSUS_RESPONSIVE{34,ii,iscutmov}=CS_US_shift_ind_cuttrail{ii,iscutmov};
            end
            for jj=[14,15,17,18,26,27,29:34]
                for ii=1:1:trial.total
                    for zz=1:trial.total
                        a=ind_all_CSUS_RESPONSIVE{jj,ii,iscutmov};
                        b=ind_all_CSUS_RESPONSIVE{jj,zz,iscutmov};
                        ind_all_emergeCSUS{jj,ii,zz,iscutmov,1}=setdiff(a,b);
                    end
                end
            end
            for jj=[14,15,17,18,26,27,29:34]
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
            for jj=[14,15,17,18,26,27,29:34]
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
            try
            save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS','ind_all_emergeCSUS_inthissession',...
                'US_down_ind','US_stable_ind','US_down_ind_all','US_stable_ind_all',...
                'US_down_ind_cuttrail','US_stable_ind_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
                'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
                'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS','-append','-v7.3');
            catch
                warning([fullfile(p,'CR_ind_summary_singletrial.mat'), 'ind_all_CSUS_RESPONSIVE','ind_all_emergeCSUS','ind_all_emergeCSUS_inthissession',...
                'US_down_ind','US_stable_ind','US_down_ind_all','US_stable_ind_all',...
                'US_down_ind_cuttrail','US_stable_ind_cuttrail','US_down_ind_all_cuttrail','US_stable_ind_all_cuttrail',...
                'CS_US_shift_ind','CS_US_shift_ind_all','CS_US_shift_ind_cuttrail','CS_US_shift_ind_all_cuttrail',...
                'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS']);
            end
            clear('A_r','trial_hab','trial_test','CS_all_hab','CS_all_tst',...
                'cut_criterion','us_win','UR_activation_ind','US_activation_response','UR_activation_ind_cuttrail','US_activation_response_cuttrail',...
                'Response_all_CSUS_RESPONSIVE','Response_all_emergeCSUS_inthissession','Response_all_emergeCSUS',...
                'US_up_ind','US_down_ind','US_stable_ind','US_up_ind_all','US_up_ind_cuttrail','US_up_ind_all_cuttrail');
            %fraction
            
            %             load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Fraction_in_region_type_emerged_inthissession','Index_in_region_emerged_inthissession','Loc_in_region_cell_emerged_inthissession','Num_in_region_emerged_inthissession','Brain_region_id_emerged_inthissession','Label_region_emerged_inthissession');
            %             load(fullfile(p,'CR_ind_summary_singletrial.mat'),'Fraction_in_region_type_emerged','Index_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged');
            %
            fraction_type=1;
            for jj=[14,15,17,18,26,27,29:34]
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
                        %Label_region{jj,iscutmov,fraction_type}(regioni,:)=Label{regioni};
                    end
                end
            end
            if false
                for jj=[14,15,17,18,26,27,29:34]
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
                for jj=[14,15,17,18,26,27,29:34]
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
            end
            try 
            save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','-append','-v7.3');%'Brain_region_id','Label_region',
            %          save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'Fraction_in_region_type_emerged_inthissession','Index_in_region_emerged_inthissession','Loc_in_region_cell_emerged_inthissession','Num_in_region_emerged_inthissession','Brain_region_id_emerged_inthissession','Label_region_emerged_inthissession','-append','-v7.3');
            %            save(fullfile(p,'CR_ind_summary_singletrial.mat'), 'Fraction_in_region_type_emerged','Index_in_region_emerged','Loc_in_region_cell_emerged','Num_in_region_emerged','Brain_region_id_emerged','Label_region_emerged','-append','-v7.3');
            catch
                warning([fullfile(p,'CR_ind_summary_singletrial.mat'), 'Fraction_in_region_type','Index_in_region','Loc_in_region_cell','Num_in_region','Brain_region_id','Label_region']);
            end
        end
    end
    
end
