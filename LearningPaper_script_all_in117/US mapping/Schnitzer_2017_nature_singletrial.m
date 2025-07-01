%% refer to Schnitzer_2017_nature,zyq,20230615
%% down-sample:1 s bin
%% construct population vector: averaging the evoked neural responses over all five presentations of each stimulus and all the time bins associated with each stimulus presentation
%% CS/US-responsive ID :restricted the changes to the CS+ population vector to those cells that were US-responsive
%% Mahalanobis population vector distance (PVD)
%% PVD��1、 CS and first US 2、CS and aftercond-CS 3、CS and US
clc;clear all;
savepath='X:\calcium data 20230224\';
load(fullfile('X:\calcium data 20230224\20220803\fish2','brain_region_related_statistic.mat'));
sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
load([savepath '\Path']);
iscutmov=1;
frame_CS=1:10;
frame_US=1:5;
bin=1;%binmum=length(frame_CS)/bin;
for batchi=1
    path=Path{batchi};
    for fishi=1
        %p='H:\1.Test US\2.Tail free！！Data from 117\20210709\fish2';
        p=path{fishi};
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'));
        pp=checkpath(fullfile(p,'singletrial'));
        load(fullfile(pp,'CR_ind_summary_singletrial.mat'), 'CRon_activation_ind','CRon_inhibition_ind');
        load(fullfile(pp,'CR_ind_summary_singletrial.mat'),'ind_all_CSUS_RESPONSIVE','CSon_activation_response','US_activation_response');
        type={'All','CS-active responsive','CS-down responsive'};
        ind_ensemble_restrict={};
        Mahalanobis_d_raw_CSUS={};Mahalanobis_mean_CSUS=[];Mahalanobis_d_raw_region_CSUS={};
        Mahalanobis_d_raw_CSUS1={};Mahalanobis_mean_CSUS1=[];Mahalanobis_d_raw_region_CSUS1={};Mahalanobis_mean_region_CSUS1=[];
        Mahalanobis_d_raw_CSCS={};Mahalanobis_mean_CSCS=[];Mahalanobis_d_raw_region_CSCS={};Mahalanobis_mean_region_CSCS=[];
        Mahalanobis_mean_region_CSUS=nan(length(type),trial.total,4,length(Label));
        CR_ds_PV={};UR_ds_PV={};
        %% CS/US-responsive ID
        ind_ensemble_restrict{1}=1:size(CRon_activation_ind,1);
        CS_up_all=find(squeeze(sum(CRon_activation_ind(:,:,iscutmov),2))>=1);
        CS_down_all=find(squeeze(sum(CRon_inhibition_ind(:,:,iscutmov),2))>=1);
        ind_ensemble_restrict{2}=CS_up_all;ind_ensemble_restrict{3}=CS_down_all;
        %         CS_up_all=find(squeeze(sum(CR_ind_up(:,:,iscutmov),2))>=1 | squeeze(sum(CR_ind_up_acq(:,:,iscutmov),2))>=1);
        %         CS_down_all=find(squeeze(sum(CR_ind_down(:,:,iscutmov),2))>=1 | squeeze(sum(CR_ind_down_acq(:,:,iscutmov),2))>=1);
        %         US_all=[];
        %         for jj=1:trial.acq_block_num+1
        %             US_all=union(ind_all_CSUS_RESPONSIVE{3,jj,iscutmov},US_all);
        %         end
        %         ind_ensemble_restrict{2}=US_all;
        %         ind_ensemble_restrict{3}=union(CS_up_all, CS_down_all);
        %         ind_ensemble_restrict{4}=union(ind_ensemble_restrict{2},ind_ensemble_restrict{3});
        %% PVD
        for typei=1
            ind=ind_ensemble_restrict{typei};
            for sessioni=1:trial.total
                CR=CSon_activation_response{sessioni,iscutmov}(frame_CS,:,:);
                x=CR;
                a=reshape(x,bin,size(x,1)/bin,size(x,2),size(x,3));%% down-sample:1 s bin
                CR_ds=mean(a,1,'omitnan');
                CR_ds_PV{typei,sessioni}=squeeze(mean(CR_ds,3,'omitnan'));  %% construct population vector
                if sessioni>=trial.acq(2) && sessioni<=trial.acq(3)
                    UR=US_activation_response{sessioni}(frame_US,:,:);x=UR;
                    a=reshape(x,bin,size(x,1)/bin,size(x,2),size(x,3));
                    UR_ds=mean(a,1,'omitnan');
                    UR_ds_PV{typei,sessioni}=squeeze(mean(UR_ds,3,'omitnan'));
                else
                    UR_ds_PV{typei,sessioni}=nan(length(frame_US)/bin,size(US_activation_response{trial.acq(2)},3));
                end
            end
            for sessioni=1:trial.total
                %% Mahalanobis population vector distance (PVD)
                reference_data=UR_ds_PV{typei,sessioni}(:,ind);%t * n
                target_data=CR_ds_PV{typei,sessioni}(:,ind);%t * n
                %                 [Mahalanobis_d_raw_CSUS{typei,sessioni}, Mahalanobis_mean_CSUS(typei,sessioni,:)]=my_mahal(reference_data, target_data); %3、CS and US
                %                 [Mahalanobis_d_raw_CSUS1{typei,sessioni}, Mahalanobis_mean_CSUS1(typei,sessioni,:)]=my_mahal(UR_ds_PV{1}(:,ind), CR_ds_PV{sessioni}(:,ind)); %1、 CS and first US
                %                 [Mahalanobis_d_raw_CSCS{typei,sessioni}, Mahalanobis_mean_CSCS(typei,sessioni,:)]=my_mahal(CR_ds_PV{end}(:,ind), CR_ds_PV{sessioni}(:,ind)); %2、 CS and aftercond-CS
                % PVD in per region
                for regioni=unique(brain_region_id(:,1))'
                    ind_regioni=find(strcmp(brain_region_id(:,1),regioni));
                    ind_regioni_typei=intersect(ind_regioni,ind);
                    if ~strcmp(regioni,"") & length(ind_regioni_typei)>=10
                        reference_data=UR_ds_PV{typei,sessioni}(:,ind_regioni_typei);%t * n
                        target_data=CR_ds_PV{typei,sessioni}(:,ind_regioni_typei);%t * n
                        [Mahalanobis_d_raw_region_CSUS{typei,sessioni,str2num(regioni)}, Mahalanobis_mean_region_CSUS(typei,sessioni,:,str2num(regioni))]=my_mahal(reference_data, target_data);
                        [Mahalanobis_d_raw_region_CSUS1{typei,sessioni,str2num(regioni)}, Mahalanobis_mean_region_CSUS1(typei,sessioni,:,str2num(regioni))]=my_mahal(UR_ds_PV{trial.acq(2)}(:,ind_regioni_typei), CR_ds_PV{sessioni}(:,ind_regioni_typei));
                        [Mahalanobis_d_raw_region_CSCS{typei,sessioni,str2num(regioni)}, Mahalanobis_mean_region_CSCS(typei,sessioni,:,str2num(regioni))]=my_mahal(CR_ds_PV{trial.test(2)}(:,ind_regioni_typei), CR_ds_PV{sessioni}(:,ind_regioni_typei));
                    end
                end
                %             y=UR;x=CR;
                %             x(find(isnan(x(:,1))),:)=[];
                %             [coeff,~,latent,~,explained] = pca(x);
                %             explained_var = cumsum(latent)/sum(latent);pcnum=min(find(explained_var>0.9));
                %             neuro_subspace = coeff(:,1:pcnum);
                %             CR_PAV=x*neuro_subspace;
                %             UR_PAV = y*neuro_subspace;
                %             MPVD_pca(ii)=mean(mahal(UR_PAV,CR_PAV));
                %             MPVD(ii)=mean(mahal(mean(CR,1)',mean(UR,1)'));
            end
        end
        save(fullfile(pp,'\MPVD_update_20230615.mat'),'iscutmov','frame_CS','frame_US','bin',...
            'CR_ds_PV','CR_ds_PV',...
            'type','ind_ensemble_restrict',...
            'Mahalanobis_d_raw_CSUS','Mahalanobis_mean_CSUS','Mahalanobis_d_raw_region_CSUS','Mahalanobis_mean_region_CSUS',...
            'Mahalanobis_d_raw_CSUS1','Mahalanobis_mean_CSUS1','Mahalanobis_d_raw_region_CSUS1','Mahalanobis_mean_region_CSUS1',...
            'Mahalanobis_d_raw_CSCS','Mahalanobis_mean_CSCS','Mahalanobis_d_raw_region_CSCS','Mahalanobis_mean_region_CSCS');
        %save(fullfile(p,'\MPVD.mat'),'MPVD','MPVD_pca','MPVD_region','MPVD_pca_region');
        %         figure,subplot(2,1,1),plot((MPVD-MPVD(1))./abs(MPVD(1))*100,'k-*','linewidth',2);ylabel({'Relative change in' 'PVD to US(%)'}); %ylim([-20 10]);
        %         xticks([1:trial.acq_block_num+2]);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off
        %         subplot(2,1,2),plot(MPVD,'k-*','linewidth',2);ylabel('PVD to US (A.U.)');xticks([1:trial.acq_block_num+2]);
        %         set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;title(p)
    end
end

%% plot together