%% refer to Schnitzer_2017_nature,zyq,20230615
%% down-sample:1 s bin
%% construct population vector: averaging the evoked neural responses over all five presentations of each stimulus and all the time bins associated with each stimulus presentation
%% CS/US-responsive ID :restricted the changes to the CS+ population vector to those cells that were US-responsive
%% Mahalanobis population vector distance (PVD)
%% PVD：1、 CS and first US 2、CS and aftercond-CS 3、CS and US
clc;clear all;
addpath(genpath('F:\DUlab\FC_analyse\Nomoto2022_NatureCommunications-main'));
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load(fullfile('H:\1.Test US\2.Tail free——Data from 117\20220803\fish2','brain_region_related_statistic.mat'));
sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
load([savepath '\Path']);
iscutmov=1;
frame_CS=1:12;
frame_US=1:3;
bin=1;%binmum=length(frame_CS)/bin;
for batchi=1:4
    path=Path{batchi};
    for fishi=1:length(path)
        %p='H:\1.Test US\2.Tail free——Data from 117\20210709\fish2';
        p=path{fishi};
        %load(fullfile(p,'CR_ind_summary.mat'));
        load(fullfile(p,'CR_ind_summary.mat'), 'CR_ind_up','CR_ind_down','CR_ind_up_acq','CR_ind_down_acq');
        load(fullfile(p,'CR_ind_summary.mat'),'ind_all_CSUS_RESPONSIVE','CS_response_acq','US_response','CS_response');
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'));
        type={'All','CS-active responsive','CS-down responsive'};%{'All','US-responsive','CS-responsive','CS&US-responsive'};
        ind_ensemble_restrict={};
        Mahalanobis_d_raw_CSUS={};Mahalanobis_mean_CSUS=[];Mahalanobis_d_raw_region_CSUS={};
        Mahalanobis_d_raw_CSUS1={};Mahalanobis_mean_CSUS1=[];Mahalanobis_d_raw_region_CSUS1={};Mahalanobis_mean_region_CSUS1=[];
        Mahalanobis_d_raw_CSCS={};Mahalanobis_mean_CSCS=[];Mahalanobis_d_raw_region_CSCS={};Mahalanobis_mean_region_CSCS=[];
        Mahalanobis_mean_region_CSUS=nan(length(type),trial.acq_block_num+2,4,length(Label));
        CR_ds_PV={};UR_ds_PV={};
        %% CS/US-responsive ID
        ind_ensemble_restrict{1}=1:size(CR_ind_up,1);
        CS_up_all=find(squeeze(sum(CR_ind_up(:,:,iscutmov),2))>=1 | squeeze(sum(CR_ind_up_acq(:,:,iscutmov),2))>=1);
        CS_down_all=find(squeeze(sum(CR_ind_down(:,:,iscutmov),2))>=1 | squeeze(sum(CR_ind_down_acq(:,:,iscutmov),2))>=1);
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
        for typei=1:length(type)
            ind=ind_ensemble_restrict{typei};
            for sessioni=1:trial.acq_block_num+2
                if sessioni>=2 & sessioni<trial.acq_block_num+2
                    CR=CS_response_acq{sessioni-1,iscutmov}(frame_CS,:,:);
                    UR=US_response{sessioni-1}(frame_US,:,:);
                elseif sessioni==trial.acq_block_num+2
                    CR=CS_response{sessioni-(trial.acq_block_num),iscutmov}(frame_CS,:,:);
                    UR=US_response{trial.acq_block_num}(frame_US,:,:);
                elseif sessioni==1
                    CR=CS_response{1,iscutmov}(frame_CS,:,:);
                    UR=US_response{1}(frame_US,:,:);
                end
                %% down-sample:1 s bin
                x=CR;
                a=reshape(x,bin,size(x,1)/bin,size(x,2),size(x,3));
                CR_ds=mean(a,1,'omitnan');
                x=UR;
                a=reshape(x,bin,size(x,1)/bin,size(x,2),size(x,3));
                UR_ds=mean(a,1,'omitnan');
                %% construct population vector
                CR_ds_PV{typei,sessioni}=squeeze(mean(CR_ds,3,'omitnan'));
                UR_ds_PV{typei,sessioni}=squeeze(mean(UR_ds,3,'omitnan'));
            end
            for sessioni=1:trial.acq_block_num+2
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
                         [Mahalanobis_d_raw_region_CSUS1{typei,sessioni,str2num(regioni)}, Mahalanobis_mean_region_CSUS1(typei,sessioni,:,str2num(regioni))]=my_mahal(UR_ds_PV{1}(:,ind_regioni_typei), CR_ds_PV{sessioni}(:,ind_regioni_typei));
                          [Mahalanobis_d_raw_region_CSCS{typei,sessioni,str2num(regioni)}, Mahalanobis_mean_region_CSCS(typei,sessioni,:,str2num(regioni))]=my_mahal(CR_ds_PV{end}(:,ind_regioni_typei), CR_ds_PV{sessioni}(:,ind_regioni_typei));
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
        save(fullfile(p,'\MPVD_update_20230615.mat'),'iscutmov','frame_CS','frame_US','bin',...
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
mat_name=[];
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);
savepath=checkpath('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\MPVD_CSCS\');
clr=[1 0 0;0.5 0.5 0.5;0 0 1;0 1 0];
MPVD_pca_all_region={};
for batchi=1:4  
   path=Path{batchi};
   %a=[];b=[];c=[];aa=[];bb=[];cc=[];aaa=[];bbb=[];ccc=[];
   for fishi=1:num{batchi}
       p=path{fishi};
       load(fullfile(p,'\MPVD_update_20230615.mat'));
       filename=Mahalanobis_mean_region_CSUS;
       for typei=1:length(type)
           for cal_type=1:size(filename,3)
               MPVD_pca_all_region{batchi,typei,cal_type}(fishi,:,:) = squeeze(filename(typei,:,cal_type,:));
           end
       end
       %         aaa(:,fishi)=CS_activation/ALL;bbb(:,fishi)=CS_inhibition/ALL;ccc(:,fishi)=US_activation/ALL;
       %         a(:,fishi)=CS_up/ALL;b(:,fishi)=CS_down/ALL;c(:,fishi)=CS_stable/ALL;
       %         aa(:,fishi)=US_up/ALL;bb(:,fishi)=US_down/ALL;cc(:,fishi)=US_stable/ALL;
   end
    %     aa=[(MPVD_all-repmat(MPVD_all(:,1),1,size(MPVD_all,2)))./repmat(MPVD_all(:,1),1,size(MPVD_all,2))*100]';
    %     figure,subplot(2,1,1),plot(aa,'k-*','linewidth',2);ylabel({'Relative change in' 'PVD to US(%)'}); %ylim([-20 10]);
    %     xticks([1:10]);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off
    %     subplot(2,1,2),plot([(MPVD_pca_all-repmat(MPVD_pca_all(:,1),1,size(MPVD_pca_all,2)))./repmat(MPVD_pca_all(:,1),1,size(MPVD_pca_all,2))*100]','k-*','linewidth',2);ylabel({'Relative change in' 'PVD to US(%)'}); %ylim([-20 10]);
    %     set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;
end
for session_cut=1
    for  typei=1:length(type)
        for cal_type=1:size(filename,3)
            h=figure('position',[1,1,1920,980]);
            for batchi=[2,3,4,1]
                for ii=1:size(Mahalanobis_mean_region_CSUS,4)
                    try
                        y=real(squeeze(MPVD_pca_all_region{batchi,typei,cal_type}(:,:,ii)));
                        yy=[(y-repmat(y(:,session_cut),1,size(y,2)))./(y+repmat(y(:,session_cut),1,size(y,2)))]';
                        x=1:10;
                        m=squeeze(mean(yy,2,'omitnan'));sd=squeeze(std(yy,[],2,'omitnan'));sem=sd./sqrt(size(yy,2));
                        subplot(5,6,ii),shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                        plot(x,m,'color',clr(batchi,:),'linewidth',2);
                        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);
                        ylabel({'Relative change in' 'PVD to US(%)'});
                        title(Label{ii});
                    catch
                        continue;
                    end
                end
            end
            saveas(h,[savepath,type{typei},'-cal_type',num2str(cal_type),'-session_cut',num2str(session_cut)],'jpg');
            close(h)
        end
    end
end
save([savepath,'MPVD.mat'],'MPVD_pca_all_region','-v7.3');
close all;

batchi=1;typei=1;cal_type=3;
clr=hsv(8);
h=figure('position',[1,1,1920,980]);colororder(clr)
for ii=1:30
y=real(squeeze(MPVD_pca_all_region{batchi,typei,cal_type}(:,:,ii)));
yy=y;% yy=[(y-repmat(y(:,session_cut),1,size(y,2)))./(y+repmat(y(:,session_cut),1,size(y,2)))]';
subplot(5,6,ii),plot(1:10,yy,'linewidth',2); title(Label{ii});
end
%% plot single fish
session_cut=1;typei=2;
h=figure('position',[1,1,1920,980]);
for ii=1:size(Mahalanobis_mean_region_CSUS,4)
    try
        a=abs(squeeze(Mahalanobis_mean_region_CSUS(typei,:,1,ii)));
        y=(squeeze(a));
        yy=[(y-repmat(y(:,session_cut),1,size(y,2)))./(repmat(y(:,session_cut),1,size(y,2)))]';
        x=1:10;
        m=squeeze(mean(yy,2,'omitnan'));sd=squeeze(std(yy,[],2,'omitnan'));sem=sd./sqrt(size(yy,2));
        subplot(5,6,ii),shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
        plot(x,m,'color',clr(batchi,:),'linewidth',2);
        set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);
        ylabel({'Relative change in' 'PVD to US(%)'});
        title(Label{ii});
    catch
        continue;
    end
end

%% plot all fish
% iscutmov=1;
% h=figure('position',[1,1,1920,980]);
% clr=[1 0 0;0.5 0.5 0.5;0 0 1;1 0 1];
% for batchi=1:2
%     path=Path{batchi};MPVD_pca_all=[];MPVD_all=[];MPVD_pca_all_region=[];MPVD_all_region=[];
%     a=[];b=[];c=[];aa=[];bb=[];cc=[];aaa=[];bbb=[];ccc=[];
%     for fishi=1:length(path)
%         p=path{fishi};
%         load(fullfile(p,'\MPVD_update_20230615.mat'));%load(fullfile(p,'\CS_US_responses_num.mat'));
%         MPVD_pca_all(fishi,:)= MPVD_pca;
%         MPVD_all(fishi,:)= MPVD;
%         MPVD_pca_all_region(fishi,:,:)= MPVD_pca_region;
%         MPVD_all_region(fishi,:,:)=MPVD_region;
% %         aaa(:,fishi)=CS_activation/ALL;bbb(:,fishi)=CS_inhibition/ALL;ccc(:,fishi)=US_activation/ALL;
% %         a(:,fishi)=CS_up/ALL;b(:,fishi)=CS_down/ALL;c(:,fishi)=CS_stable/ALL;
% %         aa(:,fishi)=US_up/ALL;bb(:,fishi)=US_down/ALL;cc(:,fishi)=US_stable/ALL;
%     end
%     %     aa=[(MPVD_all-repmat(MPVD_all(:,1),1,size(MPVD_all,2)))./repmat(MPVD_all(:,1),1,size(MPVD_all,2))*100]';
%     %     figure,subplot(2,1,1),plot(aa,'k-*','linewidth',2);ylabel({'Relative change in' 'PVD to US(%)'}); %ylim([-20 10]);
%     %     xticks([1:10]);set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off
%     %     subplot(2,1,2),plot([(MPVD_pca_all-repmat(MPVD_pca_all(:,1),1,size(MPVD_pca_all,2)))./repmat(MPVD_pca_all(:,1),1,size(MPVD_pca_all,2))*100]','k-*','linewidth',2);ylabel({'Relative change in' 'PVD to US(%)'}); %ylim([-20 10]);
%     %     set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off;
%     
%     %     for ii=1:size(MPVD_pca_region,1)
%     %         y=squeeze(MPVD_all_region(:,ii,:));
%     %         yy=[(y-repmat(y(:,2),1,size(y,2)))./repmat(y(:,2),1,size(y,2))*100]';
%     %         x=1:10;
%     %         m=squeeze(mean(yy,2,'omitnan'));sd=squeeze(std(yy,[],2,'omitnan'));sem=sd./sqrt(size(yy,2));
%     %         subplot(5,6,ii),shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.5);hold on
%     %         plot(x,m,'color',clr(batchi,:),'linewidth',2);
%     %         set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);
%     %         ylabel({'Relative change in' 'PVD to US(%)'});
%     %         title(Label{ii});
%     %     end
%     
%     %          h=figure('position',[1,1,1920,980]);
%     for sessioni=1:size(MPVD_pca_region,1)
%         y=squeeze(MPVD_pca_all_region(:,sessioni,:));
%         yy=[(y-repmat(y(:,1),1,size(y,2)))./repmat(y(:,1),1,size(y,2))*100]';
%         x=1:10;
%         m=squeeze(mean(yy,2,'omitnan'));sd=squeeze(std(yy,[],2,'omitnan'));sem=sd./sqrt(size(yy,2));
%         subplot(5,6,sessioni),shadedErrorBar(x,m,sem,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.5);hold on
%         plot(x,m,'color',clr(batchi,:),'linewidth',2);
%         set(gca,'xtick',[1:10],'xticklabel',sessionx,'XTickLabelRotation',90);xlim([1 10]);
%         ylabel({'Relative change in' 'PVD to US(%)'});
%         title(Label{sessioni});
%     end
% end

%% 统计每条鱼CS activation/inhibition US activation CS/US up/down/sable-regulate的比例
%                 CS_activation=[];CS_inhibition=[];US_activation=[];CS_up=[];CS_down=[];CS_stable=[];
%                 US_up=[];US_down=[];US_stable=[];
%                 for ii=1:trial.acq_block_num+2
%                     if ii>=2 & ii<trial.acq_block_num+2
%                         CS_activation(ii)=length(find(CR_ind_up_acq(:,ii-1,iscutmov)==1));
%                         CS_inhibition(ii)=length(find(CR_ind_down_acq(:,ii-1,iscutmov)==1));
%
%                     elseif ii==trial.acq_block_num+2
%                         CS_activation(ii)=length(find(CR_ind_up(:,2,iscutmov)==1));
%                         CS_inhibition(ii)=length(find(CR_ind_down(:,2,iscutmov)==1));
%                     elseif ii==1
%                         CS_activation(ii)=length(find(CR_ind_up(:,1,iscutmov)==1));
%                         CS_inhibition(ii)=length(find(CR_ind_down(:,1,iscutmov)==1));
%                     end
%
%                     %CS_US_overlap=
%                     CS_up(ii)=length(CS_up_ind{ii,iscutmov});
%                     CS_down(ii)=length(CS_down_ind{ii,iscutmov});
%                     CS_stable(ii)=length(CS_stable_ind{ii,iscutmov});
%                     US_up(ii)=length(US_up_ind{ii,iscutmov});
%                     US_down(ii)=length(US_down_ind{ii,iscutmov});
%                     US_stable(ii)=length(US_stable_ind{ii,iscutmov});
%                 end
%                 US_activation=[nan,NUM_UP_UR,nan];
%                 ALL=size(A_r,3);
%                 save(fullfile(p,'\CS_US_responses_num.mat'),'CS_activation','CS_inhibition','US_activation','CS_up','CS_down','CS_stable',...
%                     'US_up','US_down','US_stable','ALL')

% %CS
% x=CS_response{1};y=CS_response{2};
% p_CS_up_bef_aft=nan(size(A,2),1);p_CS_down_bef_aft=nan(size(A,2),1);
% for ii=1:size(A,2)
%     p_CS_up_bef_aft(ii)=ranksum(mean(squeeze(x(:,:,ii)),1,'omitnan'),mean(squeeze(y(:,:,ii)),1,'omitnan'),'tail','left');
%     p_CS_down_bef_aft(ii)=ranksum(mean(squeeze(x(:,:,ii)),1,'omitnan'),mean(squeeze(y(:,:,ii)),1,'omitnan'),'tail','right');
% end
% %[FDR_base,p_CS_responsive_bef_aft_adj] = mafdr( p_CS_responsive_bef_aft);
% length(find(p_CS_up_bef_aft<0.05))
% %US
% x=US_response{1};y=US_response{trial.acq_block_num};
% p_US_up_bef_aft=nan(size(A,2),1);p_US_down_bef_aft=nan(size(A,2),1);
% for ii=1:size(A,2)
%     p_US_up_bef_aft(ii)=ranksum(mean(squeeze(x(:,:,ii)),1,'omitnan'),mean(squeeze(y(:,:,ii)),1,'omitnan'),'tail','left');
%     p_US_down_bef_aft(ii)=ranksum(mean(squeeze(x(:,:,ii)),1,'omitnan'),mean(squeeze(y(:,:,ii)),1,'omitnan'),'tail','right');
% end
% %% plot pie
% CS_all=length(find(sum(CR_ind_up,2)>=1));
% CS_UP=length(find(p_CS_up_bef_aft<=0.05 & sum(CR_ind_up,2)>=1))/CS_all;
% CS_down=length(find(p_CS_down_bef_aft<=0.05 & sum(CR_ind_up,2)>=1))/CS_all;
% CS_stable=length(find((p_CS_up_bef_aft>0.05 & p_CS_down_bef_aft>0.05) & sum(CR_ind_up,2)>=1))/CS_all;
% explode = [1 1 0 ];labels = {'Up','Down','Stable'};
% figure('position',[740,615,876,262]),subplot(1,3,1),pie([CS_UP CS_down CS_stable],explode);title('CS');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2)
% US_all=length(find(sum(UR_ind,2)>=1));
% US_UP=length(find(p_US_up_bef_aft<=0.05 & sum(UR_ind,2)>=1))/CS_all;
% US_down=length(find(p_US_down_bef_aft<=0.05 & sum(UR_ind,2)>=1))/CS_all;
% US_stable=length(find((p_US_up_bef_aft>0.05 & p_US_down_bef_aft>0.05) & sum(UR_ind,2)>=1))/US_all;
% subplot(1,3,2),pie([US_UP US_down US_stable],explode);legend(labels,'position',[0.008,0.026,0.12,0.28]);title('US');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2)
% subplot(1,3,3),title('Responsive');b=bar(categorical({'CS','US','All'}),[CS_all US_all size(A,2)]);
% xtips1 = b(1).XEndPoints;ytips1 = b(1).YEndPoints;labels = string(b(1).YData);
% text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');set(gca,'fontsize',16,'FontWeight','bold','linewidth',2);box off
