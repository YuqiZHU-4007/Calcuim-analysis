%refer to Priya_2023
%correlation
%zyq,202304

%%
clc;close all;
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath,'/Path.mat']);
fc_region_all={};fc_region_all_mean={};
for batchi=1:4
    path=Path{batchi};
    for fishi=1:length(path)
        p=path{fishi};
        %p='H:\1.Test US\2.Tail free¡ª¡ªData from 117\20210709\fish2\'
        load(fullfile(p,'para.mat'));
        load(fullfile(p,'activities_aft_process.mat'));
        load(fullfile(p,'brain_region_related_statistic.mat'));
        %load(fullfile(p,'act_spon.mat'));
        %% Node
        %cell to cell whole brain
        %cell to cell within region
        %region to region
        mean_act_region=[];
        for regioni=1:30
            region_id=find(strcmp(brain_region_id,num2str(regioni)));
            mean_act_region(:,regioni)=mean(activities_preCS_dfdf_aftcorrect(:,region_id),2,'omitnan');
        end
        mean_act_region_r=reshape(mean_act_region,frame.per_cycle,trial.total,[]);
        frame_ind=frame.cs_start:frame.us_start-1 ;%1:frame.per_cycle;
        fc_region=[];
        for ii=1:10
            switch ii
                case 1
                    ind=trial.hab(2):trial.hab(3);
                case 10
                    ind=trial.test(2):trial.test(3);
                otherwise
                    ind=[trial.acq(2):trial.acq(2)+(trial.acq_block_trial-1)]+(ii-2)*trial.acq_block_trial;
            end
            for jj=1:length(ind)
                fc_region(:,:,jj,ii)=corr(squeeze(mean_act_region_r(frame_ind,ind(jj),:)));
            end
        end
        
        %figure('position',[0,0,1916,1200]),
        sessionx = {'Pre Cond','Cond 1','Cond 2','Cond 3','Cond 4','Cond 5','Cond 6','Cond 7','Cond 8','Post Cond'};
        c=[];
        for ii=1:10
            c(:,:,ii)=squeeze(mean(fc_region(:,:,:,ii),3,'omitnan'));
            %subplot(3,4,ii),
            h=figure;imagesc(abs(squeeze(c(:,:,ii))),[0,1]);colorbar;
            set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
            set(gca,'ytick',[1:30],'yticklabel',Label,'yTickLabelRotation',0);
            set(gca,'fontsize',10);title(sessionx{ii});
            saveas(h,[checkpath(fullfile(p,'fc')),'\',sessionx{ii}],'jpeg');
            close(h);
        end
        save([checkpath(fullfile(p,'fc')),'\fc_region.mat'],'mean_act_region','mean_act_region_r','fc_region','frame_ind','c','-v7.3');
        
        fc_region_all{batchi}(:,:,:,:,fishi)=fc_region;
        fc_region_all_mean{batchi}(:,:,:,fishi)=c;
    end
end

% figure,imagesc(abs(squeeze(c(:,:,10)-c(:,:,1))),[0,1]);colorbar;
% set(gca,'xtick',[1:30],'xticklabel',Label,'XTickLabelRotation',90)
% set(gca,'ytick',[1:30],'yticklabel',Label,'yTickLabelRotation',0);
% set(gca,'fontsize',10);title(sessionx{ii});

%% Edge
%% plot