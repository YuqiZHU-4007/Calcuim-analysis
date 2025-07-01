%% avg
clc;clear all;
Region_name='L H';isslidewin=1;
savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
Region_AUC_session={};Region_fraction_session={}; learning_perform={};
temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');res=[0.66,0.66,10];
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
savepath=checkpath('X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\Hb_VS_behavior\');
num{1}=[1:3];num{2}=[1:3];num{3}=[1:3];savetype='jpg';
for batchi=1:3
    for fishi=num{batchi}
        path=Path{batchi}{fishi};%'H:\1.Test US\5.fear conditioning behavioral data\20220709\fish1\';
        load(fullfile(path, '/activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(path, '/brain_region_related_statistic.mat'),'Label','brain_region_id');
        load(fullfile(path,'/para.mat'));
        load(fullfile(path,'/behavior_from_Results_of_alltheta.mat'));
        %load(fullfile(path,'singletrial','CR_ind_summary_singletrial.mat'),'ind_all_CSUS_RESPONSIVE');
        load(fullfile(path, '/env.mat'));
        

        activities_preCS_dfdf_aftcorrect=activities_preCS_dfdf_aftcorrect(Time.T_non_spon_ca,:);
        activities_preCS_dfdf_aftcorrect=zscore(activities_preCS_dfdf_aftcorrect,0,'all');
       
        nn=[path(end-13:end-6),path(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        Region_id=find(strcmp(brain_region_id(:,2),Region_name));
        ind_cell_opt_layer=find(env.supervoxel(:,3)>=opt.start_step+1 & env.supervoxel(:,3)<=opt.start_step+6-1);
       Region_id=intersect(Region_id,ind_cell_opt_layer);
        
        figure,subplot(1,3,1);scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),5,'filled');hold on;
        scatter3(env.supervoxel(Region_id,1),env.supervoxel(Region_id,2),env.supervoxel(Region_id,3),10,'filled');axis equal
        subplot(1,3,2);scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,'filled');hold on;
        scatter3(supervolxeli(Region_id,1),supervolxeli(Region_id,2),supervolxeli(Region_id,3),10,'filled');axis equal
        subplot(1,3,3);scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),5,'filled');hold on;
        scatter3(supervolxeli(Region_id,1),supervolxeli(Region_id,2),supervolxeli(Region_id,3),10,'filled');axis equal
% 
%         figure,
%         scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,'filled');hold on;
%         scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),5,'filled');hold on;
        CS_win=frame.cs_start:frame.cs_end-1;
        %% Learning performance
        [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint_sd,[],[],[],1);
        learning_perform{batchi}(:,fishi)=cat(1,bef_cond.CS_shake,dur_cond.CS192_shake,aft_cond.CS_shake);
        %% mean AUCre_startpoint_sd
        a=mean(activities_preCS_dfdf_aftcorrect(:,Region_id),2,'omitnan');
        Region_act=reshape(a,frame.per_cycle,trial.total,[]);
        h=figure;plot_rawtrace_trials(a,[],fs,frame,trial,[],1);
        saveas(h,string(fullfile(savepath,[nn,'_Meantrace'])),savetype);
        
        for sessioni=1:size(Region_act,2)
            area = cumtrapz(CS_win,Region_act(CS_win,sessioni));
            Region_AUC_session{batchi}(sessioni,fishi)= area(end) - area(1);
        end
%         %% CS_activation fraction
%         for sessioni=1:size(Region_act,2)
%             a=ind_all_CSUS_RESPONSIVE{1,sessioni,1};
%             b=union(a,Region_id);
%             Region_fraction_session{batchi}(sessioni,fishi)=length(b)/length(Region_id);
%         end
    end
end
save([savepath '\summary.mat'],'Region_AUC_session','Region_fraction_session','learning_perform');
%% plot

clr_group=[colorplus(389);colorplus(233);colorplus(253)];
group_label={'Learner(n=3)','Non-Learner(n=3)','Faded-Learner(n=3)'};
trial_ind_session=[];
if isslidewin==1
    slidwin=6;
    for sessioni=1:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
        switch sessioni
            case 1
                trial_ind_session(:,sessioni)=trial.hab(2):trial.hab(3);
            case mat2cell([2:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+1],1,ones(1,(trial.acq_block_num*trial.acq_block_trial-slidwin)+1));
                trial_ind_session(:,sessioni)=trial.hab(3)+sessioni-1 :trial.hab(3)+sessioni+slidwin-2;
            case (trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
                trial_ind_session(:,sessioni)=trial.test(2):trial.test(3);
        end
    end
else
    for sessioni=1:trial.acq_block_num+2
        switch sessioni
            case 1
                trial_ind_session(:,sessioni)=trial.hab(2):trial.hab(3);
            case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num+1));
                trial_ind_session(:,sessioni)=trial.hab(3)+(sessioni-1)*trial.acq_block_trial :trial.hab(3)+(sessioni)*trial.acq_block_trial;
            case trial.acq_block_num+2
                trial_ind_session(:,sessioni)=trial.test(2):trial.test(3);
        end
    end
end
xlab{1}='Pre Cond.';xlab{size(trial_ind_session,2)}='Post Cond.';
for sessioni=2:size(trial_ind_session,2)-1
    xlab{sessioni}=strcat('Cond.',num2str(sessioni-1));
end
%single fish（时间轴是session)
for typei=1
    for batchi=1:length(Path)
        switch typei
            case 1
                name=strcat(group_label(batchi),'_AUC of HbmeanCSact');a=Region_AUC_session{batchi};ylab='AUC';
            case 2
                name=strcat(group_label(batchi),'_Fraction of HbCSactivation');a=Region_fraction_session{batchi};ylab='Fraction';
        end
        h=figure('position',[33,30,1900,ceil(length(num{batchi})/4)*250],'color','w');t=tiledlayout(ceil(length(num{batchi})/4),4);t.TileSpacing = 'compact';t.Padding = 'compact';
        for fishi=num{batchi}
            path=Path{batchi}{fishi};nn=[path(end-13:end-6),path(end-5:end-1)];
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);
                y1(sessioni)=mean(a(ind,fishi),1,'omitnan');sd1(sessioni)=std(a(ind,fishi),[],1,'omitnan');
                y2(sessioni)=mean(learning_perform{batchi}(ind,fishi),1,'omitnan');
            end
            x=1:size(trial_ind_session,2);
            yyaxis left;shadedErrorBar(x,y1,sd1,'lineprops',{[1 0 0]},'transparent',1,'patchSaturation',0.1); hold on
            ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
            yyaxis right;
            plot(x, y2,'k-.','linewidth',1.5,'markersize',8);hold on
            ylabel('CR ratio','fontsize',16,'FontWeight','bold');ylim([0 1.2]);set(gca,'ycolor',[0 0 0],'linewidth',2);
            xlim([1 max(x)]);
            set(gca,'xtick',x,'xticklabel',xlab,'XTickLabelRotation',45,'fontsize',12);
            title(nn,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        end
        saveas(h,string(fullfile(savepath,name)),savetype);
    end
end
%single fish（时间轴是session,按第一个session，normal)
for typei=1
    for batchi=1:length(Path)
        switch typei
            case 1
                name=strcat(group_label(batchi),'_normAUC of HbmeanCSact');a=Region_AUC_session{batchi};ylab='norm. AUC';
            case 2
                name=strcat(group_label(batchi),'_normFraction of HbCSactivation');a=Region_fraction_session{batchi};ylab='norm.Fraction';
        end
        h=figure('position',[33,30,1900,ceil(length(num{batchi})/4)*250],'color','w');t=tiledlayout(ceil(length(num{batchi})/4),4);t.TileSpacing = 'compact';t.Padding = 'compact';
        for fishi=num{batchi}
            path=Path{batchi}{fishi};nn=[path(end-13:end-6),path(end-5:end-1)];
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);
                y1(sessioni)=mean(a(ind,fishi),1,'omitnan');
                y2(sessioni)=mean(learning_perform{batchi}(ind,fishi),1,'omitnan');
            end
            x=1:size(trial_ind_session,2);y1=(y1-y1(1))./(y1(1));%sd1=(sd1-sd1(1))./(sd1(1));
            yyaxis left;plot(x, y1,'r','linewidth',1.5,'markersize',8);hold on
            ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
            ylim([-1 4]);
            yyaxis right;
            plot(x, y2,'k-.','linewidth',1.5,'markersize',8);hold on
            ylabel('CR ratio','fontsize',16,'FontWeight','bold');ylim([0 1.2]);set(gca,'ycolor',[0 0 0],'linewidth',2);
            xlim([1 max(x)]);
            set(gca,'xtick',x,'xticklabel',xlab,'XTickLabelRotation',45,'fontsize',12);
            title(nn,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        end
        saveas(h,string(fullfile(savepath,name)),savetype);
    end
end
% Hbact VS performance
for typei=1
    h=figure('position',[265,126,973,769]);
    for batchi=[1:length(Path)]
        x=[];y=[];
        switch typei
            case 1
                name='AUC VS perform';a=Region_AUC_session{batchi};ylab='AUC';lim=[0 25];
            case 2
                name='Fraction VS perform';a=Region_fraction_session{batchi};ylab='Fraction';lim=[0 1];
        end
        for sessioni=1:size(trial_ind_session,2)
            ind=trial_ind_session(:,sessioni);
            x(sessioni,:)=mean(a(ind,num{batchi}),1,'omitnan');
            y(sessioni,:)=mean(learning_perform{batchi}(ind,num{batchi}),1,'omitnan');
        end
        subplot(2,2,batchi),
        scatter(reshape(x,1,[]),reshape(y,1,[]),15,clr_group(batchi,:),'filled');hold on;
        ylabel(['CR ratio','(slidewin=6)'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        title(group_label(batchi));
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
% Hbact VS performance(normal)
for typei=1
    h=figure('position',[265,126,973,769]);
    for batchi=[1:length(Path)]
        x=[];y=[];
        switch typei
            case 1
                name='normAUC VS perform';a=Region_AUC_session{batchi};ylab='norm. AUC';lim=[-1 2];
            case 2
                name='normFraction VS perform';a=Region_fraction_session{batchi};ylab='norm. Fraction';lim=[0 1];
        end
        for sessioni=1:size(trial_ind_session,2)
            ind=trial_ind_session(:,sessioni);
            x(sessioni,:)=mean(a(ind,num{batchi}),1,'omitnan');
            y(sessioni,:)=mean(learning_perform{batchi}(ind,num{batchi}),1,'omitnan');
        end
        subplot(2,2,batchi),x=(x-repmat(x(1,:),size(trial_ind_session,2),1))./(repmat(x(1,:),size(trial_ind_session,2),1));
        scatter(reshape(x,1,[]),reshape(y,1,[]),15,clr_group(batchi,:),'filled');hold on;
        ylabel(['CR ratio','(slidewin=6)'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        title(group_label(batchi));
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%Hb mean AUC in 前8个session VS learning speed
for typei=1
    h=figure('position',[265,126,900,300]);
    for indi=1:3
        ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
        learning_speed={};    
        subplot(1,3,indi)
        for batchi=[1,length(Path)]
            for fishi=num{batchi}
                path=Path{batchi}{fishi};
                load(fullfile(path, '/activities_dfdf_align.mat'),'align_win');
                learning_speed{batchi}(:,fishi)=max(align_win(:,7))-trial.hab(1);
            end
            switch typei
                case 1
                    name=['AUC VS learningspeed_session'];a=Region_AUC_session{batchi};ylab='AUC';lim=[0 25];
                case 2
                    name=['Fraction VS learningspeed_session'];a=Region_fraction_session{batchi};ylab='Fraction';lim=[0 1];
            end
            y= learning_speed{batchi}(num{batchi});
            x= mean(a(ind,num{batchi}),1,'omitnan');
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 24])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end


