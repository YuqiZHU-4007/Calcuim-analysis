%% avg
clc;clear all;
Region_name={'L H','R H'};isopto=0;savetype='jpg';
if isopto==1
    savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
    load([savepath '\Path']);
    warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
    savepath=checkpath('X:\calcium data 20230224\Hb_VS_behavior\Hb_act\');
else
    savepath='X:\calcium data 20230224\';
    load([savepath '\Path']);
    warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
    savepath=checkpath('X:\calcium data 20230224\Hb_VS_behavior\Norm\');
end
Region_AUC_session={};Region_fraction_session={};Region_AUC_session_all={}; Region_AUC_preCS_session_all={};Region_AUC_preCS_session={};
learning_speed={};learning_trial={};learning_perform={};Region_AUC_CS_session_all={};Region_AUC_CS_session={};
temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');res=[0.66,0.66,10];
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
for batchi=1:4
    for fishi=num{batchi}
        path=Path{batchi}{fishi};%'H:\1.Test US\5.fear conditioning behavioral data\20220709\fish1\';
        load(fullfile(path, '/activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
        load(fullfile(path, '/brain_region_related_statistic.mat'),'Label','brain_region_id');
        load(fullfile(path,'/behavior_from_Results_of_alltheta.mat'));
        %load(fullfile(path,'singletrial','CR_ind_summary_singletrial.mat'),'ind_all_CSUS_RESPONSIVE');
        load(fullfile(path, '/env.mat'));
        load(fullfile(path,'/para.mat'));
        load(fullfile(path, '/activities_dfdf_align.mat'),'align_win');
       
        %activities_preCS_dfdf_aftcorrect=zscore(activities_preCS_dfdf_aftcorrect,0,'all');
        nn=[path(end-14:end-7),path(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        Region_id=[];
        for ii=1:size(Region_name)
            a=find(strcmp(brain_region_id(:,2),Region_name{ii}));
            Region_id=cat(1,Region_id,a);
        end
        if isopto==1
            ind_cell_opt_layer=find(env.supervoxel(:,3)>=opt.start_step+1 & env.supervoxel(:,3)<=opt.start_step+6-1);
            Region_id=intersect(Region_id,ind_cell_opt_layer);
            Region_act=reshape(activities_preCS_dfdf_aftcorrect,frame.per_cycle,trial.total,[]);
            %             a=setdiff(1:size(activities_preCS_dfdf_aftcorrect,2),ind_cell_opt_layer);
            %             activities_preCS_dfdf_aftcorrect(find(activities_preCS_dfdf_aftcorrect==0))=nan;
            %             squeeze(mean(mean(Region_act(frame.cs_start:frame.cs_end-1,7:30,a),3,'omitnan'),2,'omitnan'));
        else
            Region_act=reshape(activities_preCS_dfdf_aftcorrect,frame.per_cycle,trial.total,[]);
        end

        trial_base=trial.hab(2):trial.hab(3);
        CS_win=frame.cs_start:frame.cs_end-1;
        pre_cs_time=[6 0];
        base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
        %% Learning performance
        [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint_sd,nn,[],savepath,1);
        learning_perform{batchi}(:,fishi)=cat(1,bef_cond.CS_shake,dur_cond.CS192_shake,aft_cond.CS_shake);
        learning_speed{batchi}(:,fishi)=(trial.acq(1)-(max(align_win(:,7))-trial.hab(1)))./trial.acq(1);
        learning_trial{batchi}(:,fishi)=max(align_win(:,7))-trial.hab(1);
        %% all AUC
        for sessioni=1:size(Region_act,2)
            a=Region_id;%ind_all_CSUS_RESPONSIVE{1,sessioni,1};
            b=intersect(a,Region_id);
            Region_AUC_session_all{batchi,fishi}(1:length(b),sessioni)=squeeze(mean(Region_act(CS_win,sessioni,b),1,'omitnan'));
            Region_AUC_preCS_session_all{batchi,fishi}(1:length(b),sessioni)=squeeze(mean(Region_act(base_win,sessioni,b),1,'omitnan'));
            Region_AUC_CS_session_all{batchi,fishi}(:,sessioni)=squeeze(mean(Region_act(CS_win,sessioni,:),1,'omitnan'));
        end
%         for sessioni=1:size(Region_act,2)
%             area = cumtrapz(squeeze(Region_act(CS_win,sessioni,:)));
%             Region_AUC_session_all{batchi,fishi}(:,sessioni)= area(end,:) - area(1,:);
%         end
        %% mean AUC
        a=mean(activities_preCS_dfdf_aftcorrect(:,Region_id),2,'omitnan');
        h=figure;plot_rawtrace_trials(a,[],fs,frame,trial,startpoint_sd,1);
        annotation(h,'textbox',[0.01 0.01 0.04 0.04],'String',nn,'Interpreter','none','fontsize',20,'EdgeColor','w'); 
        saveas(h,string(fullfile(savepath,[nn,'_Meantrace'])),savetype);
        close(h);
        %area = cumtrapz(squeeze(mean(Region_act(CS_win,:,:),3,'omitnan')));
        %Region_AUC_session{batchi}(:,fishi)= area(end,:) - area(1,:);
        %figure,histogram(activities_preCS_dfdf_aftcorrect)
        base=squeeze(mean(mean(mean(Region_act(CS_win,trial_base,:),1,'omitnan'),2,'omitnan'),3,'omitnan'));
        for sessioni=1:size(Region_act,2)
            a=Region_id;%ind_all_CSUS_RESPONSIVE{1,sessioni,1};
            b=intersect(a,Region_id);
            %Region_AUC_session{batchi}(sessioni,fishi)=squeeze(mean(mean(Region_act(CS_win,sessioni,b),1,'omitnan'),3,'omitnan'));
            Region_AUC_session{batchi}(sessioni,fishi)=squeeze(mean(mean(Region_act(CS_win,sessioni,b),1,'omitnan'),3,'omitnan'));
            Region_AUC_preCS_session{batchi}(sessioni,fishi)=squeeze(mean(mean(Region_act(base_win,sessioni,b),1,'omitnan'),3,'omitnan'));
            Region_AUC_CS_session{batchi}(sessioni,fishi)=squeeze(mean((mean(Region_act(CS_win,sessioni,b),1,'omitnan')-base)./base,3,'omitnan'));
           % Region_AUC_CS_session{batchi}(sessioni,fishi)=squeeze(sum((mean(Region_act(CS_win,sessioni,b),1,'omitnan')-base)./base,3,'omitnan'));
        end
        % CS_activation fraction
%         for sessioni=1:size(Region_act,2)
%             a=ind_all_CSUS_RESPONSIVE{1,sessioni,1};
%             b=intersect(a,Region_id);
%             Region_fraction_session{batchi}(sessioni,fishi)=length(b)/length(Region_id);
%         end

%                 figure,subplot(1,3,1);scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),5,'filled');hold on;
%                 scatter3(env.supervoxel(Region_id,1),env.supervoxel(Region_id,2),env.supervoxel(Region_id,3),10,'filled');axis equal
%                 subplot(1,3,2);scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,'filled');hold on;
%                 scatter3(supervolxeli(Region_id,1),supervolxeli(Region_id,2),supervolxeli(Region_id,3),10,'filled');axis equal
%                 subplot(1,3,3);scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),5,'filled');hold on;
%                 scatter3(supervolxeli(Region_id,1),supervolxeli(Region_id,2),supervolxeli(Region_id,3),10,'filled');axis equal
        %
        %         figure,
        %         scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,'filled');hold on;
        %         scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),5,'filled');hold on;
    end
end
save([savepath '\summary_Hbact_with_behavior.mat'],'num','Path','Region_name','Region_AUC_session','Region_AUC_session_all','Region_fraction_session',...
    'learning_speed','learning_perform','learning_trial',...
    'Region_AUC_preCS_session_all','Region_AUC_preCS_session','Region_AUC_CS_session_all','Region_AUC_CS_session','-v7.3');
%% plot
clc;clear all
savetype='jpg';
load('X:\calcium data 20230224\20210709\fish2\para.mat')
norm=load(['X:\calcium data 20230224\Hb_VS_behavior\Norm\', '\summary_Hbact_with_behavior.mat']);
norm_align=load(['X:\calcium data 20230224\',  '\align_trial.mat']);
Hb_act=load(['X:\calcium data 20230224\Hb_VS_behavior\Hb_act\', '\summary_Hbact_with_behavior.mat']);
Hb_align=load(['X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\', '\align_trial.mat']);
clr_group=[colorplus(389);colorplus(448);colorplus(233);colorplus(253)];
clr_group_bar=[clr_group(1,:);clr_group(4,:);[1 0.4 0.4];[0.4 1 0.4]];alpha=[1,1,1,1];
hhh_group_bar=categorical({'L','FL','L-Hb opto.','FL-Hb opto.'});hhh_group_bar=reordercats(hhh_group_bar,{'L','FL','L-Hb opto.','FL-Hb opto.'});
clr_cmap =addcolorplus(300);
group_label={'Learner(n=7)','Control(n=7)','Non-Learner(n=9)','Faded-Learner(n=4)'};
trial_base=trial.hab(2):trial.hab(3);
trial_ind_session=[];
isslidewin=1;
if isslidewin==1
    slidwin=3;
    sessiontype=['(slidewin=',num2str(slidwin),')'];
    for sessioni=1:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
        switch sessioni
            case 1
                trial_ind_session(:,sessioni)=trial.hab(2):trial.hab(3);
            case mat2cell([2:(trial.acq_block_num*trial.acq_block_trial-slidwin+1)+1],1,ones(1,(trial.acq_block_num*trial.acq_block_trial-slidwin)+1));
                trial_ind_session(1:slidwin,sessioni)=trial.hab(3)+sessioni-1 :trial.hab(3)+sessioni+slidwin-2;
            case (trial.acq_block_num*trial.acq_block_trial-slidwin+1)+2
                trial_ind_session(:,sessioni)=trial.test(2):trial.test(3);
        end
    end
else
    sessiontype='(block avg)';
    for sessioni=1:trial.acq_block_num+2
        switch sessioni
            case 1
                trial_ind_session(:,sessioni)=trial.hab(2):trial.hab(3)-3;
            case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num));
                trial_ind_session(1:trial.acq_block_trial,sessioni)=trial.hab(3)+(sessioni-1)*trial.acq_block_trial+1 :trial.hab(3)+(sessioni)*trial.acq_block_trial;
            case trial.acq_block_num+2
                trial_ind_session(:,sessioni)=trial.test(2):trial.test(3)-3;
        end
    end
end
norm.learning_quantity={};Hb_act.learning_quantity={};norm.learning_perform_session={};Hb_act.learning_perform_session={};
for batchi=[1,3,4]
    a=[];b=[];
    for sessioni=1:size(trial_ind_session,2)
        ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
        a(sessioni,:)=mean(norm.learning_perform{batchi}(ind,:),1,'omitnan');
        b(sessioni,:)=mean(Hb_act.learning_perform{batchi}(ind,:),1,'omitnan');
        
        norm.learning_perform_session{batchi}(sessioni,:)=mean(norm.learning_perform{batchi}(ind,:),1,'omitnan');
        Hb_act.learning_perform_session{batchi}(sessioni,:)=mean(Hb_act.learning_perform{batchi}(ind,:),1,'omitnan');
    end
    norm.learning_quantity{batchi}=max(a,[],1);Hb_act.learning_quantity{batchi}=max(b,[],1);
end
norm.learning_speed2={};Hb_act.learning_speed2={};norm.learning_trial2={};Hb_act.learning_trial2={};
for batchi=[1,3,4]
    a=norm_align.align_trial(batchi,:)+6-1;a(isnan(a))=trial.acq(1);a(find(a==0))=trial.acq(1);
    norm.learning_trial2{batchi}=a;b=(trial.acq(1)-a)./trial.acq(1);b(find(b==0))=0.05;norm.learning_speed2{batchi}=b;
    a=Hb_align.align_trial(batchi,:);a(isnan(a))=trial.acq(1);a(find(a==0))=trial.acq(1);
    Hb_act.learning_trial2{batchi}=a;b=(trial.acq(1)-a)./trial.acq(1);b(find(b==0))=0.05;Hb_act.learning_speed2{batchi}=b;
end
xlab{1}='Pre Cond.';xlab{size(trial_ind_session,2)}='Post Cond.';
for sessioni=2:size(trial_ind_session,2)-1
    xlab{sessioni}=strcat('Cond.',num2str(sessioni-1));
end



%% plot
%%%%%%%%%%%%%%check baseline%%%%%%%%%%%%%%%%%%%%%
for groupi=2
    b=nan(4,10);kk=1;
    for batchi=[1,4,3]
        switch groupi
            case 1
                c=norm.num{batchi};a=mean(norm.Region_AUC_session{batchi}(1:6,c),1,'omitnan');
            case 2
                c=Hb_act.num{batchi};a=mean(Hb_act.Region_AUC_session{batchi}(1:6,c),1,'omitnan');
        end
        b(kk,1:length(c)) =a;kk=kk+1;
    end
    h=[];p=[];
    for ii=1:4
        for jj=2:4
        [p(ii,jj),h(ii,jj)]=ranksum(b(ii,:),b(jj,:));
        end
    end
end

%%%%%%%%%%%%%%学习速率
h=figure('position',[265,126,300,300]);name=['Δ AUC VS' 'session'];
kk=1;x=nan(8,2);x2=nan(8,2);
for batchi=[1,4]
    a=norm.learning_speed2{batchi};b=Hb_act.learning_speed2{batchi};
    ylab='Learning speed';lim=[0 1];
    x(1:length(norm.num{batchi}),kk)=a(norm.num{batchi});
    x2(1:length(Hb_act.num{batchi}),kk)=b(Hb_act.num{batchi});
    kk=kk+1;
end
hh=cat(2,x,x2);%hh(find(hh==0))=0.05;
for i=1:size(hh,2)
    h=bar(i,mean(hh(:,i),1,'omitnan'),'FaceColor','flat','FaceAlpha',1,'BarWidth',1); hold on
    h.CData=clr_group_bar(i,:);h.FaceAlpha=alpha(i);
    d=hh(:,i);d(isnan(d))=[];
    scatter(i-0.15+rand(1,length(d)).*.1-.05,d,30,'filled','CData',[0 0 0],'MarkerEdgeAlpha',alpha(i)); hold on
end
ylabel(ylab,'fontsize',16,'FontWeight','bold'); set(gca,'XTick',[1:4],'Xticklabel',hhh_group_bar,'XTickLabelRotation',45);
set(gca,'XColor','k','YColor','k','linewidth',2);
set(gca,'FontName','Times New Roman', 'FontSize',13);
ylim(lim);
saveas(h,string(fullfile(savepath,name)),savetype);
%%%%%%%%%%%%%%Hb AUC in 8个session 
for typei=1
    h=figure('position',[265,126,300,300]);
    for indi=1
        ind=trial.acq(2):trial.acq(3);
        kk=1;x=nan(8,2);x2=nan(8,2);
        for batchi=[1,4]
            switch typei
                case 1
                    name=['AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
%                     pre=norm.Region_AUC_CS_session{batchi}; pre=mean(pre(trial_base,:),1,'omitnan'); a=(a-repmat(pre,size(a,1),1))./repmat(pre,size(a,1),1);
%                     pre=Hb_act.Region_AUC_CS_session{batchi};pre=mean(pre(trial_base,:),1,'omitnan');b=(b-repmat(pre,size(b,1),1))./repmat(pre,size(b,1),1);
                    ylab='AUC';lim=[0 1];
                case 2
                    name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                    ylab='Fraction';lim=[0 1];
            end
            aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
            x(1:length(norm.num{batchi}),kk)=aa;
            aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
            x2(1:length(Hb_act.num{batchi}),kk)=aa;
            kk=kk+1;
        end
        hh=cat(2,x,x2);
        for i=1:size(hh,2)
         h=bar(i,mean(hh(:,i),1,'omitnan'),'FaceColor','flat','FaceAlpha',1,'BarWidth',1); hold on
         h.CData=clr_group_bar(i,:);h.FaceAlpha=alpha(i);
         d=hh(:,i);d(isnan(d))=[];
         scatter(i-0.2+rand(1,length(d)).*.1-.05,d,30,'filled','CData',[0 0 0],'MarkerEdgeAlpha',alpha(i)); hold on
        end
        ylabel(ylab,'fontsize',16,'FontWeight','bold'); set(gca,'XTick',[1:4],'Xticklabel',hhh_group_bar,'XTickLabelRotation',45);
        set(gca,'XColor','k','YColor','k','linewidth',2);
        set(gca,'FontName','Times New Roman', 'FontSize',13);
        ylim(lim);
    end
end
for typei=1
    h=figure('position',[265,126,900,300]);
    for indi=1:3
        ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
        subplot(1,3,indi);kk=1;x=nan(8,2);x2=nan(8,2);
        for batchi=[1,4]
            switch typei
                case 1
                    name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
%                     pre=norm.Region_AUC_CS_session{batchi}; pre=mean(pre(trial_base,:),1,'omitnan'); a=(a-repmat(pre,size(a,1),1))./repmat(pre,size(a,1),1);
%                     pre=Hb_act.Region_AUC_CS_session{batchi};pre=mean(pre(trial_base,:),1,'omitnan');b=(b-repmat(pre,size(b,1),1))./repmat(pre,size(b,1),1);
                    ylab='Δ AUC';lim=[0 1];
                case 2
                    name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                    ylab='Fraction';lim=[0 1];
            end
            aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
            x(1:length(norm.num{batchi}),kk)=aa;
            aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
            x2(1:length(Hb_act.num{batchi}),kk)=aa;
            kk=kk+1;
        end
        hh=cat(2,x,x2);
        for i=1:size(hh,2)
         h=bar(i,mean(hh(:,i),1,'omitnan'),'FaceColor','flat','FaceAlpha',1,'BarWidth',1); hold on
         h.CData=clr_group_bar(i,:);h.FaceAlpha=alpha(i);
         d=hh(:,i);d(isnan(d))=[];
         scatter(i-0.15+rand(1,length(d)).*.1-.05,d,30,'filled','CData',[0 0 0],'MarkerEdgeAlpha',alpha(i)); hold on
        end
        ylabel(ylab,'fontsize',16,'FontWeight','bold'); set(gca,'XTick',[1:4],'Xticklabel',hhh_group_bar,'XTickLabelRotation',45);
        set(gca,'XColor','k','YColor','k','linewidth',2);
        set(gca,'FontName','Times New Roman', 'FontSize',13);
        ylim(lim);
        %saveas(h,string(fullfile(savepath,name)),savetype);
    end
end
%%%%%%%%%%%%%%Hb AUC in 3个align session 
for typei=1
    h=figure('position',[265,126,900,300]);
    for indi=1:3
        subplot(1,3,indi);
        kk=1;x=nan(8,2);x2=nan(8,2);
        for batchi=[1,4]
            switch typei
                case 1
                    name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                    ylab='Δ AUC';lim=[0 1];
                case 2
                    name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                    ylab='Fraction';lim=[0 1];
            end
            iii=1;aa=[];
            for ii=norm.num{batchi}
                ind=norm.learning_trial2{batchi}(ii);
                switch indi
                    case 1
                        aa(iii)=mean(a(trial.acq(2):trial.acq(2)+ind-1,ii),1,'omitnan');
                    case 2
                        aa(iii)=mean(a(trial.acq(2)+ind-1,ii),1,'omitnan');
                    case 3
                        aa(iii)=mean(a(trial.acq(2)+ind-1:trial.acq(3),ii),1,'omitnan');
                end
                iii=iii+1;
            end
            x(1:length(norm.num{batchi}),kk)=aa;
            iii=1;aa=[];
            for ii=Hb_act.num{batchi}
                ind=Hb_act.learning_trial2{batchi}(ii);
                switch indi
                    case 1
                        aa(iii)=mean(b(trial.acq(2):trial.acq(2)+ind-1,ii),1,'omitnan');
                    case 2
                        aa(iii)=mean(b(trial.acq(2)+ind-1,ii),1,'omitnan');
                    case 3
                        aa(iii)=mean(b(trial.acq(2)+ind-1:trial.acq(3),ii),1,'omitnan');
                end
                iii=iii+1;
            end
            x2(1:length(Hb_act.num{batchi}),kk)=aa;
            kk=kk+1;
        end
        hh=cat(2,x,x2);
        for i=1:size(hh,2)
         h=bar(i,mean(hh(:,i),1,'omitnan'),'FaceColor','flat','FaceAlpha',1,'BarWidth',1); hold on
         h.CData=clr_group_bar(i,:);h.FaceAlpha=alpha(i);
         d=hh(:,i);d(isnan(d))=[];
         scatter(i-0.15+rand(1,length(d)).*.1-.05,d,20,'filled','CData',[0 0 0],'MarkerEdgeAlpha',alpha(i)); hold on
        end
        ylabel(ylab,'fontsize',16,'FontWeight','bold'); set(gca,'XTick',[1:4],'Xticklabel',hhh_group_bar,'XTickLabelRotation',45);
        set(gca,'XColor','k','YColor','k','linewidth',2);
        set(gca,'FontName','Times New Roman', 'FontSize',13);
        ylim(lim);
    end
end
%%%% cdf of Hb AUC  across trial
h=figure('position',[265,126,600,600]);  name=['norm AUC cdf'];
for batchi=[1,4]
    x=[];x2=[];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi}; lim=[0 1.5];
    for sessioni=1:size(trial_ind_session,2)
        ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
        x(:,sessioni)=mean(a(ind,norm.num{batchi}),1,'omitnan');
        x2(:,sessioni)=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');
    end
    %x=a(trial.acq(2):trial.acq(3),norm.num{batchi});x2=a(trial.acq(2):trial.acq(3),Hb_act.num{batchi});
%     pre=mean(a(trial_base,norm.num{batchi}),1,'omitnan')';pre2=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan')';
%     x=(x-repmat(pre,1,size(x,2)))./(repmat(pre,1,size(x,2)));x2=(x2-repmat(pre2,1,size(x2,2)))./(repmat(pre2,1,size(x2,2)));
    %x=cumsum(x,2); x2=cumsum(x2,2);
    y=norm.learning_trial2{batchi}(norm.num{batchi});y2= Hb_act.learning_trial2{batchi}(Hb_act.num{batchi});%if batchi==1 y(y==0)=0.05; y2(y2==0)=0.05;end;
    yy= norm.learning_quantity{batchi}(norm.num{batchi}); yy2= Hb_act.learning_quantity{batchi}(Hb_act.num{batchi});
    
    plot(1:size(x,2),x,'color',clr_group(batchi,:),'linewidth',2);hold on;
    for ii=1:length(y)
    scatter(y(ii),x(ii,y(ii)),30,clr_cmap(ceil(yy(ii).*(size(clr_cmap,1)-2))+1,:),'filled');hold on;
    end
    plot(1:size(x,2),x2,'color',clr_group(batchi,:),'linewidth',2,'linestyle','--');hold on;
    for ii=1:length(y2)
        scatter(y2(ii),x2(ii,y2(ii)),30,clr_cmap(ceil(yy2(ii).*(size(clr_cmap,1)-2))+1,:),'filled');hold on;
    end
    ylabel('Δ AUC','fontsize',16,'FontWeight','bold');
    xlabel('Trial','fontsize',16,'FontWeight','bold');
    set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
    set(gca,'FontName','Times New Roman');
    xlim([1 size(x,2)]);ylim(lim)
    colormap(clr_cmap);colorbar;
    %saveas(h,string(fullfile(savepath,name)),savetype);
end
for type=1
h=figure('position',[265,126,600,600]);  name=['norm AUC cdf'];colormap(clr_cmap);%clim([0 0.8]);
for batchi=[1,4]
    x=[];x2=[];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi}; lim=[0 20];
    for sessioni=1:size(trial_ind_session,2)
        ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
        x(:,sessioni)=mean(a(ind,norm.num{batchi}),1,'omitnan');
        x2(:,sessioni)=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');
    end
    %     pre=mean(a(trial_base,norm.num{batchi}),1,'omitnan')';pre2=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan')';
    %     x=(x-repmat(pre,1,size(x,2)))./(repmat(pre,1,size(x,2)));x2=(x2-repmat(pre2,1,size(x2,2)))./(repmat(pre2,1,size(x2,2)));
    x=cumsum(x,2); x2=cumsum(x2,2);
    y=norm.learning_trial2{batchi}(norm.num{batchi});y2= Hb_act.learning_trial2{batchi}(Hb_act.num{batchi});%if batchi==1 y(y==0)=0.05; y2(y2==0)=0.05;end;
    yy= norm.learning_quantity{batchi}(norm.num{batchi}); yy2= Hb_act.learning_quantity{batchi}(Hb_act.num{batchi});
    yyy= norm.learning_speed2{batchi}(norm.num{batchi}); yyy2= Hb_act.learning_speed2{batchi}(Hb_act.num{batchi});
    for ii=1:length(y)
        plot(1:size(trial_ind_session,2),x(ii,:),'color',clr_cmap(ceil(yyy(ii).*(size(clr_cmap,1)-2))+1,:),'linewidth',2);hold on;
        scatter(y(ii),x(ii,y(ii)),30,clr_cmap(ceil(yy(ii).*(size(clr_cmap,1)-2))+1,:),'filled');hold on;
    end
    for ii=1:length(y2)
        plot(1:size(trial_ind_session,2),x2(ii,:),'color',clr_cmap(ceil(yyy2(ii).*(size(clr_cmap,1)-2))+1,:),'linewidth',2,'linestyle','--');hold on;
        scatter(y2(ii),x2(ii,y2(ii)),30,clr_cmap(ceil(yy2(ii).*(size(clr_cmap,1)-2))+1,:),'filled');hold on;
    end
end
ylabel('AUC','fontsize',16,'FontWeight','bold');
xlabel('Trial','fontsize',16,'FontWeight','bold');
set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
set(gca,'FontName','Times New Roman');
xlim([1 size(trial_ind_session,2)]);ylim(lim)
colorbar;
%saveas(h,string(fullfile(savepath,name)),savetype);
end
for type=1
    h=figure('position',[265,126,600,600]);  name=['norm AUC cdf'];colormap(clr_cmap);%clim([0 0.8]);
    kk=1;p=[];p2=[];
    for batchi=[1,4]
        x=[];x2=[];a=norm.Region_AUC_CS_session{batchi};b=Hb_act.Region_AUC_CS_session{batchi}; lim=[0 40000];
        for sessioni=1:size(trial_ind_session,2)
            ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
            x(:,sessioni)=mean(a(ind,norm.num{batchi}),1,'omitnan');
            x2(:,sessioni)=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');
        end
        x=cumsum(x,2); x2=cumsum(x2,2);
        y=norm.learning_trial2{batchi}(norm.num{batchi});y2= Hb_act.learning_trial2{batchi}(Hb_act.num{batchi});%if batchi==1 y(y==0)=0.05; y2(y2==0)=0.05;end;
        yy= norm.learning_quantity{batchi}(norm.num{batchi}); yy2= Hb_act.learning_quantity{batchi}(Hb_act.num{batchi});
        yyy= norm.learning_speed2{batchi}(norm.num{batchi}); yyy2= Hb_act.learning_speed2{batchi}(Hb_act.num{batchi});
        shadedErrorBar(1:size(trial_ind_session,2),mean(x,1),std(x,[],1,'omitnan')./sqrt(size(x,1)),'lineprops',{clr_cmap(ceil(mean(yyy).*(size(clr_cmap,1)+30))+1,:)},'transparent',1,'patchSaturation',0.1); hold on
        p(kk)=plot(1:size(trial_ind_session,2),mean(x,1),'color',clr_cmap(ceil(mean(yyy).*(size(clr_cmap,1)+30))+1,:),'linewidth',2);hold on;
        shadedErrorBar(1:size(trial_ind_session,2),mean(x2,1),std(x2,[],1,'omitnan')./sqrt(size(x2,1)),'lineprops',{clr_cmap(ceil(mean(yyy2).*(size(clr_cmap,1)+30))+1,:)},'transparent',1,'patchSaturation',0.1); hold on
        p2(kk)=plot(1:size(trial_ind_session,2),mean(x2,1),'color',clr_cmap(ceil(mean(yyy2).*(size(clr_cmap,1)+30))+1,:),'linewidth',2,'linestyle','--');hold on;
        kk=kk+1;
    end
    ylabel('Δ AUC','fontsize',16,'FontWeight','bold');
    xlabel('Trial','fontsize',16,'FontWeight','bold');
    set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
    set(gca,'FontName','Times New Roman');
    legend([p,p2],{'L','FL','L-Hb opto.','FL-Hb opto.'});
    xlim([1 size(trial_ind_session,2)]);ylim(lim)
    colorbar;
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%%%%%%%%%%%%%%Hb AUC in 8个session VS learning speed
for performtypei=1
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,300,300]);
        for indi=1
            ind=trial.acq(2):trial.acq(3);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                        xlab='AUC';lim=[0 1];
                    case 2
                        name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                        xlab='Fraction';lim=[0 1];
                end
                aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
                x(1:length(norm.num{batchi}),kk)=aa;
                aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
                x2(1:length(Hb_act.num{batchi}),kk)=aa;
                kk=kk+1;
            end
            hh=cat(2,x,x2);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for  batchi=[1,4]
                switch  performtypei
                    case 1
                        a=norm.learning_speed2{batchi};b=Hb_act.learning_speed2{batchi};ylab='Learning speed';
                    case 2
                        a=norm.learning_quantity{batchi};b=Hb_act.learning_quantity{batchi};ylab='Learning quantity';
                end
                x(1:length(norm.num{batchi}),kk)=a(norm.num{batchi});
                x2(1:length(Hb_act.num{batchi}),kk)=b(Hb_act.num{batchi});
                kk=kk+1;
            end
            hh2=cat(2,x,x2);
            for i=1:size(hh,2)
                d=hh(:,i);d(isnan(d))=[]; d2=hh2(:,i);d2(isnan(d2))=[];
                scatter(d,d2,30,'filled','CData',clr_group_bar(i,:),'MarkerFaceAlpha',alpha(i)); hold on
            end
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            ylabel(ylab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
for performtypei=1
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,300,300]);
        for indi=1
            ind=trial.acq(2):trial.acq(3);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                        xlab='AUC';lim=[0 1];
                    case 2
                        name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                        xlab='Fraction';lim=[0 1];
                end
                aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
                x(1:length(norm.num{batchi}),kk)=aa;
                aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
                x2(1:length(Hb_act.num{batchi}),kk)=aa;
                kk=kk+1;
            end
            hh=cat(2,x,x2);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for  batchi=[1,4]
                switch  performtypei
                    case 1
                        a=norm.learning_speed2{batchi};b=Hb_act.learning_speed2{batchi};ylab='Learning speed';
                    case 2
                        a=norm.learning_quantity{batchi};b=Hb_act.learning_quantity{batchi};ylab='Learning quantity';
                end
                x(1:length(norm.num{batchi}),kk)=a(norm.num{batchi});
                x2(1:length(Hb_act.num{batchi}),kk)=b(Hb_act.num{batchi});
                kk=kk+1;
            end
            hh2=cat(2,x,x2);
            for i=1:size(hh,2)
                d=hh(:,i);d(isnan(d))=[]; d2=hh2(:,i);d2(isnan(d2))=[];
                scatter(mean(d),mean(d2),30,'filled','CData',clr_group_bar(i,:),'MarkerFaceAlpha',alpha(i)); hold on
            end
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            ylabel(ylab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
for performtypei=1
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,900,300]);
        for indi=1:3
            ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
            subplot(1,3,indi);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                        xlab='AUC';lim=[0 1];
                    case 2
                        name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                        xlab='Fraction';lim=[0 1];
                end
                aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
                x(1:length(norm.num{batchi}),kk)=aa;
                aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
                x2(1:length(Hb_act.num{batchi}),kk)=aa;
                kk=kk+1;
            end
            hh=cat(2,x,x2);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for  batchi=[1,4]
                switch  performtypei
                    case 1
                        a=norm.learning_speed2{batchi};b=Hb_act.learning_speed2{batchi};ylab='Learning speed';
                    case 2
                        a=norm.learning_quantity{batchi};b=Hb_act.learning_quantity{batchi};ylab='Learning quantity';
                end
                x(1:length(norm.num{batchi}),kk)=a(norm.num{batchi});
                x2(1:length(Hb_act.num{batchi}),kk)=b(Hb_act.num{batchi});
                kk=kk+1;
            end
            hh2=cat(2,x,x2);
            for i=1:size(hh,2)
                d=hh(:,i);d(isnan(d))=[]; d2=hh2(:,i);d2(isnan(d2))=[];
                scatter(mean(d),mean(d2),30,'filled','CData',clr_group_bar(i,:),'MarkerFaceAlpha',alpha(i)); hold on
            end
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            ylabel(ylab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
for performtypei=1:2
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,900,300]);
        for indi=1:3
            ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
            subplot(1,3,indi);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                        xlab='AUC';lim=[0 1];
                    case 2
                        name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                        xlab='Fraction';lim=[0 1];
                end
                aa=mean(a(ind,norm.num{batchi}),1,'omitnan');
                x(1:length(norm.num{batchi}),kk)=aa;
                aa=mean(b(ind,Hb_act.num{batchi}),1,'omitnan');%pre=mean(b(trial_base,Hb_act.num{batchi}),1,'omitnan');aa=(aa-pre)./pre;
                x2(1:length(Hb_act.num{batchi}),kk)=aa;
                kk=kk+1;
            end
            hh=cat(2,x,x2);
            kk=1;x=nan(8,2);x2=nan(8,2);
            for  batchi=[1,4]
                switch  performtypei
                    case 1
                        a=norm.learning_speed2{batchi};b=Hb_act.learning_speed2{batchi};ylab='Learning speed';
                    case 2
                        a=norm.learning_quantity{batchi};b=Hb_act.learning_quantity{batchi};ylab='Learning quantity';
                end
                x(1:length(norm.num{batchi}),kk)=a(norm.num{batchi});
                x2(1:length(Hb_act.num{batchi}),kk)=b(Hb_act.num{batchi});
                kk=kk+1;
            end
            hh2=cat(2,x,x2);
            for i=1:size(hh,2)
                d=hh(:,i);d(isnan(d))=[]; d2=hh2(:,i);d2(isnan(d2))=[];
                scatter(d,d2,30,'filled','CData',clr_group_bar(i,:),'MarkerFaceAlpha',alpha(i)); hold on
            end
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            ylabel(ylab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
%%%%%%%%%%%%%%Hb AUC in 3个algin session VS learning speed
for performtypei=1
    switch  performtypei
        case 1
            ylab='Learning speed';
    end
    for typei=1
        h=figure('position',[265,126,900,300]);
        for indi=1
            subplot(1,3,indi);
            kk=1;x=nan(10,2);x2=nan(10,2);
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['Δ AUC VS' 'session'];a=norm.Region_AUC_session{batchi};b=Hb_act.Region_AUC_session{batchi};
                        xlab='AUC';lim=[0 1];
                    case 2
                        name=['Fracrion VS','session'];a=norm.Region_fraction_session{batchi};b=Hb_act.Region_fraction_session{batchi};
                        xlab='Fraction';lim=[0 1];
                end
                iii=1;aa=[];
                for ii=norm.num{batchi}
                    ind=norm.learning_trial2{batchi}(ii);
                    switch indi
                        case 1
                            aa(iii)=mean(a(trial.acq(2):trial.acq(2)+ind-1,ii),1,'omitnan');
                        case 2
                            aa(iii)=mean(a(trial.acq(2)+ind-1,ii),1,'omitnan');
                        case 3
                            aa(iii)=mean(a(trial.acq(2)+ind-1:trial.acq(3),ii),1,'omitnan');
                    end
                    iii=iii+1;
                end
                x(1:length(norm.num{batchi}),kk)=aa;
                iii=1;aa=[];
                for ii=Hb_act.num{batchi}
                    ind=Hb_act.learning_trial2{batchi}(ii);
                    switch indi
                        case 1
                            aa(iii)=mean(b(trial.acq(2):trial.acq(2)+ind-1,ii),1,'omitnan');
                        case 2
                            aa(iii)=mean(b(trial.acq(2)+ind-1,ii),1,'omitnan');
                        case 3
                            aa(iii)=mean(b(trial.acq(2)+ind-1:trial.acq(3),ii),1,'omitnan');
                    end
                    iii=iii+1;
                end
                x2(1:length(Hb_act.num{batchi}),kk)=aa;
                kk=kk+1;
            end
            hh=cat(2,x,x2);
            kk=1;x=nan(10,2);x2=nan(10,2);
            for  batchi=[1,4]
                x(1:length(norm.num{batchi}),kk)=norm.learning_speed2{batchi}(norm.num{batchi}).*norm.learning_quantity{batchi}(norm.num{batchi});
                x2(1:length(Hb_act.num{batchi}),kk)=Hb_act.learning_speed2{batchi}(Hb_act.num{batchi}).*Hb_act.learning_quantity{batchi}(Hb_act.num{batchi});
                kk=kk+1;
            end
            hh2=cat(2,x,x2);
            for i=1:size(hh,2)
                d=hh(:,i);d(isnan(d))=[]; d2=hh2(:,i);d2(isnan(d2))=[];
                scatter(d,d2,30,'filled','CData',clr_group_bar(i,:),'MarkerFaceAlpha',alpha(i)); hold on
            end
            a=reshape(hh,1,[]);a=a(~isnan(a))';b=reshape(hh2,1,[]);b=b(~isnan(b))';c=cat(2,a,b);
            figure,Formu=mylineplotofregression(reshape(hh,1,[]),reshape(hh2,1,[]),[],[],group_label(batchi));hold on;
            text(0.9*max(lim),1,Formu,'color','k','FontSize',12,'FontName','Times New Roman','FontWeight','bold');
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            ylabel(ylab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
%%%%%%%%%%%% Hbact VS performance
for groupi=2
    for typei=1
        h=figure('position',[265,126,973,769]);
        for batchi=[1,3,4]
            x=[];y=[];
            switch typei
                case 1
                    name='AUC VS perform';
                    switch groupi
                        case 1
                            a=norm.Region_AUC_session{batchi};b=norm.num{batchi};c=norm.learning_perform{batchi};
                        case 2
                            a=Hb_act.Region_AUC_session{batchi};b=Hb_act.num{batchi};c=Hb_act.learning_perform{batchi};
                    end
                    xlab='AUC';lim=[0 0.5];
                case 2
                    name='Fraction VS perform';
                    switch groupi
                        case 1
                            a=norm.Region_fraction_session{batchi};
                        case 2
                            a=Hb_act.Region_fraction_session{batchi};
                    end
                    xlab='Fraction';lim=[0 1];
            end
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
                x(sessioni,:)=mean(a(ind,b),1,'omitnan');
                y(sessioni,:)=mean(c(ind,b),1,'omitnan');
            end
            subplot(2,2,batchi),
            Formu=mylineplotofregression(reshape(x,1,[]),reshape(y,1,[]),[],[],group_label(batchi));hold on;
            text(0.9*max(lim),1,Formu,'color','k','FontSize',12,'FontName','Times New Roman','FontWeight','bold');
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
            ylabel(['CR ratio',sessiontype],'fontsize',16,'FontWeight','bold');
            xlabel(xlab,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
            xlim(lim);ylim([0 1])
        end
    end
end
for typei=1
    h=figure('position',[265,126,973,769]);x={};y={};
    for groupi=1:2
        for batchi=[1,4]
            switch groupi
                case 1
                    a=norm.Region_AUC_session{batchi};b=norm.num{batchi};c=norm.learning_perform{batchi};
                case 2
                    a=Hb_act.Region_AUC_session{batchi};b=Hb_act.num{batchi};c=Hb_act.learning_perform{batchi};
            end
            xlab='AUC';lim=[0 1];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
                x{groupi,batchi}(sessioni,:)=mean(a(ind,b),1,'omitnan');
                y{groupi,batchi}(sessioni,:)=mean(c(ind,b),1,'omitnan');
            end
        end
    end
    xx=cat(2,x{1,1},x{1,4},x{2,1},x{2,4}); yy=cat(2,y{1,1},y{1,4},y{2,1},y{2,4});a=reshape(xx,1,[]);b=reshape(yy,1,[]);
    Formu=mylineplotofregression(reshape(xx,1,[]),reshape(yy,1,[]),[],[],[]);hold on;
    text(0.9*max(lim),1,Formu,'color','k','FontSize',12,'FontName','Times New Roman','FontWeight','bold');
    scatter(reshape(x{1,1},1,[]),reshape(y{1,1},1,[]),30,clr_group_bar(1,:),'filled');hold on;
    scatter(reshape(x{1,4},1,[]),reshape(y{1,4},1,[]),30,clr_group_bar(2,:),'filled');hold on;
    scatter(reshape(x{2,1},1,[]),reshape(y{2,1},1,[]),30,clr_group_bar(3,:),'filled');hold on;
    scatter(reshape(x{2,4},1,[]),reshape(y{2,4},1,[]),30,clr_group_bar(4,:),'filled');hold on;
    ylabel(['CR ratio',sessiontype],'fontsize',16,'FontWeight','bold');
    xlabel(xlab,'fontsize',16,'FontWeight','bold');
    set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
    xlim(lim);ylim([0 1])
end
%%%%%%%%%single fish（时间轴是session)
for groupi=1:2
    for batchi=[1,4]
        switch groupi
            case 1
                a=norm.Region_AUC_session{batchi};b=norm.num{batchi};c=norm.learning_perform_session{batchi};d=norm.Path{batchi}; ylab='AUC';lim=[0 0.4];
                e=norm.Region_AUC_preCS_session{batchi};
            case 2
                a=Hb_act.Region_AUC_session{batchi};b=Hb_act.num{batchi};c=Hb_act.learning_perform_session{batchi};d=Hb_act.Path{batchi}; ylab='AUC';lim=[0 0.8];
        end
        h=figure('position',[33,30,1900,ceil(length(b)/4)*250],'color','w');t=tiledlayout(ceil(length(b)/4),4);t.TileSpacing = 'compact';t.Padding = 'compact';
        for fishi=b
            path=d{fishi};nn=[path(end-14:end-7),path(end-5:end-1)];
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
                y1(sessioni)=mean(a(ind,fishi),1,'omitnan');sd1(sessioni)=std(a(ind,fishi),[],1,'omitnan');
            end
            x=1:size(trial_ind_session,2);y2=c(:,fishi);
            yyaxis left;shadedErrorBar(x,y1,sd1,'lineprops',{[1 0 0]},'transparent',1,'patchSaturation',0.1); hold on
            ylim(lim);
            ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
            yyaxis right;
            plot(x, y2,'k-.','linewidth',1.5,'markersize',8);hold on
            ylabel('CR ratio','fontsize',16,'FontWeight','bold');ylim([0 1]);set(gca,'ycolor',[0 0 0],'linewidth',2);
            xlim([0 size(trial_ind_session,2)]);
            set(gca,'xtick',x,'xticklabel',xlab,'XTickLabelRotation',45,'fontsize',12);
            title(nn,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        end
    end
end
for groupi=1:2
    for batchi=[1,4]
        switch groupi
            case 1
                a=norm.Region_AUC_session{batchi};b=norm.num{batchi};c=norm.learning_perform{batchi};d=norm.Path{batchi}; ylab='AUC';lim=[0 0.4];
            case 2
                a=Hb_act.Region_AUC_session{batchi};b=Hb_act.num{batchi};c=Hb_act.learning_perform{batchi};d=Hb_act.Path{batchi}; ylab='AUC';lim=[0 0.8];
        end
        h=figure('position',[33,30,1900,ceil(length(b)/4)*250],'color','w');t=tiledlayout(ceil(length(b)/4),4);t.TileSpacing = 'compact';t.Padding = 'compact';
        for fishi=b
            path=d{fishi};nn=[path(end-14:end-7),path(end-5:end-1)];
            load(fullfile(path, '/activities_dfdf_align.mat'),'align_win');
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(align_win,2)
                ind=align_win(:,sessioni);ind(isnan(ind))=[];
                y1(sessioni)=mean(a(ind,fishi),1,'omitnan');sd1(sessioni)=std(a(ind,fishi),[],1,'omitnan');
                y2(sessioni)=mean(c(ind,fishi),1,'omitnan');
            end
            x=1:8;
            yyaxis left;shadedErrorBar(x,y1,sd1,'lineprops',{[1 0 0]},'transparent',1,'patchSaturation',0.1); hold on
            ylim(lim);
            ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
            yyaxis right;
            plot(x, y2,'k-.','linewidth',1.5,'markersize',8);hold on
            ylabel('CR ratio','fontsize',16,'FontWeight','bold');ylim([0 1]);set(gca,'ycolor',[0 0 0],'linewidth',2);
            xlim([0 9]);
            set(gca,'xtick',x,'xticklabel',{'Pre Cond.','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond.'},'XTickLabelRotation',45,'fontsize',12);
            title(nn,'fontsize',16,'FontWeight','bold');
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        end
    end
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%old version%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Avg Hb  AUC(across fish) in  align session 
for typei=1
    h=figure('position',[265,126,973,769]);
    for batchi=[1:length(Path)]
        subplot(2,2,batchi),
        switch typei
            case 1
                name=['mean AUC(acorss fish)'];a=Region_AUC_session{batchi};ylab='AUC';lim=[0 3];
            case 2
                name=['mean Fraction(acorss fish)'];a=Region_fraction_session{batchi};ylab='Fraction';lim=[0 1];
        end
        y1=[]; y2=[];kk=1;
        for fishi=num{batchi}
            path=Path{batchi}{fishi};
            load(fullfile(path, '/activities_dfdf_align.mat'),'align_win');
            for sessioni=1:size(align_win,2)
                ind=align_win(:,sessioni);ind(isnan(ind))=[];
                y1(sessioni,kk)=mean(a(ind,fishi),1,'omitnan');
                y2(sessioni,kk)=mean(learning_perform{batchi}(ind,fishi),1,'omitnan');
            end
        kk=kk+1;
        end
        x=1:size(align_win,2);
        yyaxis left;
        shadedErrorBar(x,mean(y1,2,'omitnan'),std(y1,[],2,'omitnan')./sqrt(length(num{batchi})),'lineprops',{[1 0 0]},'transparent',1,'patchSaturation',0.1); hold on
        ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
        ylim(lim);yyaxis right;
        plot(x,mean(y2,2,'omitnan'),'k--'); hold on
        ylabel('CR ratio (aligned)','fontsize',16,'FontWeight','bold');ylim([0 1.2]);set(gca,'ycolor',[0 0 0],'linewidth',2);
        ylim([0 1]);xlim([1 size(align_win,2)])
        set(gca,'XTickLabelRotation',45,'fontsize',12);
        title(group_label{batchi},'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%%%%%%%%%single fish（时间轴是session)
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
            path=Path{batchi}{fishi};nn=[path(end-14:end-7),path(end-5:end-1)];
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
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
%%%%%%%%%%single fish（时间轴是session,按第一个session，normalize)
for typei=1
    for batchi=1:length(Path)
        switch typei
            case 1
                name=strcat(group_label(batchi),'_normAUC of HbmeanCSact');a=Region_AUC_session{batchi};ylab='Δ AUC';
            case 2
                name=strcat(group_label(batchi),'_normFraction of HbCSactivation');a=Region_fraction_session{batchi};ylab='norm.Fraction';
        end
        h=figure('position',[33,30,1900,ceil(length(num{batchi})/4)*250],'color','w');t=tiledlayout(ceil(length(num{batchi})/4),4);t.TileSpacing = 'compact';t.Padding = 'compact';
        for fishi=num{batchi}
            path=Path{batchi}{fishi};nn=[path(end-14:end-7),path(end-5:end-1)];
            ax1=nexttile([1,1]);y1=[];y2=[];sd1=[];
            for sessioni=1:size(trial_ind_session,2)
                ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
                y1(sessioni)=mean(a(ind,fishi),1,'omitnan');
                y2(sessioni)=mean(learning_perform{batchi}(ind,fishi),1,'omitnan');
            end
            x=1:size(trial_ind_session,2);y1=(y1-y1(1))./(y1(1));%sd1=(sd1-sd1(1))./(sd1(1));
            yyaxis left;plot(x, y1,'r','linewidth',1.5,'markersize',8);hold on
            ylabel(ylab,'fontsize',16,'FontWeight','bold');set(gca,'ycolor',[1 0 0],'linewidth',2);
            ylim([-1 2]);
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
%%%%%%%%%%%% Hbact VS performance
for typei=1:2
    h=figure('position',[265,126,973,769]);
    for batchi=[1:length(Path)]
        x=[];y=[];
        switch typei
            case 1
                name='AUC VS perform';a=Region_AUC_session{batchi};xlab='AUC';lim=[0 500];
            case 2
                name='Fraction VS perform';a=Region_fraction_session{batchi};xlab='Fraction';lim=[0 1];
        end
        for sessioni=1:size(trial_ind_session,2)
            ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
            x(sessioni,:)=mean(a(ind,num{batchi}),1,'omitnan');
            y(sessioni,:)=mean(learning_perform{batchi}(ind,num{batchi}),1,'omitnan');
        end
        subplot(2,2,batchi),
        Formu=mylineplotofregression(reshape(x,1,[]),reshape(y,1,[]),[],[],group_label(batchi));hold on;
        text(0.9*max(lim),1,Formu,'color','k','FontSize',12,'FontName','Times New Roman','FontWeight','bold');
        scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        ylabel(['CR ratio',sessiontype],'fontsize',16,'FontWeight','bold');
        xlabel(xlab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%%%%%%%%%%%%% Hbact VS performance(normalize)
for typei=1:2
    h=figure('position',[265,126,973,769]);
    for batchi=[1:length(Path)]
        x=[];y=[];
        switch typei
            case 1
                name='normAUC VS perform';a=Region_AUC_session{batchi};xlab='Δ AUC';lim=[-1 2];
            case 2
                name='normFraction VS perform';a=Region_fraction_session{batchi};xlab='Δ Fraction';lim=[0 1];
        end
        for sessioni=1:size(trial_ind_session,2)
            ind=trial_ind_session(:,sessioni);ind(find(ind==0))=[];
            x(sessioni,:)=mean(a(ind,num{batchi}),1,'omitnan');
            y(sessioni,:)=mean(learning_perform{batchi}(ind,num{batchi}),1,'omitnan');
        end
        subplot(2,2,batchi),x=(x-repmat(x(1,:),size(trial_ind_session,2),1))./(repmat(x(1,:),size(trial_ind_session,2),1));
        Formu=mylineplotofregression(reshape(x,1,[]),reshape(y,1,[]),[],[],group_label(batchi));hold on;
        text(0.9*max(lim),1,Formu,'color','k','FontSize',12,'FontName','Times New Roman','FontWeight','bold');
        scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        ylabel(['CR ratio',sessiontype],'fontsize',16,'FontWeight','bold');
        xlabel(xlab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%%%%%%%%%%%%%%Hb AUC in 8个session VS learning speed
for performtypei=1:2
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,900,300]);
        for indi=1:3
            ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
            subplot(1,3,indi)
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['AUC VS',ylab, '_session'];a=Region_AUC_session{batchi};xlab='AUC';lim=[0 500];
                    case 2
                        name=['Fracrion VS',ylab, '_session'];a=Region_fraction_session{batchi};xlab='Fraction';lim=[0 1];
                end
                x=mean(a(ind,num{batchi}),1,'omitnan');
                switch  performtypei
                    case 1
                        y= learning_speed{batchi}(num{batchi});if batchi==1 y(y==0)=0.05; end;
                    case 2
                        y= learning_quantity{batchi}(num{batchi});
                end
                scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
                ylabel(ylab,'fontsize',16,'FontWeight','bold');
                xlabel(xlab,'fontsize',16,'FontWeight','bold');
                set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
                set(gca,'FontName','Times New Roman');
                xlim(lim);ylim([0 1])
            end
            saveas(h,string(fullfile(savepath,name)),savetype);
        end
    end
end
%%%%%%%%%%%%%%%%%%%Hb AUC(norm) in 8个session VS learning speed
for performtypei=1:2
    switch  performtypei
        case 1
            ylab='Learning speed';
        case 2
            ylab='Learning quantity';
    end
    for typei=1
        h=figure('position',[265,126,900,300]);
        for indi=1:3
            ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
            subplot(1,3,indi)
            for batchi=[1,4]
                switch typei
                    case 1
                        name=['normAUC VS',ylab, '_session'];a=Region_AUC_session{batchi};xlab='Δ AUC';lim=[-1 2];
                    case 2
                        name=['normFracrion VS',ylab, '_session'];a=Region_fraction_session{batchi};xlab='Δ Fraction';lim=[0 1];
                end
                x=mean(a(ind,num{batchi}),1,'omitnan');
                pre=mean(a(trial.hab(2):trial.hab(3),num{batchi}),1,'omitnan');
                x=(x-pre)./pre;
                switch  performtypei
                    case 1
                        y= learning_speed{batchi}(num{batchi});if batchi==1 y(y==0)=0.05; end;
                    case 2
                        y= learning_quantity{batchi}(num{batchi});
                end
                scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
                ylabel(ylab,'fontsize',16,'FontWeight','bold');
                xlabel(xlab,'fontsize',16,'FontWeight','bold');
                set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
                set(gca,'FontName','Times New Roman');
                xlim(lim);ylim([0 1])
            end
            saveas(h,string(fullfile(savepath,name)),savetype);
        end
    end
end
    %Hb mean AUC in  all session VS learning speed
for typei=1:2
    h=figure('position',[265,126,300,300]);
    for indi=1
        ind=trial.acq(2):trial.acq(3);
        for batchi=[1:3]
            switch typei
                case 1
                    name=['AUC VS learningspeed_session_all'];a=Region_AUC_session{batchi};ylab='AUC';lim=[0 25];
                case 2
                    name=['Fraction VS learningspeed_session_all'];a=Region_fraction_session{batchi};ylab='Fraction';lim=[0 1];
            end
            y= learning_speed{batchi}(num{batchi});x=a(ind,num{batchi});
            x= mean(x,1,'omitnan');
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end

%Hb mean AUC(normal) in 前8个session VS learning speed
for typei=1:2
    h=figure('position',[265,126,900,300]);
    for indi=1:3
        ind=trial.acq(2)+8*(indi-1):trial.acq(2)+8*(indi)-1;
        subplot(1,3,indi)
        for batchi=[1:length(Path)]
            switch typei
                case 1
                    name=['normAUC VS learningspeed_session'];a=Region_AUC_session{batchi};ylab='Δ AUC';lim=[-1 2];
                case 2
                    name=['normFraction VS learningspeed_session'];a=Region_fraction_session{batchi};ylab='Δ Fraction';lim=[0 1];
            end
            y= learning_speed{batchi}(num{batchi});x=a(ind,num{batchi});pre=mean(a(trial.hab(2):trial.hab(3),num{batchi}),1,'omitnan');
            x=(x-repmat(pre(1,:),size(ind,2),1))./(repmat(pre(1,:),size(ind,2),1));
            x= mean(x,1,'omitnan');
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%Hb mean AUC(normal) in  all session VS learning speed
for typei=1
    h=figure('position',[265,126,300,300]);
    for indi=1
        ind=trial.acq(2):trial.acq(3);
        for batchi=[1:length(Path)]
            switch typei
                case 1
                    name=['normAUC VS learningspeed_session_all'];a=Region_AUC_session{batchi};ylab='Δ AUC';lim=[-1 2];
                case 2
                    name=['normFraction VS learningspeed_session_all'];a=Region_fraction_session{batchi};ylab='Δ Fraction';lim=[0 1];
            end
            y= learning_speed{batchi}(num{batchi});x=a(ind,num{batchi});pre=mean(a(trial.hab(2):trial.hab(3),num{batchi}),1,'omitnan');
            x=(x-repmat(pre(1,:),size(ind,2),1))./(repmat(pre(1,:),size(ind,2),1));
            x= mean(x,1,'omitnan');
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
for typei=1
    h=figure('position',[265,126,900,300]);
    for indi=1 :2
        subplot(1,3,indi),
        for batchi=[1,length(Path)]
            for fishi=num{batchi}
                path=Path{batchi}{fishi};
                load(fullfile(path, '/activities_dfdf_align.mat'),'align_win');
                switch indi
                    case 1
                        ind=align_win(:,3);trial.acq(2):max(align_win(:,6));%trial.acq(2):trial.acq(3);
                    case 2
                        ind=max(align_win(:,6)):trial.acq(3);%ind=align_win(:,6);
                    case 3 
                end
                switch typei
                    case 1
                        name=['AUC VS learningspeed_session_all'];a=Region_AUC_session{batchi};ylab='AUC';lim=[0 25];
                    case 2
                        name=['Fraction VS learningspeed_session_all'];a=Region_fraction_session{batchi};ylab='Fraction';lim=[0 1];
                end
                y=mean(learning_perform{batchi}(ind,fishi));
                x=a(ind,fishi);x= mean(x,1,'omitnan');
                scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;   
            end    
        end
        ylabel(['Learning quantity'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
% Hb single AUC(normal) in  all session VS learning speed
for typei=1
    h=figure('position',[265,126,300,300]);name=['normAUCsingleneuron VS learningspeed_session_all'];
    for indi=1
        for batchi=[1:length(Path)]
            a=[];b=[];kk=1;
            for fishi=num{batchi}
                load(fullfile(Path{batchi}{fishi}, '/activities_dfdf_align.mat'),'align_win');
                ind=reshape(align_win,1,[]);%trial.acq(2):trial.acq(3);%reshape(align_win,1,[]);%
                x=mean(Region_AUC_session_all{batchi,fishi}(:,ind),2,'omitnan');id=find(x<15);
                a(kk)=mean(mean(Region_AUC_session_all{batchi,fishi}(id,ind),2,'omitnan'),1,'omitnan');
                b(kk)=learning_quantity{batchi}(fishi);
                kk=kk+1;
            end
            ylab='AUC';lim=[0 25];
            x=a;y=b;
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1])
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
% Hb single AUC(normal) distribution in  all align session
for typei=1
    for batchi=[1:length(Path)]
        h=figure('position',[1,1,900,900]);
        name=['normAUCsingleneuron VS learningspeed_session_all'];
        lim=[0 30];
        a=[];b=[];KK=1;
        for fishi=num{batchi}
            load(fullfile(Path{batchi}{fishi}, '/activities_dfdf_align.mat'),'align_win');
            ind=reshape(align_win,1,[]);%trial.acq(2):trial.acq(3);%reshape(align_win,1,[]);%
            a=mean(Region_AUC_session_all{batchi,fishi}(:,ind),2,'omitnan');
            subplot(3,3,KK),histogram(a,[0:1:25],'Normalization','probability','FaceColor',clr_group(batchi,:));
            ylabel(['Fraction'],'fontsize',16,'FontWeight','bold');
            xlabel('AUC','fontsize',16,'FontWeight','bold');title(num2str(learning_speed{batchi}(fishi)))
            set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
            set(gca,'FontName','Times New Roman');KK=KK+1;
            xlim(lim);ylim([0 0.2]);
        end
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
for typei=1
    h=figure('position',[1,1,900,900]);
    name=['normAUCsingleneuron VS learningspeed_session_all'];
    for batchi=[1,length(Path)]
        lim=[0 30];
        a=[];b=[];KK=1;
        for fishi=num{batchi}
            load(fullfile(Path{batchi}{fishi}, '/activities_dfdf_align.mat'),'align_win');
            ind=reshape(align_win,1,[]);%trial.acq(2):trial.acq(3);%reshape(align_win,1,[]);%
            a=cat(1,a,mean(Region_AUC_session_all{batchi,fishi}(:,ind),2,'omitnan'));
        end
        histogram(a,[0:1:25],'Normalization','probability','FaceColor',clr_group(batchi,:));hold on;
        ylabel(['Fraction'],'fontsize',16,'FontWeight','bold');
        xlabel('AUC','fontsize',16,'FontWeight','bold');title(num2str(learning_speed{batchi}(fishi)))
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');KK=KK+1;
        xlim(lim);ylim([0 0.2]);
        
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
% Hb single AUC(normal) distribution in  all session VS learning speed
for typei=1
    h=figure('position',[265,126,300,300]);
    name=['normAUCsingleneuron VS learningspeed_session_all'];
    for indi=1
        ind=trial.acq(2):trial.acq(3);
        for batchi=[1:length(Path)]
            a=[];b=[];
            for fishi=num{batchi}
                a=cat(1,a,Region_AUC_session_all{batchi,fishi});
                b=cat(1,b,repmat(learning_speed{batchi}(fishi),1,size(Region_AUC_session_all{batchi,fishi},1))');
            end
            
            ylab='Δ AUC';lim=[-1 2];
            x=a(:,ind);
            pre=mean(a(:,trial.hab(2):trial.hab(3)),2,'omitnan');
            x=(x-repmat(pre(:,1),1,size(ind,2)))./(repmat(pre(:,1),1,size(ind,2)));
            x= mean(x,2,'omitnan');y=b;
            scatter(reshape(x,1,[]),reshape(y,1,[]),30,clr_group(batchi,:),'filled');hold on;
        end
        ylabel(['Learning speed'],'fontsize',16,'FontWeight','bold');
        xlabel(ylab,'fontsize',16,'FontWeight','bold');
        set(gca,'XColor','k','YColor','k','linewidth',2);set(gca, 'FontSize',13);
        set(gca,'FontName','Times New Roman');
        xlim(lim);ylim([0 1]);
    end
    saveas(h,string(fullfile(savepath,name)),savetype);
end
%group avg （时间轴是session)
