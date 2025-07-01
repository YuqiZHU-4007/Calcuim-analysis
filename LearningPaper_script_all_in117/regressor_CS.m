set(0,'defaultfigurecolor','w')
clc;clear all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
US_related_stimcorr={};
CS_related_stimcorr={};
CSUS_related_stimcorr={};
Thr_related={};
for batchi=1:4
    path=Path{batchi};
    for fishi=1:length(path)
        %path='X:\calcium data 20230224\\20220709\fish1\'
        load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align');
        load(fullfile(path{fishi}, '/para.mat'));
        nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
        A=activities_dfdf_align;
        %thr
        frame_indCS=frame.cs_start:frame.cs_end;
        frame_indUS=frame.us_start:frame.us_start+5;
        a=squeeze(nanmean(A(frame_indUS,:,:),1));
        b= squeeze(nanmean(A(frame_indCS,:,:),1));
        B1=(a-b)>0 & abs(a-b)>b*0.3;B2=(a-b)<0 & abs(a-b)>b*0.1;
        Thr_related{batchi,fishi,1}=B1;Thr_related{batchi,fishi,2}=B2;
        %ind_thr=find(sum(B1(2:4,:),1)>=2 & sum(B2(5:7,:),1)>=2);
        %regress
        frame_indCS=1:frame.cs_end;
        frame_indUS=frame.us_start:frame.per_cycle;
        sessions=size(A,2);
        regressor1=[];regressor2=[];regressor3=[];regress1=zeros(5,sessions*frame.per_cycle); regress2=zeros(5,(sessions-2)*frame.per_cycle);regress3=zeros(5,(sessions)*frame.per_cycle);
        figure,
        for ii=1:sessions
            regressor1(:,ii)=zeros(1,frame.per_cycle);
            regressor1(frame.cs_start:frame.cs_end,ii)=(ii-1)/sessions*2;
            subplot(1,sessions,ii),plot(regressor1(:,ii));
            ylim([0 2])
        end
        figure,
        for ii=1:sessions-2
            regressor2(:,ii)=zeros(1,frame.per_cycle);
            regressor2(frame.us_start:frame.us_start+2,ii)=(sessions-2-ii)/(sessions-2)*2;
            subplot(1,sessions,ii),plot(regressor2(:,ii));
            ylim([0 2])
        end
        figure,
        for ii=1:sessions
            regressor3(:,ii)=zeros(1,frame.per_cycle);
            regressor3(frame.cs_start:frame.cs_end,ii)=(ii-1)/sessions;
            if ii>=2 && ii<=sessions-1
            regressor3(frame.us_start:frame.us_start+2,ii)=(sessions-ii+1)/(sessions)*2;
            end
            subplot(1,sessions,ii),plot(regressor3(:,ii));
            ylim([0 2])
        end
        regress1(1,:)=reshape(regressor1,1,[]);
        regress2(1,:)=reshape(regressor2,1,[]);
        regress3(1,:)=reshape(regressor3,1,[]);
        close all;
        
        [regressors,~,~,~] = GetMotorRegressor(regress2);output=[regressors(1,1).im]./max(regressors(1,1).im);
        output= reshape(output,frame.per_cycle,sessions-2);output= output(frame_indUS,:);output= reshape(output,[],1);
        AA=reshape(A(frame_indUS,2:7,:),length(frame_indUS)*(sessions-2),[]);
        a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
        [stimcorr,~] = MotorSourceCorrelation(AA', output',[]);
        US_related_stimcorr{batchi,fishi}=stimcorr;
        
        [regressors,~,~,~] = GetMotorRegressor(regress1);output=[regressors(1,1).im]./max(regressors(1,1).im);
        output= reshape(output,frame.per_cycle,sessions);output= output(frame_indCS,:);output= reshape(output,[],1);
        %figure,plot(output)
        AA=reshape(A(frame_indCS,:,:),length(frame_indCS)*(sessions),[]);
         a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
        [stimcorr,~] = MotorSourceCorrelation(AA',output',[]);
        CS_related_stimcorr{batchi,fishi}=stimcorr;
        
        [regressors,~,~,~] = GetMotorRegressor(regress3);output=[regressors(1,1).im]./max(regressors(1,1).im);
        output= reshape(output,frame.per_cycle,sessions);output= output(:,2:sessions-1);output= reshape(output,[],1)';
        %figure,plot(output)
        AA=reshape(A(:,2:7,:),frame.per_cycle*(sessions-2),[]);
         a=find(isnan(AA(:,1)));AA(a,:)=0;output(a)=0;
        [stimcorr,~] = MotorSourceCorrelation(AA',output,[]);
        CSUS_related_stimcorr{batchi,fishi}=stimcorr;
    end
end
save([savepath '\regressorCSUSshift.mat'],'US_related_stimcorr','CS_related_stimcorr','CSUS_related_stimcorr','Thr_related');

%%
clc;clear all
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
load([savepath '\regressorCSUSshift.mat']);
group={'Learner','Unpair','Non-learner','Faded-learner'};
columns = {'t','x', 'y', 'z'};
coloredge=[colorplus(389);colorplus(448);colorplus(233);colorplus(253)];
colorface=[[1,1,1];[1,1,1];[1,1,1];[1,1,1]];
colorerrorbar=coloredge;
sap='X:\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';

figure('position',[5,5,900,800]),
for batchi=1:4
    path=Path{batchi};
     subplot(2,2,batchi)
    for fishi=num{batchi}
        stimcorr=CS_related_stimcorr{batchi,fishi};
        histogram(stimcorr,[-1:0.04:1],'Normalization', 'probability');hold on;
        title(group{batchi})
    end
end
figure('position',[5,5,900,800]),
for batchi=1:4
    path=Path{batchi};
    subplot(2,2,batchi)
    for fishi=num{batchi}
        stimcorr=US_related_stimcorr{batchi,fishi};
       histogram(stimcorr,[-1:0.04:1],'Normalization', 'probability');hold on;
    end
end
%trace for each fish
col=1;
for typei=4
    for batchi=1
        path=Path{batchi};
        figure('position',[5,5,1800,800]),
        for fishi=1:length(Path{batchi})
            %           pp=checkpath(fullfile(path{fishi},'singletrial'));
            %         load(fullfile(pp,'CR_ind_summary_singletrial.mat'), 'ind_all_CSUS_RESPONSIVE');
            load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align');
            load(fullfile(path{fishi}, '/para.mat'));
            nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
            A=reshape(activities_dfdf_align,size(activities_dfdf_align,1)*size(activities_dfdf_align,2),[]);
            ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei);
            CS_related_act=A(:,ind_type);
            if col==1
                a=1.5;
                subplot(1,10,fishi)
                B=reshape(CS_related_act,[],size(activities_dfdf_align,2),length(ind_type));
                patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[0,0,20,20],'r','facealpha',0.3,'edgealpha',0);hold on;
                for jj=size(activities_dfdf_align,2):-1:1
                    b=squeeze(B(:,jj,:))+a*(size(activities_dfdf_align,2)+0.5-jj);
                    x=1:size(B,1);
                    shadedErrorBar(x,nanmean(b,2),nanstd(squeeze(B(:,jj,:)),[],2),'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on;%./sqrt(size(b,2)-1)
                    plot(nanmean(b,2),'k','linewidth',2);hold on
                    title(nn);ylim([0 14]);
                end
            else
                subplot(10,1,fishi)
                patch([[frame.cs_start:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
                    [frame.cs_end:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
                    [frame.cs_end:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]',...
                    [frame.cs_start:frame.per_cycle:size(activities_dfdf_align,2)*frame.per_cycle]']',...
                    [0*ones(size(activities_dfdf_align,2),1), 0*ones(size(activities_dfdf_align,2),1) 1*ones(size(activities_dfdf_align,2),1) 1*ones(size(activities_dfdf_align,2),1)]',...
                    'r','facealpha',0.5,'edgealpha',0);hold on;
                shadedErrorBar(1:size(CS_related_act,1),nanmean(CS_related_act,2),nanstd(CS_related_act,[],2)./sqrt(length(ind_type)-1),'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                plot(mean(CS_related_act,2),'k','linewidth',2);%legend('regressor','corr > thr');
                ylim([-1,1])
                hold on;
                title(nn)
            end
        end
    end
end
%% typical case
for typei=4
    for batchi=1
        path=Path{batchi};
        for fishi=4
            load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align');
            load(fullfile(path{fishi}, '/para.mat'));
            nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
            A=reshape(activities_dfdf_align,size(activities_dfdf_align,1)*size(activities_dfdf_align,2),[]);
            ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei);
            CS_related_act=A(:,ind_type);
            CS_related_act_r=reshape(CS_related_act,75,8,[]);
            for ss=1:8
                h=figure('position',[5,5,400,400]);
                a=squeeze(CS_related_act_r(:,ss,:));
                patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[-0.2, -0.2,0.6,0.6], 'r','facealpha',0.3,'edgealpha',0);hold on;
                shadedErrorBar(1:size(a,1),nanmean(a,2),nanstd(a,[],2)./sqrt(length(ind_type)-1)*2,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                plot(mean(a,2),'k','linewidth',2);%legend('regressor','corr > thr');
                ylim([-0.3,0.6])
                hold on;
                title(num2str(ss));
                h.InvertHardcopy = 'off';
                savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\CSUSshift_regressors');
                name=['avg_',num2str(ss)];
                saveas(h,[savepath_eps,'\',name],'png');
                savefig(h,[savepath_eps,'\',name]);
            end
        end
    end
end
for ss=1:8
    h=figure('position',[5,5,400,400]);
    a=squeeze(regressor3(:,ss));
    patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[-0.2, -0.2,2,2], 'r','facealpha',0.3,'edgealpha',0);hold on;
    plot(mean(a,2),'k','linewidth',2);%legend('regressor','corr > thr');
    ylim([-0.3,2])
    hold on;
    title(num2str(ss));
    h.InvertHardcopy = 'off';
    savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\CSUSshift_regressors');
    name=['Regressors_',num2str(ss)];
    saveas(h,[savepath_eps,'\',name],'png');
    savefig(h,[savepath_eps,'\',name]);
end
[regressors,~,~,~] = GetMotorRegressor(regress3);output=[regressors(1,1).im]./max(regressors(1,1).im);
output= reshape(output,frame.per_cycle,sessions);%output= output(:,2:sessions-1);output= reshape(output,[],1)';
for ss=1:8
    h=figure('position',[5,5,400,400]);
    a=squeeze(output(:,ss));
    patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[-0.2, -0.2,2,2], 'r','facealpha',0.3,'edgealpha',0);hold on;
    plot(mean(a,2),'k','linewidth',2);%legend('regressor','corr > thr');
    ylim([-0.3,2])
    hold on;
    title(num2str(ss));
    h.InvertHardcopy = 'off';
    savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign\CSUSshift_regressors');
    name=['Regressors_conv',num2str(ss)];
    saveas(h,[savepath_eps,'\',name],'png');
    savefig(h,[savepath_eps,'\',name]);
end
%% mapback
figure('position',[5,5,900,800]),
for typei=4
    for batchi=1:4
        path=Path{batchi};
        subplot(2,2,batchi);
        supervolxeli=readmatrix([warped_SyN_csv_path ,'20210709fish2','vol_env_spatialloc_warped_SyN.csv']);
        scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),0.5,[0.5 0.5 0.5]);hold on;
        for fishi=num{batchi}()
            warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
            nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
            supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
            
            ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei);
            scatter3(supervolxeli(ind_type,1),supervolxeli(ind_type,2),supervolxeli(ind_type,3),16,'r','filled');grid off;view([0,-90]);
            set(gca,'visible','off')
        end
        title(group{batchi})
    end
end
for typei=4
    for batchi=1:4
        path=Path{batchi};
        for ii=1:2
            figure('position',[5,5,120,1200]),
            for fishi=1:length(Path{batchi})
                warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
                nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
                supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                subplot(10,1,fishi);
            ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei);
                scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),0.5,[0.5 0.5 0.5]);hold on;
                scatter3(supervolxeli( ind_type,1),supervolxeli( ind_type,2),supervolxeli( ind_type,3),16,'r','filled');grid off;
                if ii==1 view([0,-90]);else view([-90,0]);end
                title([num2str(length(ind_type))]);%set(gca,'visible','off')
            end
        end
    end
end
%%
Brain_region_id={};
Index_in_region={};
Fraction_in_region_type={};
Loc_in_region_cell={};
Num_in_region={};
Brain_region_id_shuffled={};
Index_in_region_shuffled={};
Fraction_in_region_type_shuffled={};
Loc_in_region_cell_shuffled={};
Num_in_region_shuffled={};
for batchi=1:4
    path=Path{batchi};
    for fishi=num{batchi}
        warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
        nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        
        load(fullfile(path{fishi},'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
        for typei=1:4
            ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,typei);
            ind_shuffled=randi(size(supervolxeli,1),1,length(ind_type));
            %Brain_region_id{batchi}{fishi,ii}=brain_region_id(ind_type,:);
            for regioni=1:length(Label)
                ind_inregion=intersect(ind_type,index_in_region_in_clust{regioni,1});
                Index_in_region{batchi}{fishi,typei}{regioni}=ind_inregion;
                Fraction_in_region_type{batchi}{fishi,typei}(regioni)=length(ind_inregion)./length(index_in_region_in_clust{regioni,1});
                Loc_in_region_cell{batchi}{fishi,typei}{regioni}=supervolxeli(ind_inregion,:);
                Num_in_region{batchi}{fishi,typei}(regioni)= length(ind_inregion)./length(ind_type);
                 ind_inregion_shuffled=intersect(ind_shuffled,index_in_region_in_clust{regioni,1});
                Index_in_region_shuffled{batchi}{fishi,typei}{regioni}=ind_inregion_shuffled;
                Fraction_in_region_type_shuffled{batchi}{fishi,typei}(regioni)=length(ind_inregion_shuffled)./length(index_in_region_in_clust{regioni,1});
                Loc_in_region_cell_shuffled{batchi}{fishi,typei}{regioni}=supervolxeli(ind_inregion_shuffled,:);
                Num_in_region_shuffled{batchi}{fishi,typei}(regioni)= length(ind_inregion_shuffled)./length(ind_shuffled);
                %Num_in_region{batchi}{fishi,ii}(regioni)= length(ind_inregion);
                %Num_in_region{batchi}{fishi,ii}(regioni)= length(ind_type);
                %Label_region{fishi}(regioni)=Label{regioni};
            end
        end
    end
end
lim=[0.6,0.15,0.5 0.05];
yb=[0.05,0.15];
for typei=4
    figure('position',[5,5,1200,1200]);
    for batchi=1:4;
        subplot(4,1,batchi),
        barData1=zeros(30,length(num{batchi}))*NaN;
        a=[];b=[];k=1;
        for j=num{batchi}
            a(:,k)=Fraction_in_region_type{batchi}{j,typei};
            b(:,k)=Fraction_in_region_type_shuffled{batchi}{j,typei};
            k=k+1;
        end
        barData1=a'; barData2=b';
        barHdl=[];
        for i=1:size(barData1,2)
            barHdl(i)=bar(i,nanmean(barData1(:,i)),0.75,'FaceColor',colorface(batchi,:),'LineWidth',2,'EdgeColor',coloredge(batchi,:));hold on
        end
        % 设置0基线属性
        for i=1:size(barData1,2)% 绘制抖动散点图
            scatter(i+rand(1,length(find(~isnan(barData1(:,i))))).*.6-.3,barData1(find(~isnan(barData1(:,i))),i),20,coloredge(batchi,:),"filled");hold on
        end
        for i=1:size(barData1,2)% 绘制误差条
            errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
                'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(batchi,:),'CapSize',8);hold on
        end
        hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=12;ax.XTick=1:30;
        ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
        ax.YTick=0:0.05:1;ax.YLim=[0,lim(typei)];% ax.Title.String='Pericyte response ratio';
        ylabel('Fraction','FontName','Arial','Fontsize',16);box off
       % breakyaxis(gca,yb);
        xtickangle(30);
        set(ax(1),'XTickLabel',Label);box off;% set(gca,'Position', [.15 .18 .75 .75]);
    end
end

lim=[0.3,0.35,0.4 0.4];
for typei=4
    figure('position',[5,5,1200,1200]);
    for batchi=1:4;
        subplot(4,1,batchi),
        barData1=zeros(30,length(num{batchi}))*NaN;
        a=[];k=1;
        for j=num{batchi}
            a(:,k)=Num_in_region{batchi}{j,typei};% Fraction_in_region_type{batchi}{j,typei};
            k=k+1;
        end
        barData1=a';
        barHdl=[];
        for i=1:size(barData1,2)
            barHdl(i)=bar(i,nanmean(barData1(:,i)),0.75,'FaceColor',colorface(batchi,:),'LineWidth',2,'EdgeColor',coloredge(batchi,:));hold on
        end
        % 设置0基线属性
        for i=1:size(barData1,2)% 绘制抖动散点图,老版本方案
            scatter(i+rand(1,length(find(~isnan(barData1(:,i))))).*.6-.3,barData1(find(~isnan(barData1(:,i))),i),20,coloredge(batchi,:),"filled");hold on
        end
        for i=1:size(barData1,2)% 绘制误差条
            errorbar(i,nanmean(barData1(:,i)),nanstd(barData1(:,i),1)./sqrt(length(find(~isnan(barData1(:,i))))),...
                'vertical','LineStyle','none','LineWidth',2,'Color',colorerrorbar(batchi,:),'CapSize',8);hold on
        end
        hold on;ax=gca,ax.Color='none';ax.LineWidth=1.5;ax.TickDir='in';ax.FontSize=12;ax.XTick=1:30;
        ax.Title.FontWeight='bold';ax.Title.FontSize=18;ax.YDir='normal';
        ax.YTick=0:0.1:1;ax.YLim=[0,lim(typei)];% ax.Title.String='Pericyte response ratio';
        ylabel('Fraction','FontName','Arial','Fontsize',16);box off
        xtickangle(30);
        set(ax(1),'XTickLabel',Label);box off;% set(gca,'Position', [.15 .18 .75 .75]);
    end
end

%
typei =4;batchi=3;
a=[];b=[];k=1;
for j=num{batchi}
    a(:,k)=Fraction_in_region_type{batchi}{j,typei};
    b(:,k)=Fraction_in_region_type_shuffled{batchi}{j,typei};
    k=k+1;
end
%% cb,hb trace
load(fullfile(Path{1}{1},'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
region_label={{'L H'}};%{{'L P','R P'},{'L H'},{'R H'},{'L TeO','R TeO'},{'Va','CCe'}};
region_label1=cell(1,length(region_label));
for ii=1:length(region_label)
    for jj=1:length(region_label{ii})
        region_label1{ii}(jj)=find(strcmp(region_label{ii}(jj),Label));
    end
end
for align_type=1
    for regioni=1:length(region_label1)
        for batchi=1
            path=Path{batchi};
            figure('position',[5,5,1800,800]),
            for fishi=1:length(Path{batchi})
                switch align_type
                    case 1
                        a=load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align');
                        activities_dfdf_align=a.activities_dfdf_align;
                    case 2
                        a=load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align_go');
                        activities_dfdf_align=a.activities_dfdf_align_go;
                        
                    case 3
                        a=load(fullfile(path{fishi}, '/activities_dfdf_align.mat'), 'activities_dfdf_align_nogo');
                        activities_dfdf_align=a.activities_dfdf_align_nogo;       
                end
                load(fullfile(path{fishi},'\brain_region_related_statistic.mat'),'index_in_region_in_clust','brain_region_id','Label');
                nn=[path{fishi}(end-14:end-7),path{fishi}(end-5:end-1)];
                A=reshape(activities_dfdf_align,size(activities_dfdf_align,1)*size(activities_dfdf_align,2),[]);
                ind_type=find_type_regressorCS(CS_related_stimcorr,US_related_stimcorr,CSUS_related_stimcorr,Thr_related,batchi,fishi,4);
                ind_inregion=[];
                for jj=1:length(region_label1{regioni})
                    a= intersect(ind_type,index_in_region_in_clust{region_label1{regioni}(jj),1});
                    ind_inregion=cat(2,ind_inregion,a);
                end
                B=A(:,ind_inregion);
                B=reshape(B,[],size(activities_dfdf_align,2),length(ind_inregion));
                
                subplot(1,10,fishi),patch([frame.cs_start,frame.cs_end,frame.cs_end,frame.cs_start],[0,0,20,20],'r','facealpha',0.3,'edgealpha',0);hold on;
                a=2.5;
                for jj=size(activities_dfdf_align,2):-1:1
                    b=squeeze(B(:,jj,:))+a*(size(activities_dfdf_align,2)+0.5-jj);
                    x=1:size(B,1);
                    if ~isempty(b)
                        shadedErrorBar(x,nanmean(b,2),nanstd(b,[],2)./sqrt(length(ind_inregion)-1),'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                        plot(nanmean(b,2),'k','linewidth',2);hold on
                    end
                end
                    title([nn,region_label{regioni}]);ylim([0 20]);
            end
        end
    end
end
%% spatial location
a={};loc={};columns = {'t','x', 'y', 'z'};

for batchi=1:4
    path=Path{batchi};
    %figure,
    loc{batchi}=[nan,nan,nan,nan,nan];typei=4;
    for fishi=num{batchi}
        for regioni=1:30
            if ~isempty(Loc_in_region_cell{batchi}{fishi,typei}{regioni})
                loc{batchi}=cat(1,loc{batchi},Loc_in_region_cell{batchi}{fishi,typei}{regioni});
            end
            a{batchi}(fishi,regioni)= Fraction_in_region_type{batchi}{fishi,typei}(regioni);
            % scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),2,[0.5 0.5 0.5]);hold on;
            %scatter3(loc(:,1),loc(:,2),loc(:,3),16,'r','filled');grid off;view([0,-90]);hold on;
        end
    end
    t=ones(size(loc{batchi},1),1);
    M=table(t,loc{batchi}(:,1),loc{batchi}(:,2),loc{batchi}(:,3),'VariableNames', columns);
    name=strcat(group(batchi),'_',num2str(typei),'.txt');
    outpath=fullfile(sap,name);
    writetable(M,string(outpath),'Delimiter','space','WriteVariableNames',false);
end



