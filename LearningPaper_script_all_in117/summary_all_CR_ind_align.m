clc;
clear all;
%% load all 
iscutmov=1;fraction_type=1;
savepath='X:\calcium data 20230224\';
load([savepath,'/Path.mat'])
A=struct;
path_all={};kk=1;corr_CS_up_all=cell(4,1);corr_CS_up_fish={};
A=setfield(A,'NUM_CR_UR',[]);A=setfield(A,'fraction_CR_UR',[]);
for groupi=[1:4]
    path=Path{groupi};
    for fishi=1:length(path)
        p=path{fishi};path_all{1,kk}=p;kk=kk+1;
        load(fullfile(p,'CR_ind_summary_align.mat'),'ind_all_CSUS_RESPONSIVE','Fraction_in_region_type');
        load(fullfile(p,'CR_ind_summary_align.mat') ,'corr_CS_up','corr_CS_down');
        b=cat(1,corr_CS_up_all{groupi},corr_CS_up);
        corr_CS_up_all{groupi}=b;corr_CS_up_fish{groupi}{fishi}=corr_CS_up;
        NUM_CR_UR=[];fraction_CR_UR=[];
        for typei=1:10
            for sessioni=1:min(size(corr_CS_up,2),size(Fraction_in_region_type{typei,iscutmov,fraction_type},2))
                NUM_CR_UR(sessioni,typei)=length(ind_all_CSUS_RESPONSIVE{typei,sessioni,iscutmov})./size(corr_CS_up,1);
                fraction_CR_UR(sessioni,:,typei)=Fraction_in_region_type{typei,iscutmov,fraction_type}(:,sessioni);
            end
        end
        a=getfield(A,'NUM_CR_UR');a=cat(3,a,NUM_CR_UR);A=setfield(A,'NUM_CR_UR',a);
        a=getfield(A,'fraction_CR_UR');a=cat(4,a,fraction_CR_UR);A=setfield(A,'fraction_CR_UR',a);
    end 
end
A=setfield(A,'corr_CS_up',corr_CS_up_all);A=setfield(A,'corr_CS_up_fish',corr_CS_up_fish);

save([savepath,'/Summary_all_CRUR_align.mat'],'A','-v7.3');
%% select fish
numL=[1,3,4,6,7,8,9,10];
numC=[12,13,14,16,17,18,19,20];numNL=[21:29];numFL=[30:35];
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
group_label={'L','FL','NL','C'};
labels_typei_ADJ={'CS_avtivation','CS_inhibition','US_activation','CSupregulate','CSdownregulate','CSstableregulate','USupregulate','USdownregulate','USstableregulate','CSnewemergedneuron'};
load(fullfile(path{1},'CR_ind_summary_align.mat'),'Label_region');
regionxx=categorical(Label_region{1,1}(:,1));regionxx = reordercats(regionxx,Label_region{1,1}(:,1));
mean1= @(x)(mean(x,ndims(x),'omitnan'));
error1= @(x)(std(x,[],ndims(x),'omitnan')./sqrt(length(x)));
error2= @(x)(mean(x,ndims(x),'omitnan')-std(x,[],ndims(x),'omitnan'));
clr_group=[colorplus(389);colorplus(253);colorplus(233);colorplus(448)];    
clr_cmap = addcolorplus(336);%cmap = GetColormap('jet',81);
clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
clr_CStypes = addcolorplus(275);%clr_CStypes= ColorMap(clr_CStypes,3);clr_CStypes = flipud(clr_CStypes);
clr_UStypes = addcolorplus(272);%clr_UStypes= ColorMap(clr_UStypes,3);clr_UStypes = flipud(clr_UStypes);
%% plot CS corr
savepath_eps='X:\calcium data 20230224\FIGURE\';
savetype='jpg';
for groupi=1:4;
    stAvrcorr_all=getfield(A,'corr_CS_up');
    stAvrcorr_all=stAvrcorr_all{groupi};
    name=[group_label{groupi},'_hist_of_StimAveCorr_autothr','_ALL'];
    h=figure('name',name);
%     hh=[];
%     for sessioni=1:size(stAvrcorr_all,2)
%         hh(sessioni,1)=histogram( stAvrcorr_all(:,sessioni),'BinEdges',-1:0.03:1,'FaceColor',clr_session(sessioni,:),'Normalization','probability');hold on
%     end
%     xlabel('CS corr.');ylabel('Percentage');set(gca,'fontsize',16);legend([hh],sessionx);hold off;
    stAvrcorr={};
    for sessioni=1:size(stAvrcorr_all,2)
        stAvrcorr{sessioni}=stAvrcorr_all(:,sessioni);
    end
    JP=joyPlot(stAvrcorr,'ColorMode','Order','Scatter','on','ColorList',clr_session ,'MedLine','on','Kdensity','Sep',30);
    JP=JP.draw();
    legendHdl=JP.getLegendHdl();
    for sessioni=1:size(stAvrcorr_all,2)
        JP.setMedLine(sessioni,'Color',clr_session(sessioni,:))
    end
    %legend(sessionx)
    saveas(h,[savepath_eps,name],savetype);
    name=[group_label{groupi},'_histcount_of_StimAveCorr_autothr','_ALL_errorbar'];
    h=figure('name',name);
    x=[];y=[];
    for fishi=1:size(corr_CS_up_fish{groupi},2)
        kk=1;
        for sessioni=1:size(corr_CS_up_fish{groupi}{fishi},2)
            X=corr_CS_up_fish{groupi}{fishi}(:,sessioni);
            hh=histogram( X,'BinEdges',-1:0.01:1,'Normalization','cdf','visible','off');
            y(sessioni,:,fishi)=hh.Values;
            x(sessioni,:,fishi)=hh.BinEdges(2:end);
            kk=kk+1;
        end
    end

    xx=mean(x,3,'omitnan');yy=mean(y,3,'omitnan');errBar=std(y,0,3,'omitnan')/size(y,3);kk=1;
    for sessioni=1:size(corr_CS_up_fish{groupi}{fishi},2)
        shadedErrorBar(xx(sessioni,:),yy(sessioni,:),errBar(sessioni,:),'lineProps',clr_session( sessioni,:));hold on;
        hh(kk)=plot(mean(xx(sessioni,:),1,'omitnan'),mean(yy(sessioni,:),1,'omitnan'),'linewidth',1.5,'color',clr_session( sessioni,:),'linestyle','-');hold on
        kk=kk+1;
    end
    legend([hh],sessionx,'Location','southeast');xlabel('CS corr.');ylabel('Cdf.');set(gca,'fontsize',16);xlim([-1 1]);ylim([0 1]);
    saveas(h,[savepath_eps,name],savetype);
    name=[group_label{groupi},'_boxplot_of_StimAveCorr_autothr','_ALL'];
    h=figure('name',name);
    boxplot(stAvrcorr_all,'Labels',string(sessionx),'Symbol','wo','OutlierSize',2,'Orientation','horizontal');
    xlabel('CS corr.');set(gca,'fontsize',16); xlim([-1 1]);
    saveas(h,[savepath_eps,name],savetype);
end
%% plot CS activation
data=getfield(A,'NUM_CR_UR');
typei=[1:3];
data=squeeze(data(:,typei,:));
name=[group_label{groupi},'_Num_of',labels_typei_ADJ{typei(1)}];
h=figure('position',[33,129,1721,403],'name',name,'color','w');
tiledlayout(1,4);
for kk=1:4
    switch kk
        case 1
            ind=numL;
        case 4
            ind=numC;
        case 3
            ind=numNL;
        case 2
            ind=numFL;
    end
    dataL=squeeze(data(:,:,ind));
    ax1 = nexttile;
    b=bar(squeeze(mean1(dataL)),'FaceColor','flat','FaceAlpha',1,'BarWidth',1);hold on;
    for ii=1:length(typei)
        b(ii).CData = clr_group(typei(ii),:);
        errorbar(b(ii).XEndPoints,squeeze(mean1(dataL(:,ii,:))),squeeze(error1(dataL(:,ii,:))),'LineStyle','none','Color','k', 'LineWidth',.8);hold on;
    end
    set(gca,'fontsize',16,'XTickLabel',sessionx,'XTickLabelRotation',45);title(group_label{kk});
    ylim([0 0.2]);if kk==4 legend(b,{labels_typei_ADJ{typei}},'Interpreter','none') ;end
end
saveas(h,[savepath_eps,name],savetype);

data=getfield(A,'fraction_CR_UR');
for typei=1:10;
    for kk=1:4
        switch kk
            case 1
                ind=numL;
            case 4
                ind=numC;
            case 3
                ind=numNL;
            case 2
                ind=numFL;
        end
        dataL=mean1(data(:,:,typei,ind))*100;
        cdata=discretize(dataL,[1:1:100]);
        name=[group_label{kk},'_Fraction_of',labels_typei_ADJ{typei}];
        h=bubbleheatmap(dataL',dataL',sessionx,regionxx,clr_cmap,[1 80]);
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end

Label={'CS-upregulate','CS-downregulate','CS-stableregulate'};
typei=4:6;
data=getfield(A,'fraction_CR_UR');
for groupi=1:4
    switch groupi
        case 1
            ind=numL;
        case 4
            ind=numC;
        case 3
            ind=numNL;
        case 2
            ind=numFL;
    end
    dataL=mean1(data(:,:,typei,ind))*100;
    figure;kk=1;
    for sessioni=2:size(dataL,1)
        for regioni=1:size(dataL,2)
            subplot(size(dataL,1),size(dataL,2),kk),
            mypieplot(squeeze(dataL(sessioni,regioni,:)),Label,clr_CStypes);kk=kk+1;
            title(regionxx(regioni));
        end
    end
end