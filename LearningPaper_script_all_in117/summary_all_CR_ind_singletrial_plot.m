
%% demo
num=9;res=[0.66,0.66,10];v=[90,90];radium=floor(15/0.66);
load(fullfile(path_all{num},'para.mat'));load(fullfile(path_all{num},'env.mat'));
load(fullfile(path_all{num}, 'activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo')
load(fullfile(path_all{num},'singletrial','CR_ind_summary_singletrial.mat'),'CR_ind_up');
supervoxel_raw=env.supervoxel;
for ii=unique(supervoxel_raw(:,3))'
    id=find(supervoxel_raw(:,3)==ii);
    supervoxel_raw(id,3)=(ii-1)*10/0.66+1;
end
supervoxel_raw(:,1)=env.width-supervoxel_raw(:,1);
for ii=1:size(align_win,2)
    a=align_win(:,ii);a(isnan(a))=[];
    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((length(a)/3)),3);t.TileSpacing = 'compact';t.Padding = 'compact';
    for jj=1:length(a)
        ind=find(CR_ind_up(:,a(jj),iscutmov)==1);leg=[num2str(a(jj))];
        nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
        set(gca,'fontsize',16)
    end
end



%%
clc;clear all
addpath(genpath('X:\calcium data 20230224\script'));
savepath='X:\calcium data 20230224\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load([savepath,'/Path.mat'])
%% select fish
numL=num{1};
numC=num{2}+10;numNL=num{3}+10+10;numFL=num{4}+10+10+9;
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
group_label={'Learner (n=8)','Faded-Learner (n=6)','Non-Learner (n=9)','Control (n=8)'};
align_label={'All trials','Acting trials','Non-acting trials'};
load(fullfile(Path{1}{1},'CR_ind_summary_align.mat'),'Label_region');
load(fullfile(Path{1}{1},'para.mat'));
regionxx=categorical(Label_region{1,1}(:,1));regionxx = reordercats(regionxx,Label_region{1,1}(:,1));
mean1= @(x)(mean(x,ndims(x),'omitnan'));
error1= @(x)(std(x,[],ndims(x),'omitnan'));
error2= @(x)(std(x,[],ndims(x),'omitnan')./sqrt(length(x)));
error2= @(x)(mean(x,ndims(x),'omitnan')-std(x,[],ndims(x),'omitnan'));
clr_group=[colorplus(389);colorplus(253);colorplus(233);colorplus(448)];
clr_cmap =addcolorplus(300);addcolorplus(336);%cmap = GetColormap('jet',81);
clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
clr_CStypes = addcolorplus(275);%clr_CStypes= ColorMap(clr_CStypes,3);clr_CStypes = flipud(clr_CStypes);
clr_UStypes = addcolorplus(272);%clr_UStypes= ColorMap(clr_UStypes,3);clr_UStypes = flipud(clr_UStypes);
region_id=setdiff(1:length(regionxx),[18:20,27,28]);
%region_id=1:30;
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
%% loadpath
res=[0.66,0.66,10];
savepath_eps=checkpath('X:\calcium data 20230224\Summary_all_CRUR_singletrial_3sd_newalign');
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial.mat',...
    'NUM_CR_UR_align_all','fraction_CR_UR_align_all','region_CS_response_amp_align_all',...
    'Response_amp_region_align_all');
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_append.mat','Response_amp_regddion_align_all');
%load('X:\calcium data 20230224\Summary_all_CRUR_singletrial_append.mat','Response_amp_all_session_align_all');
load('X:\calcium data 20230224\Summary_all_CRUR_singletrial.mat','Response_amp_align_all');
%load('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial_2sd.mat','supervoltemp_type_align_all');
savetype='jpg';
temp_env=load('X:\calcium data 20230224\�����ָ�\env.mat');
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);


%% p-value
Pval_pre_region={};Reject_pre_region={};
data=fraction_CR_UR_align_all;%fraction_CR_UR_align_all(session_align,:,typei,align_typei,fishi)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx);clim=[1 80];
        end
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
            data1=squeeze(data(ss,:,typei,align_typei,ind));%data1(setdiff(1:size(data,1),ss),:,:)=nan;
            for regioni=1:size(data1,2)
                for sessioni=1:length(ss)
                    a=squeeze(data1(1,regioni,:));b=squeeze(data1(sessioni,regioni,:));
                    if ~(sum(isnan(a))==size(data1,3) | sum(isnan(b))==size(data1,3))
                        [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                        Pval_pre_region{typei}(ss(sessioni),regioni,align_typei,groupi)=p;
                        Reject_pre_region{typei}(ss(sessioni),regioni,align_typei,groupi)=p;
                    end
                end
            end
        end
    end
end
Pval_group_region={};Reject_group_region={};
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:34;
  if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        for groupi=1:4
            switch groupi
                case 1
                    indi=numL;
                case 4
                    indi=numC;
                case 3
                    indi=numNL;
                case 2
                    indi=numFL;
            end
            for groupj=1:4
                switch groupj
                    case 1
                        indj=numL;
                    case 4
                        indj=numC;
                    case 3
                        indj=numNL;
                    case 2
                        indj=numFL;
                end
                data1=squeeze(data(ss,:,typei,align_typei,indi));
                data2=squeeze(data(ss,:,typei,align_typei,indj));
                for regioni=1:size(data1,2)
                    for sessioni=1:length(ss)
                        a=squeeze(data1(sessioni,regioni,:));b=squeeze(data2(sessioni,regioni,:));
                        if ~(sum(isnan(a))==size(data1,3) | sum(isnan(b))==size(data2,3))
                            [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                            Pval_group_region{typei,groupi}(ss(sessioni),regioni,align_typei,groupj)=p;
                            Reject_group_region{typei,groupi}(ss(sessioni),regioni,align_typei,groupj)=p;
                        end
                    end
                end
            end
        end
    end
end

Pval_pre={};Reject_pre={};Pval_group={};Reject_group={};
data=NUM_CR_UR_align_all;%NUM_CR_UR_align(session_align,typei,align_typei,fishi)
for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
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
            data1=squeeze(data(ss,typei,align_typei,ind));%data1(setdiff(1:size(data,1),ss),:,:)=nan;
            for sessioni=1:length(ss)
                a=squeeze(data1(1,:));b=squeeze(data1(sessioni,:));
                if ~(sum(isnan(a))==size(data1,2) | sum(isnan(b))==size(data1,2))
                    [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                    Pval_pre{typei}(ss(sessioni),align_typei,groupi)=p;
                    Reject_pre{typei}(ss(sessioni),align_typei,groupi)=p;
                end
            end
        end
    end
end
for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        for groupi=1:4
            switch groupi
                case 1
                    indi=numL;
                case 4
                    indi=numC;
                case 3
                    indi=numNL;
                case 2
                    indi=numFL;
            end
            for groupj=1:4
                switch groupj
                    case 1
                        indj=numL;
                    case 4
                        indj=numC;
                    case 3
                        indj=numNL;
                    case 2
                        indj=numFL;
                end
                data1=squeeze(data(ss,typei,align_typei,indi));
                data2=squeeze(data(ss,typei,align_typei,indj));
                for sessioni=1:length(ss)
                    a=squeeze(data1(sessioni,:));b=squeeze(data2(sessioni,:));
                    if ~(sum(isnan(a))==size(data1,2) | sum(isnan(b))==size(data2,2))
                        [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                        Pval_group{typei,groupi}(ss(sessioni),align_typei,groupj)=p;
                        Reject_group{typei,groupi}(ss(sessioni),align_typei,groupj)=p;
                    end
                end
            end
        end
    end
end
save('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial.mat','Pval_pre','Reject_pre','Pval_group','Reject_group',...
    'Pval_pre_region','Reject_pre_region','Pval_group_region','Reject_group_region','-append','-v7.3');
%% check pval
typei=1;groupi=1;regioni=5;align_typei=1;
a=squeeze(NUM_CR_UR_align_all(1,typei,align_typei,numL));mean(a)
b=squeeze(NUM_CR_UR_align_all(7,typei,align_typei,numL));mean(b)
[p,~]=ranksum(a,b)
a=(squeeze(Pval_pre{typei}(:,align_typei,groupi))<0.05);
b=squeeze(data(:,regioni,typei,numL));

typei=1;groupi=3;groupj=4;align_typei=1;ss=1;
p=squeeze(Pval_group_region{typei,groupi}(:,:,align_typei,groupj));
h=regionxx(p<0.05);
%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% typical case
typei=1;align_type=1;sessioni=8;
data=Response_amp_all_session_align_all;size(data)
a=mean1(data(:,1,sessioni,typei,align_type,numL))+0.1;sd=error1(data(:,1,sessioni,typei,align_type,numL));
b=mean1(data(:,sessioni,sessioni,typei,align_type,numL));sd2=error1(data(:,sessioni,sessioni,typei,align_type,numL));
h=figure('position',[33,432,300,300],'color','k');
yl=[-0.2 1]; x=[1:size(a,1)]*fs.ca;
patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(yl) min(yl) max(yl)  max(yl)],'b','Facealpha',0.4,'edgealpha',0);
shadedErrorBar(x,a,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on
shadedErrorBar(x,b,sd2,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on

plot(x,a,'b','linewidth',4);hold on;
plot(x,b,'r','linewidth',4);hold on;
set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
ylim(yl);xlim([15,60]*fs.ca);
%% baseline CR
region_idd=1:30;
data=fraction_CR_UR_align_all;
for align_typei=1
    for typei=6
        ss=1;
        xl=[-1,1];
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
            ss=2:length(sessionx)-1;
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            m=squeeze(mean(dataL,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        xl=[0 60];
        h=figure('position',[10,50,1200,900],'name',name,'color','w');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            ss=2;
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            sd=std(dataL,[],2,'omitnan')./sqrt(size(m,2));m=mean1(dataL);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));  hold on
            title(group_label{groupi},'color','k');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Baseline_CR')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Baseline_CR')),'\',name]);
        %close(h)
    end
end

%% trace
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
data=Response_amp_align_all;size(data)
for zzi=1
    for align_typei=1:3
        for typei=1:34;
            if (typei>=1 & typei<=4)
                ss=1:length(sessionx);clim=[1 60];
            elseif (typei>=5 & typei<=6)
                ss=2:length(sessionx)-1;clim=[1 60];
            elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
                ss=2:length(sessionx);clim=[1 60];
            elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
                ss=3:length(sessionx)-1;clim=[1 60];
            elseif typei>=31 & typei<=34
                ss=3:length(sessionx)-1;clim=[1 60];
            end
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
                name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
                h=figure('position',[33,432,1900,300],'name',name,'color','k');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
                for sessioni=1:length(ss)
                    dataL=squeeze(mean1(data(:,ss(sessioni),typei,align_typei,ind)));
                    sd=squeeze(error1(data(:,ss(sessioni),typei,align_typei,ind)));
                    x=[1:size(dataL,1)]*fs.ca;
                    %[start_index,y,~,~]=findonset(dataL,base_win,onset_win);
                    nexttile;
                    a=[-1 2];
                    patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(a) min(a) max(a)  max(a)],'b','Facealpha',0.4,'edgealpha',0);
                    shadedErrorBar(x,dataL,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on
                    plot(x,dataL,'r','linewidth',4);hold on;
                    %line(([frame.cs_start frame.cs_start]*fs.ca),[min(a) max(a)],'Color','k','LineStyle','--','linewidth',2);hold on;
                    %line(([frame.cs_end frame.cs_end]*fs.ca),[min(a) max(a)],'Color','k','LineStyle','--','linewidth',2);hold on;
                    title(sessionx(ss(sessioni)),'color','w');
                    set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
                    ylim(a);xlim([15,60]*fs.ca);
                end
                            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
                text(0.8,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
                h.InvertHardcopy = 'off';
                saveas(h,[checkpath(fullfile(savepath_eps,'Avg. Trace')),'\',name],savetype);
                savefig(h,[checkpath(fullfile(savepath_eps,'Avg. Trace')),'\',name]);
                close(h)
            end
        end
    end
end

data=Response_amp_all_session_align_all;size(data)
%clr_cmap =addcolorplus(277);clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 60];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
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
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
            h=figure('position',[33,50,1500,1300],'name',name,'color','k');p=[];
            t=tiledlayout(length(ss),length(sessionx));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                for sessionj=1:length(sessionx)
                    dataL=squeeze(mean1(data(:,sessionj,ss(sessioni),typei,align_typei,ind)));
                    sd=squeeze(error1(data(:,sessionj,ss(sessioni),typei,align_typei,ind)));
                    x=[1:size(dataL,1)]*fs.ca;
                    a=[-1 2];
                    nexttile;%xlabel(string(sessionx(sessionj)));
                    patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(a) min(a) max(a)  max(a)],'b','Facealpha',0.4,'edgealpha',0);
                    shadedErrorBar(x,dataL,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                    p(sessioni)=plot(x,dataL,'color','r','linewidth',4);hold on;      
                    set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
                    ylim(a);xlim([15,60]*fs.ca);     
                    if sessionj==1; ylabel(string(sessionx(ss(sessioni)))) ; end
                end
            end
            %legend(p,sessionx(ss),'color','w','location','NorthEastOutside','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none');
            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(0.8,0.2,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name]);
            close(h)
        end
    end
end
%% bubbleheatmap plot
clr_cmapp =addcolorplus(300);
data=fraction_CR_UR_align_all;size(data)
kk=1;
%fraction_CR_UR_align(session_align,:,typei,align_typei);CS_response_amp_align(sessioni,typei,align_typei)
for align_typei=1
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        clim=[1 2];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            clim=[1 max(ceil(max(mean(dataL,1,'omitnan'))/5)*5,clim(2))];
        end
        if isempty(clim) clim=[1 2];end
        for groupi=1:4
            name=[num2str(kk),'_',labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            %y=max(squeeze(max(max(data(:,:,typei,align_typei,:)))));
            %cdata=discretize(dataL,[1:1:100]);
            %imagesc(dataL',[0 50]);colormap(clr_cmap);colorbar;
            %             a=reshape(repmat([1:length(sessionx)],length(regionxx),1)',1,[]);
            %             b=reshape(repmat([1:length(regionxx)],length(sessionx),1),1,[]);
            %             c=reshape(mat2str(dataL),[],1);
            %             text(a,b,c)
            %set(gca,'ytick',[1:length(regionxx(region_id))],'yticklabel',regionxx(region_id),'xtick',[1:length(sessionx(ss))],'xticklabel',sessionx(ss),'XTickLabelRotation',45)
            h=bubbleheatmap(dataL',dataL',sessionx(ss),regionxx(region_id),clr_cmapp,clim,clim);
            a=[num2str(kk);string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(1.2,0.5,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot')),'\',name]);
            close(h)
            kk=kk+1;
        end
    end
end
%% bubbleheatmap plot with CR
clr_cmapp =addcolorplus(301);
data=fraction_CR_UR_align_all;size(data)%fraction_CR_UR_align(session_align,:,typei,align_typei);
cdata=region_CS_response_amp_align_all;%region_CS_response_amp_align(sessioni,:,typei,align_typei)
kk=1;
for align_typei=1
    for typei=1;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        clim=[1 2];clim2=[0.2,1];
        for groupi=[1,2,3,4]
            name=[num2str(kk),'_',labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            %a=reshape(squeeze(cdata(ss,region_id,typei,align_typei,ind)),[],length(ind));
            %dataL2=a;%dataL2=normalize(a,1);
            %dataL2=mean1(reshape(dataL2,length(ss),[],length(ind)));
            dataL2=squeeze(mean1(cdata(ss,region_id,typei,align_typei,ind)));
            clim=[1 max(ceil(max(mean(dataL,1,'omitnan'))/5)*5,clim(2))];
            clim2=[min(floor(min(mean(dataL2,1))/0.5)*0.5,clim2(1)) max(ceil(max(mean(dataL2,1))/0.5)*0.5,clim2(2))];
        end
        if isempty(clim) clim=[1 2];end
        if isempty(clim2) clim2=[0.2 1];end
%         clim=[1 15];
%         clim2=[-0.8 0];
        %clr_cmapp2=flip(clr_cmapp);%
        for groupi=[1]
            %name=[num2str(kk),'_',labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
            name=[labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            %a=reshape(squeeze(cdata(ss,region_id,typei,align_typei,ind)),[],length(ind));
            %dataL2=a;%dataL2=normalize(a,1);
            %dataL2=mean1(reshape(dataL2,length(ss),[],length(ind)));
             dataL2=squeeze(mean1(cdata(ss,region_id,typei,align_typei,ind)));
            %dataL2=(dataL2-repmat(dataL2(1,:),length(ss),1))./(dataL2+repmat(dataL2(1,:),length(ss),1));
            %y=max(squeeze(max(max(data(:,:,typei,align_typei,:)))));
            %cdata=discretize(dataL,[1:1:100]);
            %set(gca,'ytick',[1:length(regionxx(region_id))],'yticklabel',regionxx(region_id),'xtick',[1:length(sessionx(ss))],'xticklabel',sessionx(ss),'XTickLabelRotation',45)

            h=bubbleheatmap(dataL',dataL2',sessionx(ss),regionxx(region_id),clr_cmapp2,clim,clim2);
            a=[num2str(kk);string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(1.2,0.5,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
%             saveas(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor')),'\',name],savetype);
%             savefig(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor')),'\',name]);
%             %close(h)
%             kk=kk+1;
        end
    end
end

clr_cmapp =addcolorplus(301);
data=fraction_CR_UR_align_all;size(data)%fraction_CR_UR_align(session_align,:,typei,align_typei);
cdata=region_CS_response_amp_align_all;%region_CS_response_amp_align(sessioni,:,typei,align_typei)
kk=1;
for align_typei=1
    for typei=1;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        clr_cmapp2=clr_cmapp;
        for groupi=[1]
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
            for regioni=1:30
                name=strcat(labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi},'_',string(regionxx(regioni)));
                dataL=mean1(squeeze(data(ss,regioni,typei,align_typei,ind)))*100;
                if sum(~isnan(dataL))~=0
                    dataL2=squeeze(mean1(cdata(ss,regioni,typei,align_typei,ind)));
                    clim=[1,ceil(max(dataL)*1.1)];
                    clim2=[0,ceil(max(dataL2)*1.1/0.2)*0.2];
                    h=bubbleheatmap2(dataL',dataL2',sessionx(ss),regionxx(regioni),clr_cmapp2,clim,clim2);
                    %             a=[num2str(kk);string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei});string(regionxx(regioni))];
                    %             text(1.2,0.5,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
                    h.InvertHardcopy = 'off';
                    saveas(h,strcat(checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor_region')),'\',name),savetype);
                    savefig(h,strcat(checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor_region')),'\',name));
                    close(h)
                end
            end
        end
    end
end

%% increased index
region_idd=1:30;
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','w');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));  hold on
            title(group_label{groupi},'color','k');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_mean')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_mean')),'\',name]);
        close(h)
    end
end
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        xl=[-10 10];
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];xx={};
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            xx{groupi}=m;
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2))/1;m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));            hold on
            title(group_label{groupi},'color','k');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum')),'\',name]);
                close(h)

    end
end
%分group�?
clr_map=[1,0,0;0,0,1;0.5,0.5,0.5];
for align_typei=1
    for typei=17
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xx={};sd=[];m=[];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);
            xx{groupi}=squeeze(sum(x,1,'omitnan')); sd(:,groupi)=nanstd(xx{groupi},[],2)./sqrt(size(xx{groupi},2)*1);m(:,groupi)=mean1(xx{groupi});
        end
        xl=[-1 3];
        name=[labels_typei{typei},'_',align_label{align_typei}];
        h=figure('position',[10,50,736,1000],'name',name);t=tiledlayout(5,6);t.TileSpacing = 'compact';t.Padding = 'compact';
        groupi=[1,3];lab=categorical({'Learner','Non-learner'});lab = reordercats(lab,{'Learner','Non-learner'});
       for regioni=1:length(regionxx)
            nexttile; 
            myStackedBarwithErrorbar(lab,m(regioni,groupi),sd(regioni,groupi),clr_cmap,[],[]);
            title(regionxx(regioni),'color','k');
            ylim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        hh=[];p=[];
        for ii=1:30
            [p(ii,:),hh(ii,:)]=ranksum(xx{1}(ii,:),xx{4}(ii,:));
        end
        regionxx(find(hh==1))
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum_group')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum_group')),'\',name]);
        %close(h)
    end
end
data=region_CS_response_amp_align_all;
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,[],[]);            hold on
            title(group_label{groupi},'color','k');
            %set(gca,'YTick',[1:length(region_id)],'Xticklabel',regionxx(region_id));
            xlim(xl);
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_mean')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_mean')),'\',name]);
        close(h)
    end
end
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));            hold on
            title(group_label{groupi},'color','k');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','k','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_sum')),'\',name],savetype);
       % savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_sum')),'\',name]);
                close(h)

    end
end
%% onset heatmap(先按鱼avg trace再算onset)
issort=1;region_idd=region_id;
data=Response_amp_region_align_all;%Response_amp_regddion_align_all;
size(data)
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        for groupi=1
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
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
            h=figure('position',[33,432,1900,500],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                dataL=squeeze(mean1(data(:,ss(sessioni),region_idd,typei,align_typei,ind)));
                [start_index,y,m,sd]=findonset(dataL,base_win,onset_win);
                %sepplot([1:frame.per_cycle],y,addcolorplus(332));hold on;
                nexttile;
                if issort==1
                    [start_index_s,I]=sort(start_index,'ascend');
                    imagesc(y(:,I)',[-1,2]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                    a=regionxx(region_idd); set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                    set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
                else
                    imagesc(y',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index) scatter(start_index,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[],'yticklabel',[]);
                    
                    if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                    else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
                end
                if sessioni==length(ss) colorbar; end
                line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
            end
                    h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Onset_of_regiontrace')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Onset_of_regiontrace')),'\',name]);
        close(h)
        end
    end
end
%% onset heatmap singlefish(先算onset再按鱼avg trace)
issort=1;
data=Response_amp_region_align_all;size(data)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 |(typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
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
            name=[align_label{align_typei},'_',group_label{groupi},'_Onset_of_singlefish_',labels_typei_ADJ{typei}];
            h=figure('position',[1,1,1900,400],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                start_index=[]; y=[];
                dataL=squeeze(data(:,ss(sessioni),region_idd,typei,align_typei,ind));
                for fishi=1:length(ind)
                    [a,b,~,~]=findonset(squeeze(dataL(:,:,fishi)),base_win,onset_win);
                    start_index(:,fishi)=a;
                    y(:,:,fishi)=b;
                end
                start_index_m=mean1(start_index);y_m=mean1(y);
                nexttile;
                if issort==1
                    [start_index_s,I]=sort(start_index_m,'ascend');
                    imagesc(y_m(:,I)',[0,0.2]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index_m) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                    set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
                else
                    imagesc(y_m',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index_m) scatter(start_index_m,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[],'yticklabel',[]);
                    if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                    else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
                end
                title(num2str(sessioni));
                if sessioni==length(ss) colorbar; end
                line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% onset heatmap singletrial
data=Response_amp_region_all;size(data)
for typei=1:16;
    ss=1:trial.total;
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
        name=[group_label{groupi},'_Onset_of_singletrial_',labels_typei_ADJ{typei}];
        h=figure('position',[1,1,1900,1000],'name',name,'color','w');t=tiledlayout(3,length(ss)/3);t.TileSpacing = 'compact';t.Padding = 'compact';
        for sessioni=1:length(ss)
            dataL=squeeze(mean1(data(:,ss(sessioni),region_idd,typei,ind)));
            [start_index,y,m,sd]=findonset(dataL,base_win,onset_win);
            nexttile;
            if issort==1
                [start_index_s,I]=sort(start_index,'ascend');
                imagesc(y(:,I)',[0,0.2]);colormap(addcolorplus(332));hold on;
                if ~isempty(start_index) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                a=regionxx(region_idd); set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
            else
                imagesc(y',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                if ~isempty(start_index) scatter(start_index,1:length(region_idd),15,'r','filled');end
                set(gca,'ytick',[],'yticklabel',[]);
                
                if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
            end
            title(num2str(sessioni));
            if sessioni==length(ss) colorbar; end
            line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
            line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
            set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
        end
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end

%% flatuate before CS
data=Response_amp_region_align_all;size(data)
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 |(typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
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
            name=[align_label{align_typei},'_',group_label{groupi},'_Fluctucate_befCS_of_',labels_typei{typei}];
            h=figure('position',[33,432,1900,500],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                dataL=squeeze(data(:,ss(sessioni),region_idd,typei,align_typei,ind));
                fluctucate=[];
                for fishi=1:length(ind)
                    [~,~,~,sd]=findonset(squeeze(dataL(:,:,fishi)),base_win,onset_win);
                    fluctucate(:,fishi)=sd;
                end
                sd=error1(fluctucate);m=mean1(fluctucate);
                [mm,I]=sort(m,'ascend');
                nexttile,
                myStackedBarwithErrorbarh(regionxx(region_idd(I)),mm',sd(I)',clr_cmap,sessionx(ss),regionxx(region_idd(I)));
                xlim([0 0.1]);
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% fraction_with Pval
data=NUM_CR_UR_align_all;
for align_typei=1
    for typei=[1,6,7,8,16,17];
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        name=[align_label{align_typei},'_BarwithPval_',labels_typei{typei(1)}];
        h=figure('position',[33,129,1721,803],'name',name,'color','w');
        tiledlayout(2,1);x=[];y=[];sig=[];sig2=[];
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
            dataL=squeeze(data(ss,typei,align_typei,ind))*100;
           a=dataL;% 
           a=(dataL-repmat(dataL(1,:),8,1))./repmat(dataL(1,:),8,1);
            x=cat(1,x,mean1(a)'); y=cat(1,y,error1(a)');
            p=Pval_pre{typei}(ss,align_typei,groupi);
            ind=find(p<0.05);
            a=[ones(length(ind),1)*groupi,ones(length(ind),1),ones(length(ind),1)*groupi,ind];
            sig=cat(1,sig,a);
            p=squeeze(Pval_group{typei,groupi}(ss,align_typei,:));
            [i,j]=find(p<0.05); a=[ones(length(i),1)*groupi,ones(length(i),1),j,i];
            sig2=cat(1,sig2,a);
        end
        a=[x(2,:);y(2,:);length(numL)*ones(1,8);x(3,:);y(3,:);length(numNL)*ones(1,8);x(4,:);y(4,:);length(numC)*ones(1,8)]';
        nexttile;mygroupbarwithsigline2(x,y,clr_group,sig,group_label);
        nexttile;mygroupbarwithsigline(x,y,clr_group,sig2,group_label);
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end
data=fraction_CR_UR_align_all;
for align_typei=1
    for typei=[1,6,7,8,16,17];
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        for regioni=1:length(region_id)
            name=[align_label{align_typei},'_',char(regionxx(regioni)),'_',labels_typei{typei(1)}];
            h=figure('position',[33,129,1721,803],'name',name,'color','w');
            tiledlayout(2,1);x=[];y=[];sig=[];sig2=[];
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
                dataL=squeeze(data(ss,region_id(regioni),typei,align_typei,ind));
                x=cat(1,x,mean1(dataL)'); y=cat(1,y,error1(dataL)'./sqrt(length(ind)));
                p=Pval_pre_region{typei}(ss,region_id(regioni),align_typei,groupi);
                ind=find(p<0.05);
                a=[ones(length(ind),1)*groupi,ones(length(ind),1),ones(length(ind),1)*groupi,ind];
                sig=cat(1,sig,a);
                p=squeeze(Pval_group_region{typei,groupi}(ss,region_id(regioni),align_typei,:));
                [i,j]=find(p<0.05); a=[ones(length(i),1)*groupi,ones(length(i),1),j,i];
                sig2=cat(1,sig2,a);
            end
            nexttile;mygroupbarwithsigline(x,y,clr_group,sig,group_label);
            nexttile;mygroupbarwithsigline(x,y,clr_group,sig2,group_label);
            saveas(h, [checkpath(fullfile(savepath_eps,'BarwithPval')),'\',name],savetype);
            close(h)
        end
    end
end
%% fraction
typei=[1,3];
data=NUM_CR_UR_align_all;
for align_typei=1:3;
    dataL2=squeeze(data(:,typei,align_typei,:));
    name=[align_label{align_typei},'_Num_of_',labels_typei_ADJ{typei(1)}];
    h=figure('position',[33,129,1721,403],'name',name,'color','w');
    tiledlayout(1,4);
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
        dataL=squeeze(dataL2(:,:,ind));
        ax1 = nexttile;
        b=bar(squeeze(mean1(dataL)),'FaceColor','flat','FaceAlpha',1,'BarWidth',1);hold on;
        for ii=1:length(typei)
            b(ii).CData = clr_group(typei(ii),:);
            errorbar(b(ii).XEndPoints,squeeze(mean1(dataL(:,ii,:))),squeeze(error1(dataL(:,ii,:))),'LineStyle','none','Color','k', 'LineWidth',.8);hold on;
        end
        set(gca,'fontsize',16,'XTickLabel',sessionx,'XTickLabelRotation',45);title(group_label{groupi});
        ylim([0 0.5]);if groupi==4; legend(b,{labels_typei_ADJ{typei}},'Interpreter','none') ;end
    end
    saveas(h,[savepath_eps,name],savetype);
end
close all;
%% Num plot
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2;clim=[1 60];
        elseif typei<=2
            ss=1;clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Fraction_of pre cond_',labels_typei_ADJ{typei}];
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
            
            %h=figure('position',[1040,90,431,856]);
            %dataL=squeeze(mean(data(ss,region_id,typei,align_typei,ind),1,'omitnan'))*100;
            dataL=squeeze(data(ss,region_id,typei,align_typei,ind))*100;
            m=mean1(dataL)';sd=error1(dataL)';
            h=myStackedBarwithErrorbarh(regionxx(region_id),m,sd,clr_cmap,[],regionxx(region_id));
            xlim([0,70]);saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Fraction_of all cond_',labels_typei_ADJ{typei}];
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
            dataL=squeeze(mean(data(ss,region_id,typei,align_typei,ind),1,'omitnan'))*100;
            %dataL=squeeze(data(ss,region_id,typei,align_typei,ind))*100;
            m=mean1(dataL)';sd=error1(dataL)';
            h=myStackedBarwithErrorbarh(regionxx(region_id),m,sd,clr_cmap,[],regionxx(region_id));
            xlim([0,60]);saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% linear fit
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Linearfit_of',labels_typei_ADJ{typei}];
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
            ss=2:length(sessionx);
            dataL=squeeze(data(:,:,typei,align_typei,ind))*1;
            h=figure('position',[1,41,1920,962]);t=tiledlayout(h,length(regionxx)/6,6);t.TileSpacing = 'compact';t.Padding = 'compact';
            for regioni=1:length(regionxx)
                x=repmat([ss-1]*0.1,length(ind),1)';
                y=squeeze(dataL(ss,regioni,:));
                x(isnan(y))=[];y(isnan(y))=[];
                nexttile;
                if ~isempty(x)
                    mylineplotofregression(reshape(x,1,[]),reshape(y,1,[]),sessionx(ss),'Fraction',regionxx(regioni))
                end
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% US trend index
data=fraction_CR_UR_align_all;
for align_typei=3
    for typei=3;
        name=[align_label{align_typei},'_boxplot_sum_of_',labels_typei_ADJ{typei}];
        h=figure('position',[1,41,1920,962]);t=tiledlayout(h,length(regionxx)/6,6);t.TileSpacing = 'compact';t.Padding = 'compact';
        for regioni=1:30;
            dataL={};
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
                ss=2:length(sessionx)-1;
                dataL{groupi}=squeeze(data(ss,regioni,typei,align_typei,ind))'*100;%dataL=mean1(dataL);
            end
            p1=[0:6:(length(ss)-1)*6]+0.5;
            nexttile;myGroupedFilledBoxPlot(dataL{1},p1,dataL{2},p1+1,dataL{3},p1+2,dataL{4},p1+3,[clr_group([4,3,2,1],:)],sessionx(ss),group_label);
            title(regionxx(regioni))
            %h=bubbleheatmap(dataL',dataL',sessionx(ss),regionxx(region_id),clr_cmap,[1 50]);
        end
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end
%% mapback
v=[90,90];
radium=floor(15/0.66);
data=supervoltemp_type_align_all;size(data)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Location_of_',labels_typei_ADJ{typei}];
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
            for sessioni=1:length(ss)
                dataL2=[];
%% demo
num=9;res=[0.66,0.66,10];v=[90,90];radium=floor(15/0.66);
load(fullfile(path_all{num},'para.mat'));load(fullfile(path_all{num},'env.mat'));
load(fullfile(path_all{num}, 'activities_dfdf_align.mat'),'align_win','align_win_go','align_win_nogo')
load(fullfile(path_all{num},'singletrial','CR_ind_summary_singletrial.mat'),'CR_ind_up');
supervoxel_raw=env.supervoxel;
for ii=unique(supervoxel_raw(:,3))'
    id=find(supervoxel_raw(:,3)==ii);
    supervoxel_raw(id,3)=(ii-1)*10/0.66+1;
end
supervoxel_raw(:,1)=env.width-supervoxel_raw(:,1);
for ii=1:size(align_win,2)
    a=align_win(:,ii);a(isnan(a))=[];
    h=figure('position',[7,154,1910,824]);t=tiledlayout(ceil((length(a)/3)),3);t.TileSpacing = 'compact';t.Padding = 'compact';
    for jj=1:length(a)
        ind=find(CR_ind_up(:,a(jj),iscutmov)==1);leg=[num2str(a(jj))];
        nexttile,count_spatial_location_density(supervoxel_raw(ind,:),supervoxel_raw(:,:), radium);view(v(1),v(2));title(leg);
        set(gca,'fontsize',16)
    end
end



%%
clc;clear all
savepath='I:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
load([savepath,'/Path.mat'])
%% select fish
numL=num{1};
numC=num{2}+10;numNL=num{3}+10+10;numFL=num{4}+10+10+9;
sessionx = categorical({'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
sessionx = reordercats(sessionx,{'Pre Cond','Cond.1','Cond.2','Cond.3','Cond.4','Cond.5','Cond.6','Post Cond'});
group_label={'Learner (n=8)','Faded-Learner (n=6)','Non-Learner (n=9)','Control (n=8)'};
align_label={'All trials','Acting trials','Non-acting trials'};
load(fullfile(Path{1}{1},'CR_ind_summary_align.mat'),'Label_region');
load(fullfile(Path{1}{1},'para.mat'));
regionxx=categorical(Label_region{1,1}(:,1));regionxx = reordercats(regionxx,Label_region{1,1}(:,1));
mean1= @(x)(mean(x,ndims(x),'omitnan'));
error1= @(x)(std(x,[],ndims(x),'omitnan'));%error1= @(x)(std(x,[],ndims(x),'omitnan')./sqrt(length(x)));
error2= @(x)(mean(x,ndims(x),'omitnan')-std(x,[],ndims(x),'omitnan'));
clr_group=[colorplus(389);colorplus(253);colorplus(233);colorplus(448)];
clr_cmap =addcolorplus(300);addcolorplus(336);%cmap = GetColormap('jet',81);
clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
clr_CStypes = addcolorplus(275);%clr_CStypes= ColorMap(clr_CStypes,3);clr_CStypes = flipud(clr_CStypes);
clr_UStypes = addcolorplus(272);%clr_UStypes= ColorMap(clr_UStypes,3);clr_UStypes = flipud(clr_UStypes);
region_id=setdiff(1:length(regionxx),[1,2,18:20,23,24,27,28]);
%region_id=1:30;
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
%% loadpath
res=[0.66,0.66,10];
savepath_eps=checkpath('I:\2.paper_fig_eps_2024\singletrial_plus_3sd\');
load('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial.mat',...
    'NUM_CR_UR_align_all','fraction_CR_UR_align_all','region_CS_response_amp_align_all',...
    'Response_amp_region_align_all');
load('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial_append.mat','Response_amp_all_session_align_all');
%load('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial.mat','Response_amp_align_all','Response_amp_all_session_align_all');
%load('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial_2sd.mat','supervoltemp_type_align_all');
savetype='jpg';
temp_env=load('I:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);


%% p-value
Pval_pre_region={};Reject_pre_region={};
data=fraction_CR_UR_align_all;%fraction_CR_UR_align_all(session_align,:,typei,align_typei,fishi)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx);clim=[1 80];
        end
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
            data1=squeeze(data(ss,:,typei,align_typei,ind));%data1(setdiff(1:size(data,1),ss),:,:)=nan;
            for regioni=1:size(data1,2)
                for sessioni=1:length(ss)
                    a=squeeze(data1(1,regioni,:));b=squeeze(data1(sessioni,regioni,:));
                    if ~(sum(isnan(a))==size(data1,3) | sum(isnan(b))==size(data1,3))
                        [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                        Pval_pre_region{typei}(ss(sessioni),regioni,align_typei,groupi)=p;
                        Reject_pre_region{typei}(ss(sessioni),regioni,align_typei,groupi)=p;
                    end
                end
            end
        end
    end
end
Pval_group_region={};Reject_group_region={};
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx);clim=[1 80];
        end
        for groupi=1:4
            switch groupi
                case 1
                    indi=numL;
                case 4
                    indi=numC;
                case 3
                    indi=numNL;
                case 2
                    indi=numFL;
            end
            for groupj=1:4
                switch groupj
                    case 1
                        indj=numL;
                    case 4
                        indj=numC;
                    case 3
                        indj=numNL;
                    case 2
                        indj=numFL;
                end
                data1=squeeze(data(ss,:,typei,align_typei,indi));
                data2=squeeze(data(ss,:,typei,align_typei,indj));
                for regioni=1:size(data1,2)
                    for sessioni=1:length(ss)
                        a=squeeze(data1(sessioni,regioni,:));b=squeeze(data2(sessioni,regioni,:));
                        if ~(sum(isnan(a))==size(data1,3) | sum(isnan(b))==size(data2,3))
                            [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                            Pval_group_region{typei,groupi}(ss(sessioni),regioni,align_typei,groupj)=p;
                            Reject_group_region{typei,groupi}(ss(sessioni),regioni,align_typei,groupj)=p;
                        end
                    end
                end
            end
        end
    end
end

Pval_pre={};Reject_pre={};Pval_group={};Reject_group={};
data=NUM_CR_UR_align_all;%NUM_CR_UR_align(session_align,typei,align_typei,fishi)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx);clim=[1 80];
        end
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
            data1=squeeze(data(ss,typei,align_typei,ind));%data1(setdiff(1:size(data,1),ss),:,:)=nan;
            for sessioni=1:length(ss)
                a=squeeze(data1(1,:));b=squeeze(data1(sessioni,:));
                if ~(sum(isnan(a))==size(data1,2) | sum(isnan(b))==size(data1,2))
                    [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                    Pval_pre{typei}(ss(sessioni),align_typei,groupi)=p;
                    Reject_pre{typei}(ss(sessioni),align_typei,groupi)=p;
                end
            end
        end
    end
end
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx);clim=[1 80];
        end
        for groupi=1:4
            switch groupi
                case 1
                    indi=numL;
                case 4
                    indi=numC;
                case 3
                    indi=numNL;
                case 2
                    indi=numFL;
            end
            for groupj=1:4
                switch groupj
                    case 1
                        indj=numL;
                    case 4
                        indj=numC;
                    case 3
                        indj=numNL;
                    case 2
                        indj=numFL;
                end
                data1=squeeze(data(ss,typei,align_typei,indi));
                data2=squeeze(data(ss,typei,align_typei,indj));
                for sessioni=1:length(ss)
                    a=squeeze(data1(sessioni,:));b=squeeze(data2(sessioni,:));
                    if ~(sum(isnan(a))==size(data1,2) | sum(isnan(b))==size(data2,2))
                        [p,h]=ranksum(a(~isnan(a)),b(~isnan(b)));
                        Pval_group{typei,groupi}(ss(sessioni),align_typei,groupj)=p;
                        Reject_group{typei,groupi}(ss(sessioni),align_typei,groupj)=p;
                    end
                end
            end
        end
    end
end
% save('H:\1.Test US\5.fear conditioning behavioral data\Summary_all_CRUR_singletrial_2sd.mat','Pval_pre','Reject_pre','Pval_group','Reject_group',...
%     'Pval_pre_region','Reject_pre_region','Pval_group_region','Reject_group_region','-append','-v7.3');
%% check pval
typei=1;groupi=1;regioni=5;align_typei=1;
a=squeeze(NUM_CR_UR_align_all(1,typei,align_typei,numL));mean(a)
b=squeeze(NUM_CR_UR_align_all(7,typei,align_typei,numL));mean(b)
[p,~]=ranksum(a,b)
a=(squeeze(Pval_pre{typei}(:,align_typei,groupi))<0.05);
b=squeeze(data(:,regioni,typei,numL));

typei=1;groupi=3;groupj=4;align_typei=1;ss=1;
p=squeeze(Pval_group_region{typei,groupi}(:,:,align_typei,groupj));
h=regionxx(p<0.05);
%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% typical case
typei=1;align_type=1;sessioni=8;
data=Response_amp_all_session_align_all;size(data)
a=mean1(data(:,1,sessioni,typei,align_type,numL))+0.1;sd=error1(data(:,1,sessioni,typei,align_type,numL));
b=mean1(data(:,sessioni,sessioni,typei,align_type,numL));sd2=error1(data(:,sessioni,sessioni,typei,align_type,numL));
h=figure('position',[33,432,300,300],'color','k');
yl=[-0.2 1]; x=[1:size(a,1)]*fs.ca;
patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(yl) min(yl) max(yl)  max(yl)],'b','Facealpha',0.4,'edgealpha',0);
shadedErrorBar(x,a,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on
shadedErrorBar(x,b,sd2,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on

plot(x,a,'b','linewidth',4);hold on;
plot(x,b,'r','linewidth',4);hold on;
set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
ylim(yl);xlim([15,60]*fs.ca);

%% trace
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
data=Response_amp_align_all;size(data)
for zzi=[]
    for align_typei=1:3
        for typei=1:34;
            if (typei>=1 & typei<=4)
                ss=1:length(sessionx);clim=[1 60];
            elseif (typei>=5 & typei<=6)
                ss=2:length(sessionx)-1;clim=[1 60];
            elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
                ss=2:length(sessionx);clim=[1 60];
            elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
                ss=3:length(sessionx)-1;clim=[1 60];
            elseif typei>=31 & typei<=34
                ss=3:length(sessionx)-1;clim=[1 60];
            end
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
                name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
                h=figure('position',[33,432,1900,300],'name',name,'color','k');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
                for sessioni=1:length(ss)
                    dataL=squeeze(mean1(data(:,ss(sessioni),typei,align_typei,ind)));
                    sd=squeeze(error1(data(:,ss(sessioni),typei,align_typei,ind)));
                    x=[1:size(dataL,1)]*fs.ca;
                    %[start_index,y,~,~]=findonset(dataL,base_win,onset_win);
                    nexttile;
                    a=[-1 2];
                    patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(a) min(a) max(a)  max(a)],'b','Facealpha',0.4,'edgealpha',0);
                    shadedErrorBar(x,dataL,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.4);hold on
                    plot(x,dataL,'r','linewidth',4);hold on;
                    %line(([frame.cs_start frame.cs_start]*fs.ca),[min(a) max(a)],'Color','k','LineStyle','--','linewidth',2);hold on;
                    %line(([frame.cs_end frame.cs_end]*fs.ca),[min(a) max(a)],'Color','k','LineStyle','--','linewidth',2);hold on;
                    title(sessionx(ss(sessioni)),'color','w');
                    set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
                    ylim(a);xlim([15,60]*fs.ca);
                end
                            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
                text(0.8,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
                h.InvertHardcopy = 'off';
                saveas(h,[checkpath(fullfile(savepath_eps,'Avg. Trace')),'\',name],savetype);
                savefig(h,[checkpath(fullfile(savepath_eps,'Avg. Trace')),'\',name]);
                close(h)
            end
        end
    end
end

data=Response_amp_all_session_align_all;size(data)
%clr_cmap =addcolorplus(277);clr_session = clr_cmap(1:size(clr_cmap,1)/length(sessionx):end,:);
for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 60];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
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
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
            h=figure('position',[33,50,1500,1300],'name',name,'color','k');p=[];
            t=tiledlayout(length(ss),length(sessionx));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                for sessionj=1:length(sessionx)
                    dataL=squeeze(mean1(data(:,sessionj,ss(sessioni),typei,align_typei,ind)));
                    sd=squeeze(error1(data(:,sessionj,ss(sessioni),typei,align_typei,ind)));
                    x=[1:size(dataL,1)]*fs.ca;
                    a=[-1 2];
                    nexttile;%xlabel(string(sessionx(sessionj)));
                    patch([frame.cs_start frame.cs_end frame.cs_end frame.cs_start ]*fs.ca,[min(a) min(a) max(a)  max(a)],'b','Facealpha',0.4,'edgealpha',0);
                    shadedErrorBar(x,dataL,sd,'lineprops',{[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.2);hold on
                    p(sessioni)=plot(x,dataL,'color','r','linewidth',4);hold on;      
                    set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');
                    ylim(a);xlim([15,60]*fs.ca);     
                    if sessionj==1; ylabel(string(sessionx(ss(sessioni)))) ; end
                end
            end
            %legend(p,sessionx(ss),'color','w','location','NorthEastOutside','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none');
            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(0.8,0.2,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Avg. Trace_singlepanel')),'\',name]);
            close(h)
        end
    end
end
%% bubbleheatmap plot
clr_cmapp =addcolorplus(301);
data=fraction_CR_UR_align_all;size(data)
kk=1;
%fraction_CR_UR_align(session_align,:,typei,align_typei);CS_response_amp_align(sessioni,typei,align_typei)
for align_typei=1
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        clim=[1 2];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            clim=[1 max(ceil(max(mean(dataL,1,'omitnan'))/5)*5,clim(2))];
        end
        if isempty(clim) clim=[1 2];end
        for groupi=1:4
            name=[num2str(kk),'_',labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            %y=max(squeeze(max(max(data(:,:,typei,align_typei,:)))));
            %cdata=discretize(dataL,[1:1:100]);
            %imagesc(dataL',[0 50]);colormap(clr_cmap);colorbar;
            %             a=reshape(repmat([1:length(sessionx)],length(regionxx),1)',1,[]);
            %             b=reshape(repmat([1:length(regionxx)],length(sessionx),1),1,[]);
            %             c=reshape(mat2str(dataL),[],1);
            %             text(a,b,c)
            %set(gca,'ytick',[1:length(regionxx(region_id))],'yticklabel',regionxx(region_id),'xtick',[1:length(sessionx(ss))],'xticklabel',sessionx(ss),'XTickLabelRotation',45)
            h=bubbleheatmap(dataL',dataL',sessionx(ss),regionxx(region_id),clr_cmapp,clim,clim);
            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(1.2,0.5,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot')),'\',name]);
            close(h)
        end
    end
end
%% bubbleheatmap plot with CR
data=fraction_CR_UR_align_all;size(data)%fraction_CR_UR_align(session_align,:,typei,align_typei);
cdata=region_CS_response_amp_align_all;%region_CS_response_amp_align(sessioni,:,typei,align_typei)
for align_typei=1:3
    for typei=1:34;
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        clim=[1 2];clim2=[0.2,1];
        for groupi=[1,2,3]
            name=[labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            a=reshape(squeeze(cdata(ss,region_id,typei,align_typei,ind)),[],length(ind));
            dataL2=a;%dataL2=normalize(a,1);
            dataL2=mean1(reshape(dataL2,length(ss),[],length(ind)));
            clim=[1 max(ceil(max(mean(dataL,1,'omitnan'))/5)*5,clim(2))];
            clim2=[min(floor(min(mean(dataL2,1))/0.5)*0.5,clim2(1)) max(ceil(max(mean(dataL2,1))/0.5)*0.5,clim2(2))];
        end
         if isempty(clim) clim=[1 2];end
          if isempty(clim2) clim2=[0.2 1];end
        for groupi=[1,2,3]
            name=[labels_typei{typei},'_',align_label{align_typei},'_',group_label{groupi}];
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
            dataL=mean1(squeeze(data(ss,region_id,typei,align_typei,ind)))*100;
            a=reshape(squeeze(cdata(ss,region_id,typei,align_typei,ind)),[],length(ind));
            dataL2=a;%dataL2=normalize(a,1);
            dataL2=mean1(reshape(dataL2,length(ss),[],length(ind)));
            %dataL2=(dataL2-repmat(dataL2(1,:),length(ss),1))./(dataL2+repmat(dataL2(1,:),length(ss),1));
            %y=max(squeeze(max(max(data(:,:,typei,align_typei,:)))));
            %cdata=discretize(dataL,[1:1:100]);
            %set(gca,'ytick',[1:length(regionxx(region_id))],'yticklabel',regionxx(region_id),'xtick',[1:length(sessionx(ss))],'xticklabel',sessionx(ss),'XTickLabelRotation',45)
            
            h=bubbleheatmap(dataL',dataL2',sessionx(ss),regionxx(region_id),clr_cmapp,clim,clim2);
            a=[string(labels_typei{typei});string(group_label{groupi});string(align_label{align_typei})];
            text(1.2,0.5,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
            h.InvertHardcopy = 'off';
            saveas(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor')),'\',name],savetype);
            savefig(h,[checkpath(fullfile(savepath_eps,'Fraction_Bubbleplot_CRcolor')),'\',name]);
            close(h)
        end
    end
end
%% increased index
region_idd=1:30;
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));            hold on
            title(group_label{groupi},'color','w');
            set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13,'FontName','Times New Roman');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_mean')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_mean')),'\',name]);
        close(h)
    end
end
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));            hold on
            title(group_label{groupi},'color','w');
            set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13,'FontName','Times New Roman');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_sum')),'\',name]);
                close(h)

    end
end
data=region_CS_response_amp_align_all;
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(mean(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,[],[]);            hold on
            title(group_label{groupi},'color','w');
            %set(gca,'YTick',[1:length(region_id)],'Xticklabel',regionxx(region_id));
            set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13,'FontName','Times New Roman');
            xlim(xl);
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_mean')),'\',name],savetype);
        savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_mean')),'\',name]);
        close(h)
    end
end
for align_typei=1:3
    for typei=1:34
        if (typei>=1 & typei<=4)
            ss=1:length(sessionx);clim=[1 40];
        elseif (typei>=5 & typei<=6)
            ss=2:length(sessionx)-1;clim=[1 60];
        elseif (typei>=7 & typei<=12) | (typei>=19 & typei<=24)
            ss=2:length(sessionx);clim=[1 60];
        elseif (typei>=13 & typei<=18) | (typei>=25 & typei<=30)
            ss=3:length(sessionx)-1;clim=[1 60];
        elseif typei>=31 & typei<=34
            ss=3:length(sessionx)-1;clim=[1 60];
        end
        xl=[-1,1];
        for groupi=1:4
            name=[labels_typei{typei},'_',group_label{groupi},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            m=mean1(m);
            xl=[min(floor(min(m)/0.5)*0.5,xl(1)) max(ceil(max(m)/0.5)*0.5,xl(2))];
        end
        h=figure('position',[10,50,1200,900],'name',name,'color','k');t=tiledlayout(1,4);t.TileSpacing = 'compact';t.Padding = 'compact';
        name=[labels_typei{typei},'_',align_label{align_typei}];
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
            %ss=2:length(sessionx);
            dataL=squeeze(data(ss,region_idd,typei,align_typei,ind))*100;%dataL=mean1(dataL);
            dataL2=squeeze(data([ss(2:end) ss(end)],region_idd,typei,align_typei,ind))*100;%dataL2=mean1(dataL2);
            x=dataL2-dataL;x=x(1:end-1,:,:);m=squeeze(sum(x,1,'omitnan'));%sd=squeeze(std(x,[],3,'omitnan'))./size(x,3);
            sd=std(m,[],2,'omitnan')./sqrt(size(m,2));m=mean1(m);
            nexttile; myStackedBarwithErrorbarh(regionxx(region_idd),m',sd',clr_cmap,sessionx(ss(2:end)),regionxx(region_idd));            hold on
            title(group_label{groupi},'color','w');
            set(gca,'XColor','w','YColor','w','linewidth',2,'color','k');set(gca, 'FontSize',13,'FontName','Times New Roman');
            xlim(xl);
            %close(h)
        end
        a=[strcat(num2str(typei),'.',string(labels_typei{typei}));strcat([],string(align_label{align_typei}))];
        text(0.75,1,a,'color','w','FontSize',10,'FontName','Times New Roman','FontWeight','bold','Interpreter','none','Units','normalized');
        
        h.InvertHardcopy = 'off';
        saveas(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_sum')),'\',name],savetype);
       % savefig(h,[checkpath(fullfile(savepath_eps,'Increaseindex_amp_sum')),'\',name]);
                close(h)

    end
end
%% onset heatmap(先按鱼avg trace再算onset)
issort=1;
data=Response_amp_region_align_all;size(data)
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 |(typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
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
            name=[align_label{align_typei},'_',group_label{groupi},'_Onset_of_',labels_typei_ADJ{typei}];
            h=figure('position',[33,432,1900,500],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                dataL=squeeze(mean1(data(:,ss(sessioni),region_idd,typei,align_typei,ind)));
                [start_index,y,m,sd]=findonset(dataL,base_win,onset_win);
                %sepplot([1:frame.per_cycle],y,addcolorplus(332));hold on;
                nexttile;
                if issort==1
                    [start_index_s,I]=sort(start_index,'ascend');
                    imagesc(y(:,I)',[-0,0.2]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                    a=regionxx(region_idd); set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                    set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
                else
                    imagesc(y',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index) scatter(start_index,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[],'yticklabel',[]);
                    
                    if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                    else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
                end
                if sessioni==length(ss) colorbar; end
                line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% onset heatmap singlefish(先算onset再按鱼avg trace)
issort=1;
data=Response_amp_region_align_all;size(data)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 |(typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
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
            name=[align_label{align_typei},'_',group_label{groupi},'_Onset_of_singlefish_',labels_typei_ADJ{typei}];
            h=figure('position',[1,1,1900,400],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                start_index=[]; y=[];
                dataL=squeeze(data(:,ss(sessioni),region_idd,typei,align_typei,ind));
                for fishi=1:length(ind)
                    [a,b,~,~]=findonset(squeeze(dataL(:,:,fishi)),base_win,onset_win);
                    start_index(:,fishi)=a;
                    y(:,:,fishi)=b;
                end
                start_index_m=mean1(start_index);y_m=mean1(y);
                nexttile;
                if issort==1
                    [start_index_s,I]=sort(start_index_m,'ascend');
                    imagesc(y_m(:,I)',[0,0.2]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index_m) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                    set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
                else
                    imagesc(y_m',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                    if ~isempty(start_index_m) scatter(start_index_m,1:length(region_idd),15,'r','filled');end
                    set(gca,'ytick',[],'yticklabel',[]);
                    if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                    else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
                end
                title(num2str(sessioni));
                if sessioni==length(ss) colorbar; end
                line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
                set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% onset heatmap singletrial
data=Response_amp_region_all;size(data)
for typei=1:16;
    ss=1:trial.total;
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
        name=[group_label{groupi},'_Onset_of_singletrial_',labels_typei_ADJ{typei}];
        h=figure('position',[1,1,1900,1000],'name',name,'color','w');t=tiledlayout(3,length(ss)/3);t.TileSpacing = 'compact';t.Padding = 'compact';
        for sessioni=1:length(ss)
            dataL=squeeze(mean1(data(:,ss(sessioni),region_idd,typei,ind)));
            [start_index,y,m,sd]=findonset(dataL,base_win,onset_win);
            nexttile;
            if issort==1
                [start_index_s,I]=sort(start_index,'ascend');
                imagesc(y(:,I)',[0,0.2]);colormap(addcolorplus(332));hold on;
                if ~isempty(start_index) scatter(start_index_s,1:length(region_idd),15,'r','filled');end
                a=regionxx(region_idd); set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd(I)));
                set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca);
            else
                imagesc(y',[-0.1,0.3]);colormap(addcolorplus(332));hold on;
                if ~isempty(start_index) scatter(start_index,1:length(region_idd),15,'r','filled');end
                set(gca,'ytick',[],'yticklabel',[]);
                
                if sessioni==1 set(gca,'ytick',[1:length(region_idd)],'yticklabel',regionxx(region_idd));
                else  set(gca,'xtick',[0:25:frame.per_cycle],'xticklabel',[0:25:frame.per_cycle]*fs.ca); end
            end
            title(num2str(sessioni));
            if sessioni==length(ss) colorbar; end
            line(([frame.cs_start frame.cs_start]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
            line(([frame.cs_end frame.cs_end]),[0 length(region_idd)],'Color','w','LineStyle','--','linewidth',2);hold on;
            set(gca,'XColor','k','YColor','k','linewidth',1);set(gca, 'FontSize',13);set(gca,'FontName','Times New Roman');xlim([5,70]);
        end
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end

%% flatuate before CS
data=Response_amp_region_align_all;size(data)
pre_cs_time=[5 0];
base_win=floor(frame.cs_start-pre_cs_time(1)/fs.ca:frame.cs_start-pre_cs_time(2)/fs.ca-1);
onset_win=frame.cs_start-4:frame.per_cycle;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 |(typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
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
            name=[align_label{align_typei},'_',group_label{groupi},'_Fluctucate_befCS_of_',labels_typei_ADJ{typei}];
            h=figure('position',[33,432,1900,500],'name',name,'color','w');t=tiledlayout(1,length(ss));t.TileSpacing = 'compact';t.Padding = 'compact';
            for sessioni=1:length(ss)
                dataL=squeeze(data(:,ss(sessioni),region_idd,typei,align_typei,ind));
                fluctucate=[];
                for fishi=1:length(ind)
                    [~,~,~,sd]=findonset(squeeze(dataL(:,:,fishi)),base_win,onset_win);
                    fluctucate(:,fishi)=sd;
                end
                sd=error1(fluctucate);m=mean1(fluctucate);
                [mm,I]=sort(m,'ascend');
                nexttile,
                myStackedBarwithErrorbarh(regionxx(region_idd(I)),mm',sd(I)',clr_cmap,sessionx(ss),regionxx(region_idd(I)));
                xlim([0 0.1]);
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% fraction_with Pval
data=NUM_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        name=[align_label{align_typei},'_BarwithPval_',labels_typei_ADJ{typei(1)}];
        h=figure('position',[33,129,1721,803],'name',name,'color','w');
        tiledlayout(2,1);x=[];y=[];sig=[];sig2=[];
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
            dataL=squeeze(data(ss,typei,align_typei,ind));
            x=cat(1,x,mean1(dataL)'); y=cat(1,y,error1(dataL)');
            p=Pval_pre{typei}(ss,align_typei,groupi);
            ind=find(p<0.05);
            a=[ones(length(ind),1)*groupi,ones(length(ind),1),ones(length(ind),1)*groupi,ind];
            sig=cat(1,sig,a);
            p=squeeze(Pval_group{typei,groupi}(ss,align_typei,:));
            [i,j]=find(p<0.05); a=[ones(length(i),1)*groupi,ones(length(i),1),j,i];
            sig2=cat(1,sig2,a);
        end
        nexttile;mygroupbarwithsigline(x,y,clr_group,sig,group_label);
        nexttile;mygroupbarwithsigline(x,y,clr_group,sig2,group_label);
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        for regioni=1:length(region_id)
            name=[align_label{align_typei},'_',char(regionxx(regioni)),'_',labels_typei_ADJ{typei(1)}];
            h=figure('position',[33,129,1721,803],'name',name,'color','w');
            tiledlayout(2,1);x=[];y=[];sig=[];sig2=[];
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
                dataL=squeeze(data(ss,region_id(regioni),typei,align_typei,ind));
                x=cat(1,x,mean1(dataL)'); y=cat(1,y,error1(dataL)');
                p=Pval_pre_region{typei}(ss,region_id(regioni),align_typei,groupi);
                ind=find(p<0.05);
                a=[ones(length(ind),1)*groupi,ones(length(ind),1),ones(length(ind),1)*groupi,ind];
                sig=cat(1,sig,a);
                p=squeeze(Pval_group_region{typei,groupi}(ss,region_id(regioni),align_typei,:));
                [i,j]=find(p<0.05); a=[ones(length(i),1)*groupi,ones(length(i),1),j,i];
                sig2=cat(1,sig2,a);
            end
            nexttile;mygroupbarwithsigline(x,y,clr_group,sig,group_label);
            nexttile;mygroupbarwithsigline(x,y,clr_group,sig2,group_label);
            saveas(h, [checkpath(fullfile(savepath_eps,'BarwithPval')),'\',name],savetype);
            close(h)
        end
    end
end
%% fraction
typei=[1,3];
data=NUM_CR_UR_align_all;
for align_typei=1:3;
    dataL2=squeeze(data(:,typei,align_typei,:));
    name=[align_label{align_typei},'_Num_of_',labels_typei_ADJ{typei(1)}];
    h=figure('position',[33,129,1721,403],'name',name,'color','w');
    tiledlayout(1,4);
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
        dataL=squeeze(dataL2(:,:,ind));
        ax1 = nexttile;
        b=bar(squeeze(mean1(dataL)),'FaceColor','flat','FaceAlpha',1,'BarWidth',1);hold on;
        for ii=1:length(typei)
            b(ii).CData = clr_group(typei(ii),:);
            errorbar(b(ii).XEndPoints,squeeze(mean1(dataL(:,ii,:))),squeeze(error1(dataL(:,ii,:))),'LineStyle','none','Color','k', 'LineWidth',.8);hold on;
        end
        set(gca,'fontsize',16,'XTickLabel',sessionx,'XTickLabelRotation',45);title(group_label{groupi});
        ylim([0 0.5]);if groupi==4; legend(b,{labels_typei_ADJ{typei}},'Interpreter','none') ;end
    end
    saveas(h,[savepath_eps,name],savetype);
end
close all;
%% Num plot
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2;clim=[1 60];
        elseif typei<=2
            ss=1;clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Fraction_of pre cond_',labels_typei_ADJ{typei}];
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
            
            %h=figure('position',[1040,90,431,856]);
            %dataL=squeeze(mean(data(ss,region_id,typei,align_typei,ind),1,'omitnan'))*100;
            dataL=squeeze(data(ss,region_id,typei,align_typei,ind))*100;
            m=mean1(dataL)';sd=error1(dataL)';
            h=myStackedBarwithErrorbarh(regionxx(region_id),m,sd,clr_cmap,[],regionxx(region_id));
            xlim([0,70]);saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Fraction_of all cond_',labels_typei_ADJ{typei}];
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
            dataL=squeeze(mean(data(ss,region_id,typei,align_typei,ind),1,'omitnan'))*100;
            %dataL=squeeze(data(ss,region_id,typei,align_typei,ind))*100;
            m=mean1(dataL)';sd=error1(dataL)';
            h=myStackedBarwithErrorbarh(regionxx(region_id),m,sd,clr_cmap,[],regionxx(region_id));
            xlim([0,60]);saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% linear fit
data=fraction_CR_UR_align_all;
for align_typei=1:3
    for typei=1:16;
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Linearfit_of',labels_typei_ADJ{typei}];
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
            ss=2:length(sessionx);
            dataL=squeeze(data(:,:,typei,align_typei,ind))*1;
            h=figure('position',[1,41,1920,962]);t=tiledlayout(h,length(regionxx)/6,6);t.TileSpacing = 'compact';t.Padding = 'compact';
            for regioni=1:length(regionxx)
                x=repmat([ss-1]*0.1,length(ind),1)';
                y=squeeze(dataL(ss,regioni,:));
                x(isnan(y))=[];y(isnan(y))=[];
                nexttile;
                if ~isempty(x)
                    mylineplotofregression(reshape(x,1,[]),reshape(y,1,[]),sessionx(ss),'Fraction',regionxx(regioni))
                end
            end
            saveas(h,[savepath_eps,name],savetype);
            close(h)
        end
    end
end
%% US trend index
data=fraction_CR_UR_align_all;
for align_typei=3
    for typei=3;
        name=[align_label{align_typei},'_boxplot_sum_of_',labels_typei_ADJ{typei}];
        h=figure('position',[1,41,1920,962]);t=tiledlayout(h,length(regionxx)/6,6);t.TileSpacing = 'compact';t.Padding = 'compact';
        for regioni=1:30;
            dataL={};
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
                ss=2:length(sessionx)-1;
                dataL{groupi}=squeeze(data(ss,regioni,typei,align_typei,ind))'*100;%dataL=mean1(dataL);
            end
            p1=[0:6:(length(ss)-1)*6]+0.5;
            nexttile;myGroupedFilledBoxPlot(dataL{1},p1,dataL{2},p1+1,dataL{3},p1+2,dataL{4},p1+3,[clr_group([4,3,2,1],:)],sessionx(ss),group_label);
            title(regionxx(regioni))
            %h=bubbleheatmap(dataL',dataL',sessionx(ss),regionxx(region_id),clr_cmap,[1 50]);
        end
        saveas(h,[savepath_eps,name],savetype);
        close(h)
    end
end
%% mapback
v=[90,90];
radium=floor(15/0.66);
data=supervoltemp_type_align_all;size(data)
for align_typei=1:3
    for typei=1:16;
        if (typei>=4 & typei<=6) | (typei>=10 & typei<=13)
            ss=2:length(sessionx);clim=[1 60];
        elseif typei<=2
            ss=1:length(sessionx);clim=[1 60];
        elseif typei==3 | (typei>=7 & typei<=9) | (typei>=14 & typei<=16)
            ss=2:length(sessionx)-1;clim=[1 80];
        end
        for groupi=1:4
            name=[align_label{align_typei},'_',group_label{groupi},'_Location_of_',labels_typei_ADJ{typei}];
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
            for sessioni=1:length(ss)
                dataL2=[];
                for ii=1:length(ind)
                    dataL=data(ss(sessioni),typei,align_typei,ind(ii));
                    for kk=1:length(dataL{1})
                        dataL2=cat(1,dataL{1}{kk},dataL2);
                    end
                end
                kk=1;
                
                nexttile(2+kk,[2,1]);kk=kk+1;
                scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),8,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
                scatter3(dataL2(:,1),dataL2(:,2),dataL2(:,3),12,[1 0 0],'filled');
                view([90 -90]);
                
                
            end
        end
    end
end
                for ii=1:length(ind)
                    dataL=data(ss(sessioni),typei,align_typei,ind(ii));
                    for kk=1:length(dataL{1})
                        dataL2=cat(1,dataL{1}{kk},dataL2);
                    end
                end
                kk=1;
                
                nexttile(2+kk,[2,1]);kk=kk+1;
                scatter3(temp_supervoxel(:,1),temp_supervoxel(:,2),temp_supervoxel(:,3),8,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1);hold on;
                scatter3(dataL2(:,1),dataL2(:,2),dataL2(:,3),12,[1 0 0],'filled');
                view([90 -90]);
                
                
            end
        end
    end
end