%clustring of hab最后一个trial，acq block第一个trial，test第一个trial
set(0,'defaultfigurecolor','w');

%%%%all fish
load('E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\path.mat')%correct actpath
name='activities_aft_process';
for ii=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{ii});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{ii}=fullfile(fpath,[fname,ext]);
    end
end
trace=struct;trace.hab=[];index_all=[];
for ii=1:length(actbatch)
    load(parabatch{ii});
    load(actbatch{ii});
    index_all=[index_all;[ones(size(activities_preCS_dfdf_aftcorrect,2),1)*ii (1:size(activities_preCS_dfdf_aftcorrect,2))']];
    %每个trial最后一个block和后一个trial第一个block
    frame_ind=[1,frame.per_cycle];%%%%%%%%%%%%%%%%%%%!!!!!!!
    ind=(trial.hab(3)-1)*frame.per_cycle+frame_ind(1):(trial.hab(3)-1)*frame.per_cycle+frame_ind(2);
    trace.hab=[trace.hab activities_preCS_dfdf_aftcorrect(ind,:)];
end
result.trace_hab_end=trace;
trace=struct;trace.hab=[];trace.acq=[];trace.test=[];index_all=[];
for ii=1:length(actbatch)
    load(parabatch{ii});
    load(actbatch{ii});
    index_all=[index_all;[ones(size(activities_preCS_dfdf_aftcorrect,2),1)*ii (1:size(activities_preCS_dfdf_aftcorrect,2))']];
    %每个trial最后一个block和后一个trial第一个block
    frame_ind=[frame.cs_start,frame.us_start-1];%%%%%%%%%%%%%%%%%%%!!!!!!!
    ind=(trial.hab(2)-1)*frame.per_cycle+frame_ind(1):(trial.hab(2)-1)*frame.per_cycle+frame_ind(2);
    trace.hab=[trace.hab activities_preCS_dfdf_aftcorrect(ind,:)];
    ind=(trial.test(2)-1)*frame.per_cycle+frame_ind(1):(trial.test(2)-1)*frame.per_cycle+frame_ind(2);
    trace.test=[trace.test activities_preCS_dfdf_aftcorrect(ind,:)];
    indd=trial.acq(2):trial.acq_block_trial:trial.acq(3);ind=[];
    for ii=1:length(indd)
        ind=[ind (indd(ii)-1)*frame.per_cycle+frame_ind(1):(indd(ii)-1)*frame.per_cycle+frame_ind(2)];
    end
    trace.acq=[trace.acq activities_preCS_dfdf_aftcorrect(ind,:)];
end
result.trace=trace;%trace.hab(:,1)=[];trace.acq(:,1)=[];trace.test(:,1)=[];

trace=struct;trace.hab=[];trace.acq=[];trace.test=[];index_all=[];
for ii=1:length(actbatch)
    load(parabatch{ii});
    load(actbatch{ii});
    index_all=[index_all;[ones(size(activities_preCS_dfdf_aftcorrect,2),1)*ii (1:size(activities_preCS_dfdf_aftcorrect,2))']];
    %每个trial最后一个block和后一个trial第一个block
    frame_ind=[1,frame.per_cycle];%%%%%%%%%%%%%%%%%%%!!!!!!!
    ind=(trial.hab(2)-1)*frame.per_cycle+frame_ind(1):(trial.hab(2)-1)*frame.per_cycle+frame_ind(2);
    trace.hab=[trace.hab activities_preCS_dfdf_aftcorrect(ind,:)];
    ind=(trial.test(2)-1)*frame.per_cycle+frame_ind(1):(trial.test(2)-1)*frame.per_cycle+frame_ind(2);
    trace.test=[trace.test activities_preCS_dfdf_aftcorrect(ind,:)];
    indd=trial.acq(2):trial.acq_block_trial:trial.acq(3);ind=[];
    for ii=1:length(indd)
        ind=[ind (indd(ii)-1)*frame.per_cycle+frame_ind(1):(indd(ii)-1)*frame.per_cycle+frame_ind(2)];
    end
    trace.acq=[trace.acq activities_preCS_dfdf_aftcorrect(ind,:)];
end
result.trace_all=trace;
%acq block间变化
trace=struct;trace.hab=[];trace.acq=[];trace.test=[];
for ii=1:length(actbatch)
    load(parabatch{ii});
    load(actbatch{ii});
    %每个trial最后一个block和后一个trial第一个block
    frame_ind=[1 frame.per_cycle];%%%%%%%%%%%%%%%%%%%!!!!!!!
    ind=(trial.hab(2)-1)*frame.per_cycle+frame_ind(1):(trial.hab(2)-1)*frame.per_cycle+frame_ind(2);
    trace.hab=[trace.hab activities_preCS_dfdf_aftcorrect(ind,:)];
    ind=(trial.test(2)-1)*frame.per_cycle+frame_ind(1):(trial.test(2)-1)*frame.per_cycle+frame_ind(2);
    trace.test=[trace.test activities_preCS_dfdf_aftcorrect(ind,:)];
    indd=trial.acq(2):trial.acq(3);ind=[];
    for ii=1:length(indd)
        ind=[ind (indd(ii)-1)*frame.per_cycle+frame_ind(1):(indd(ii)-1)*frame.per_cycle+frame_ind(2)];
    end
    trace.acq=[trace.acq activities_preCS_dfdf_aftcorrect(ind,:)];
end
result.index_all=index_all;
result.trace_all_block=trace;

%% classification
trace=result.trace;
trace.hab(:,find(index_all(:,1)==9))=[];trace.acq(:,find(index_all(:,1)==9))=[];trace.test(:,find(index_all(:,1)==9))=[];
addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
M_0=cat(1,trace.hab,trace.acq,trace.test)';
M_norm=normalize(M_0,2,'zscore');M_0=M_norm;
%%%clustering
numK=[];
ind=1:length(M_0);
ind=ind([1:50:length(ind)]);
[gIX,~,numK,~]=find_best_k_in_range(M_0(ind,:),7:4:80);
cIX=ind;cIX_reg = ind;%(1:size(M_0,1))';
if isempty(numK)
    numK=3;
end
numK
masterthres=0.8;
para=struct;
para=setfield(para,'k1',numK);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);para=setfield(para,'minSize',20);
%cIX=(1:size(M_0,1))';gIX=(1:size(M_0,1))'; cIX_reg = (1:size(M_0,1))';
for ii=1:10
    %[cIX,gIX] = AutoClustering(randi(size(M_0,1),[1,50000]),[1:50000],M_0,cIX_reg,1,para,1,0.7);
    [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
    cIX_reg = (1:size(M_0,1))';
end
result.cIX=cIX;result.gIX=gIX;%!!!!!!!!!!!!!!!!!!!!!!!!
length(unique(gIX))
M_1=[];
for ii=1:length(unique(gIX))
    M_1(ii,:)=mean(M_0(cIX(find(gIX==ii)),:),1);
end
[gIX2,~,numK2,~]=find_best_k_in_range(M_1,1:2:length(unique(gIX)));
cIX2=1:length(unique(gIX));para=setfield(para,'k1',numK2);
para=setfield(para,'minSize',1);
for ii=1:5
    [cIX2,gIX2] = AutoClustering(cIX2,gIX2,M_1,cIX2,1,para,1,masterthres);%[cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
end
length(unique(gIX2))
result.cIX2=cIX2;result.gIX2=gIX2;%!!!!!!!!!!!!!!!!!!!!!!!!
gIX_rep=gIX;ind_rep=[];
for ii=1:length(unique(gIX2))
    ind=find(gIX2==ii);
    for jj=1:length(ind)
        ind_rep=[ind_rep;find(gIX==ind(jj))];
        gIX_rep(find(gIX==ind(jj)))=ii;
    end
end
ind_cut=setdiff(1:length(cIX),ind_rep);cIX_rep=cIX;
gIX_rep(ind_cut)=[];cIX_rep(ind_cut)=[];
length(unique(gIX_rep))
result.cIX_rep=cIX_rep;result.gIX_rep=gIX_rep;%!!!!!!!!!!!!!!!!!!!!!!!!
save('E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\results_20190731mat','result','-v7.3');

%% some plot
cIX=result.cIX_rep;gIX=result.gIX_rep;%!!!!!!!!!!!!!!!!!!!!!!!!

k = get(0, 'screensize');
clrmap_name = 'hsv_new';%getappdata(hfig,'clrmap_name');
clrmap = GetColormap(clrmap_name,max(gIX));
ind_control=setdiff(1:length(index_all),cIX)';%length(ind)+length(cIX)== length(index_all)
clr_control=[0,0,0];
cIX_all=[cIX;ind_control];
gIX_all=[gIX;(max(gIX)+1)*ones(size(ind_control))];
clrmap_all=cat(1,clrmap,clr_control);
length(unique(gIX))
trace=result.trace_all;
%cIX(35670:35679)=[];gIX(35670:35679)=[];
vi='off';
colorCS=[0.5 0.5 0.9];
for iii=1
    switch iii
        case 1
            cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;
        case 2
            cIX_iii=cIX_all;gIX_iii=gIX_all;clrmap_iii=clrmap_all;
    end
    %outputpath=checkpath([savepath '\figures_masterthres' num2str(masterthres) '_' num2str(iii,'%02d')]);
    for jj=1:3
        figure,
        a=[-0.05 0.15];
        switch jj
            case 1
                x=trace.hab;isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
            case 2
                x=trace.acq;isus=true;set(gcf,'position',[50 50 0.6*k(3) 0.85*k(4)]);
            case 3
                x=trace.test;isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
        end
        for ii=unique(gIX_iii)'
            ind=find(gIX_iii==ii);
            step=(1/(length(unique(gIX_iii))+1)*1)/(length(unique(gIX_iii))-1);
            high=1/(length(unique(gIX_iii))+3);
            ax=subplot('position',[0.1 0.05+(ii-1)*(1.5*step+high) 0.8 high]);%length(unique(gIX_iii)),1,ii,
            patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
                [frame.cs_end:frame.per_cycle:size(x,1)]'...
                [frame.cs_end:frame.per_cycle:size(x,1)]'...
                [frame.cs_start:frame.per_cycle:size(x,1)]']',...
                repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
                colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
            plot(x(:,cIX_iii(ind)),'linewidth',0.5,'color',[0.5 0.5 0.5]);hold on
            if isus
                l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
            end
            plot(mean(x(:,cIX_iii(ind)),2),'linewidth',2,'color',clrmap(ii,:,:));hold on
            ylim(a);box on;
        end
    end
end

%% check
num=11;%11;16;17
trace=result.trace_all;
loc=result.cIX_rep(find(result.gIX_rep==num));
figure,
a=[-0.05 0.15];trace_seed=[];
for ii=1:3
   
    switch ii
        case 1
            subplot(1,7,1), x=trace.hab(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
        case 2
            subplot(1,7,2:6),x=trace.acq(:,loc);isus=true;set(gcf,'position',[50 50 0.6*k(3) 0.85*k(4)]);
        case 3
            subplot(1,7,7),x=trace.test(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
    end
    patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_start:frame.per_cycle:size(x,1)]']',...
        repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    %plot(x(:,cIX_iii(ind)),'linewidth',0.5,'color',[0.5 0.5 0.5]);hold on
    if isus
        l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
    end
    plot(mean(x,2),'linewidth',2,'color',clrmap(num,:,:));hold on
    ylim(a);box on;
    trace_seed=[trace_seed;mean(x,2)];
end
trace=result.trace;
trace_seed=[];
for ii=1:3
    switch ii
        case 1
            x=trace.hab(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
        case 2
            x=trace.acq(:,loc);isus=true;set(gcf,'position',[50 50 0.6*k(3) 0.85*k(4)]);
        case 3
            x=trace.test(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
    end
    trace_seed=[trace_seed;mean(x,2)];
end
%top corr
M=[result.trace.hab;result.trace.acq;result.trace.test];
nCells_total = size(M,2);prct_const = 2;topN = round(prct_const/100 * nCells_total);%para
Corr=corr(trace_seed,M);
[~,IX] = sort(Corr,'descend');
x0 = Corr(IX(topN));
figure,hist(Corr);hold on;
line([x0 x0],[0 20000],'linewidth',2,'color','r');
IX_passX = find(Corr>=x0);
inter=length(intersect(result.cIX_rep(find(result.gIX_rep==num)),IX_passX))
dif=length(setdiff(result.cIX_rep(find(result.gIX_rep==num)),IX_passX))
IX_passX=union(result.cIX_rep(find(result.gIX_rep==num)),IX_passX);IX_passX(find(index_all(IX_passX,1)==9))=[];
loc=IX_passX;trace=result.trace_all;
number_IX_passX=length(IX_passX)
figure,
a=[-0.05 0.15];
area_win_acq=frame.cs_start:frame.us_start-1;%frame.us_start
type_cal_area='event';%type={'lowest','lowest_all','polyfit','first','mean','sum','event'};
ref_win=ceil(4.8/fs.ca+1):frame.cs_start-1;%取取时间段的baseline作为参考判断event
area_win_hab_tst=frame.cs_start:frame.cs_end-1;
area=struct;area_sum=struct;
for ii=1:3
    switch ii
        case 1
            subplot(1,7,1),
            x=trace.hab(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
            area.hab=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
            area_sum.hab(:,1)=mean(area.hab);area_sum.hab(:,2)=std(area.hab)/sqrt(length(loc));
        case 2
            subplot(1,7,2:6),
            x=trace.acq(:,loc);isus=true;set(gcf,'position',[50 50 0.6*k(3) 0.85*k(4)]);
            y=reshape(x,frame.per_cycle,5,[]);%y(:,:,1);
            for j=1:5;size(y,3)
                xx=reshape(y(:,j,:),frame.per_cycle,[]);
                area.acq{j,1}=calculate_integrate_dfdf(xx,area_win_acq,type_cal_area,ref_win);
                area_sum.acq(j,1)=mean(area.acq{j,1});area_sum.acq(j,2)=std(area.acq{j,1})/sqrt(length(loc));
            end
        case 3
            subplot(1,7,7),
            x=trace.test(:,loc);isus=false;set(gcf,'position',[50 50 0.2*k(3) 0.85*k(4)]);
            area.test=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
            area_sum.test(:,1)=mean(area.test);area_sum.test(:,2)=std(area.test)/sqrt(length(loc));
    end
    patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_end:frame.per_cycle:size(x,1)]'...
        [frame.cs_start:frame.per_cycle:size(x,1)]']',...
        repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
        colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    %plot(x(:,cIX_iii(ind)),'linewidth',0.5,'color',[0.5 0.5 0.5]);hold on
    if isus
        l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
    end
    plot(mean(x,2),'linewidth',2,'color',clrmap(num,:,:));hold on
    ylim(a);xlim([1 length(mean(x,2))]);box on;
end
figure,
a=[-0.02 0.06];
x=result.trace_all_block.acq(:,IX_passX);
patch1=patch([[frame.cs_start:frame.per_cycle:size(x,1)]'...
    [frame.cs_end:frame.per_cycle:size(x,1)]'...
    [frame.cs_end:frame.per_cycle:size(x,1)]'...
    [frame.cs_start:frame.per_cycle:size(x,1)]']',...
    repmat([min(a) min(a) max(a) max(a)],length([frame.cs_end:frame.per_cycle:size(x,1)]),1)',...
    colorCS,'edgecolor',colorCS,'facecolor',colorCS,'edgealpha',0.2,'facealpha',0.25);hold on
    l3=plot(repmat([(frame.us_start):frame.per_cycle:size(x,1)],2,1)',[min(a) max(a)],'r','LineStyle','--','linewidth',1.2);hold on
    l3=plot(repmat([1:frame.per_cycle*trial.acq_block_trial:size(x,1)],2,1)',[min(a) max(a)],'b','LineStyle','-','linewidth',1.2);hold on
plot(mean(x,2),'linewidth',2,'color',clrmap(num,:,:));hold on
ylim(a);box on;

area_mean=struct;area_sum_mean=struct;
for ii=1:3
    switch ii
        case 1
            x=result.trace_all.hab(:,loc);
            area_mean.hab=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
            area_sum_mean.hab(:,1)=mean(area_mean.hab);area_sum_mean.hab(:,2)=std(area_mean.hab)/sqrt(length(loc));
        case 2
            x=result.trace_all_block.acq(:,loc);
            y=reshape(x,frame.per_cycle,trial.acq_block_trial,trial.acq_block_num,[]);%y(:,:,:,1);
            y=mean(y,2);y=reshape(y,frame.per_cycle,trial.acq_block_num,[]);
            for j=1:5;size(y,3)
                xx=reshape(y(:,j,:),frame.per_cycle,[]);
                area_mean.acq{j,1}=calculate_integrate_dfdf(xx,area_win_acq,type_cal_area,ref_win);
                area_sum_mean.acq(j,1)=mean(area_mean.acq{j,1});area_sum_mean.acq(j,2)=std(area_mean.acq{j,1})/sqrt(length(loc));
            end
        case 3
            x=result.trace_all.test(:,loc);
            area_mean.test=calculate_integrate_dfdf(x,area_win_acq,type_cal_area,ref_win);
            area_sum_mean.test(:,1)=mean(area_mean.test);area_sum_mean.test(:,2)=std(area_mean.test)/sqrt(length(loc));
    end
end
[h,p] = ttest(area_mean.acq{3},area_mean.acq{5})
%% location
M=[result.trace_all_block.hab;result.trace_all_block.acq;result.trace_all_block.test];
for ii=1:2
    switch ii
        case 1
            loc=index_all(IX_passX,:);
            figure,hist(loc(:,1),[1:8]);
%             for kk=unique(loc(:,1))'
%             figure,plot(mean(M(:,IX_passX(find(loc(:,1)==kk))),2));
%             end
        case 2
            loc=index_all(result.cIX_rep(find(result.gIX_rep==num)),:);
    end
    figure,
    subplot(1,2,1);hist(loc(:,1),unique(index_all(:,1)));
    [counts,C]=hist(loc(:,1),unique(index_all(:,1)));
    subplot(1,2,2);hist(index_all(:,1),unique(index_all(:,1)));
    [counts_all,C]=hist(index_all(:,1),unique(index_all(:,1)));
    figure, bar(C,counts./counts_all);
    for kk=unique(loc(:,1))'
        ind_event_in_this_fish=loc(find(loc(:,1)==kk),:);
        load(envbatch{kk})
        %         load(actbatch{kk})
        %          figure,
        %          plot(mean(activities_preCS_dfdf_aftcorrect(:,ind_event_in_this_fish(:,2)),2))
        h=DrawTiledPics_zyq_20190530(ind_event_in_this_fish(:,2),num*ones(size(ind_event_in_this_fish,1),1),[1:size(env.supervoxel,1)],[env.supervoxel(:,2) env.supervoxel(:,1) env.supervoxel(:,3)],env.vol,clrmap);
    end
end


