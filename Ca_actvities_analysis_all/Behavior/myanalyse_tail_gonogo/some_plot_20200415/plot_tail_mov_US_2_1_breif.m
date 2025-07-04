%updata from F:\DUlab\FC
%analyse\Ca_actvities_analyze\Behavior\myanalyse_tail_gonogo\some_plot_20200415\plot_tail_mov_US.m
%updata: setpara_testUS & getUSstim
%zyq,20200430

clc;clear all;close all;
set(0,'defaultfigurecolor','w');

[name path]=uigetfile('H:\1.Test US\2.Tail free����Data from 117\*.mat','MultiSelect', 'on','alltailpos,movement detection;US response');
if ~iscell(name)==1
    load([path name]);
else
    for ii=1:length(name)
        load([path name{ii}]);
    end
end
disp(path)
%% setpara
[frame,trial]=setpara_testUS(3);
[stimUS, ind_US,T_spon,T_US,l_spon,l_us]=getUSstim(trial,frame);%figure,plot( stimUS);
% figure('position',[50,50,1750,950]),
% for ii=1:length(ind_US)
%     subaxis(trial.acq_block_trial,trial.acq_block_number,ii,'SpacingVertical',0.01,'SpacingHorizontal',0.01);
%     imshow(US_I_m(:,:,ii),[min(min(US_I_m(:,:,ii))) max(max(US_I_m(:,:,ii)))]);
% end
%% alltheta
alltheta_n=normalize(alltheta,'zscore');
d_alltheta=diff([alltheta(1);alltheta]);
d_alltheta_n=diff([alltheta_n(1);alltheta_n]);

%% plot
T=[1:trial.total*frame.per_cycle];
% avg angle
Results_of_alltheta=struct;
%spon
a=50;max(alltheta);
a=a*[-1 1];
Results_of_alltheta.Theta_avg_spon=[];Results_of_alltheta.Theta_sum_spon=[];
Results_of_alltheta.Theta_avg_spon_n=[];Results_of_alltheta.Theta_sum_spon_n=[];
Results_of_alltheta.avg_env_dur_spon=[];Results_of_alltheta.avg_inter_env_interv_spon=[];Results_of_alltheta.M_num_spon=[];
Results_of_alltheta.spon_env_loc=nan(trial.acq_block_number+1,100,2);
figure('position',[50,1,60*trial.hab,125*trial.acq_block_number+1]);
for ii=1:trial.acq_block_number+1
    t=T_spon{ii};
    %subplot(trial.acq_block_number+1,1,ii);
    subaxis(trial.acq_block_number+1,1,ii,'SpacingVertical',0.04,'SpacingHorizontal',0.01);
    plot(T(t),stimUS(t)*180,'r','linewidth',1.5);hold on;plot(T(t),alltheta(t),'linewidth',1);
    [env_num,env_loc,env_end,~]=env_detect(alltheta(t),intermovinterv_spon,1);
    scatter(T(t(env_loc)),alltheta(t(env_loc)),10,'r','filled');
    scatter(T(t(env_end)),alltheta(t(env_end)),10,'k','filled');
    ss=[];ss_n=[];
    %figure,
    for jj=1:length(env_loc)
        ss_n(jj)=sum(abs(d_alltheta_n(t(env_loc(jj)):t(env_end(jj)))));
        ss(jj)=sum(abs(d_alltheta(t(env_loc(jj)):t(env_end(jj)))));
        %          subplot(length(env_loc),1,jj);plot(abs(d_alltheta(t(env_loc(jj)):t(env_end(jj)))));hold on
        %          plot(abs(alltheta(t(env_loc(jj)):t(env_end(jj)))));
    end
    Results_of_alltheta.Theta_avg_spon_n(ii)=mean(ss_n);%sum(abs(d_alltheta_n(t)),'omitnan')/env_num;
    Results_of_alltheta.Theta_sum_spon_n(ii)=sum(ss_n);%sum(abs(d_alltheta_n(t)),'omitnan')/(length(t)/60/60);
    Results_of_alltheta.Theta_avg_spon(ii)=mean(ss);%sum(abs(d_alltheta(t)),'omitnan')/env_num;
    Results_of_alltheta.Theta_sum_spon(ii)=sum(ss);%sum(abs(d_alltheta(t)),'omitnan')/(length(t)/60/60);
    Results_of_alltheta.M_num_spon(ii)=env_num/(length(t)/60/60);
    Results_of_alltheta.avg_env_dur_spon(ii)=mean(env_end-env_loc)/60;
    Results_of_alltheta.avg_inter_env_interv_spon(ii)=mean(diff(env_loc))/60;
    Results_of_alltheta.spon_env_loc(ii,1:length(env_loc),1)=t(env_loc);
    Results_of_alltheta.spon_env_loc(ii,1:length(env_loc),2)=t(env_end);
    SS{ii}=ss;
    xlim([min(T(t)) max(T(t))]);ylim(a);title(l_spon{ii},'fontsize',15);box on;
    xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',14,'linewidth',1.5);
end
% us
Results_of_alltheta.Theta_avg_US=[];Results_of_alltheta.Theta_sum_US=[];
Results_of_alltheta.Theta_avg_US_n=[];Results_of_alltheta.Theta_sum_US_n=[];
Results_of_alltheta.avg_env_dur_us=[];Results_of_alltheta.avg_inter_env_interv_us=[];Results_of_alltheta.M_num_US=[];
Results_of_alltheta.Theta_US=nan(trial.acq_block_number,trial.acq_block_trial);Results_of_alltheta.Theta_US_n=nan(trial.acq_block_number,trial.acq_block_trial);Results_of_alltheta.env_dur_us=nan(trial.acq_block_number,trial.acq_block_trial);
Results_of_alltheta.env_loc=nan(trial.acq_block_number,trial.acq_block_trial,2);
figure('position',[50,1,200*trial.acq_block_trial,120*trial.acq_block_number]);
for ii=1:trial.acq_block_number
    t=T_US{ii};
    subaxis(trial.acq_block_number,1,ii,'SpacingVertical',0.05,'SpacingHorizontal',0.01);
    %subplot(trial.acq_block_number,1,ii);%subplot(ceil(trial.acq_block_number/2),2,ii);
    plot(T(t),stimUS(t)*180,'r-','linewidth',1.5);hold on;
    plot(T(t),alltheta(t),'linewidth',1);hold on;    %plot(T(t),abs([d_alltheta(t)]));hold on;
    [env_num,env_loc,env_end,~]=env_detect(alltheta(t),intermovinterv_US,1);
    env_loc=max(env_loc-1,1);
    scatter(T(t(env_loc)),alltheta(t(env_loc)),12,'g','filled');hold on;
    ss=nan(trial.acq_block_trial,1);ss_n=nan(trial.acq_block_trial,1);
    for jj=1:length(env_loc)
        ss_n(jj,:)=sum(abs(d_alltheta_n(t(env_loc(jj)):t(env_end(jj)))));
        ss(jj,:)=sum(abs(d_alltheta(t(env_loc(jj)):t(env_end(jj)))));
    end
    
    Results_of_alltheta.Theta_sum_US(ii)=sum(ss);%sum(abs(d_alltheta(t)),'omitnan')/(length(t)/60/60);
    Results_of_alltheta.Theta_sum_US_n(ii)=sum(ss_n);%sum(abs(d_alltheta_n(t)),'omitnan')/(length(t)/60/60);
    %
    tt=[];
    if ~isempty(env_loc)
        scatter(T(t(env_end)),alltheta(t(env_end)),12,'k','filled');hold on;
        aa=min(abs(repmat(t(env_loc),length(ind_US),1)-repmat(ind_US,1,length(env_loc))),[],1);
        ind=find(aa<=merge_ind_US);
        env_loc=env_loc(ind);env_num=length(ind);env_end=env_end(ind);
        for jj=1:length(env_loc)
            tt=[tt,env_loc(jj):env_end(jj)];
            plot(T(t(env_loc(jj):env_end(jj))),alltheta(t(env_loc(jj):env_end(jj))),'color',[1 0.5 0],'linewidth',1.2);hold on;
            %plot(abs(d_alltheta(t(env_loc(jj):env_end(jj)))),'color',[1 0.5 0],'linewidth',1.2);hold on;
        end
        tt=t(tt);
        scatter(T(t(env_loc)),alltheta(t(env_loc)),12,'m','filled');hold on;
        Results_of_alltheta.Theta_avg_US(ii)=sum(ss(ind,:))/env_num;trial.acq_block_trial;%sum(abs(d_alltheta(tt)),'omitnan')/env_num;
        Results_of_alltheta.Theta_avg_US_n(ii)=sum(ss_n(ind,:))/env_num;trial.acq_block_trial;%sum(abs(d_alltheta_n(tt)),'omitnan')/env_num;
        Results_of_alltheta.avg_env_dur_us(ii)=mean(env_end-env_loc)/60;
        Results_of_alltheta.M_num_US(ii)=env_num;
        if ~isempty(env_loc)
            [~,idx]=min(abs(repmat(t(env_loc),length(ind_US),1)-repmat(ind_US,1,length(env_loc))),[],1);
            idx=mod(idx,trial.acq_block_trial);idx(find(idx==0))=trial.acq_block_trial;
            Results_of_alltheta.Theta_US(ii,idx)=ss(ind,:);
            Results_of_alltheta.Theta_US_n(ii,idx)=ss_n(ind,:);
            Results_of_alltheta.env_dur_us(ii,idx)=(env_end-env_loc)/60;
            Results_of_alltheta.env_loc(ii,idx,1)=t(env_loc); Results_of_alltheta.env_loc(ii,idx,2)=t(env_end);
        else
            Results_of_alltheta.Theta_US(ii,:)=nan(1,trial.acq_block_trial);
            Results_of_alltheta.Theta_US_n(ii,:)=nan(1,trial.acq_block_trial);
            Results_of_alltheta.env_dur_us(ii,:)=nan(1,trial.acq_block_trial);
            Results_of_alltheta.env_loc(ii,:,1)=nan(1,trial.acq_block_trial);
            Results_of_alltheta.env_loc(ii,:,2)=nan(1,trial.acq_block_trial);
        end
    else
        Results_of_alltheta.Theta_avg_US(ii)=nan;
        Results_of_alltheta.Theta_avg_US_n(ii)=nan;
        Results_of_alltheta.avg_env_dur_us(ii)=nan;
        Results_of_alltheta.M_num_US(ii)=0;
        Results_of_alltheta.Theta_US(ii,:)=nan;
        Results_of_alltheta.Theta_US_n(ii,:)=nan;
        Results_of_alltheta.env_loc(ii,:,:)=nan;
    end
    %figure,plot(alltheta(tt));hold on;plot(abs(d_alltheta(tt)));
    xlim([min(T(t)) max(T(t))]);ylim(a);title(['Session ' num2str(ii)],'fontsize',15);box on;
    xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',14,'linewidth',1.2);
end
%statistic
figure('position',[50 1 275*2 255*3]),
subplot(3,2,1),
bar(categorical(l_spon),Results_of_alltheta.M_num_spon);colormap('parula')
ylabel('Mov. Rate/Min');set(gca,'fontsize',14);
subplot(3,2,3),
bar(categorical(l_spon), Results_of_alltheta.avg_env_dur_spon);colormap('parula')
ylabel('Avg. Mov Dur(s)');set(gca,'fontsize',14);
subplot(3,2,5),
bar(categorical(l_spon),Results_of_alltheta.avg_inter_env_interv_spon);colormap('parula')
ylabel('Avg. Inter Mov Interv(s)');set(gca,'fontsize',14);
subplot(3,2,2),
bar(categorical(l_spon),Results_of_alltheta.Theta_avg_spon);colormap('parula')
ylabel('Avg. Theta');set(gca,'fontsize',14);
subplot(3,2,4),
bar(categorical(l_spon),Results_of_alltheta.Theta_sum_spon);colormap('parula')
ylabel('Sum Theta');set(gca,'fontsize',14);

figure('position',[50 50 275 255*3]),
subplot(3,1,1),
bar(categorical(l_us),Results_of_alltheta.M_num_US);
ylabel('UR. Num');set(gca,'fontsize',14);ylim([0 trial.acq_block_trial]);
% bar(categorical(l_us),Results_of_alltheta.Theta_sum_US);
% ylabel('Sum Theta');set(gca,'fontsize',14);
subplot(3,1,2),
bar(categorical(l_us),Results_of_alltheta.Theta_avg_US);hold on;
scatter(reshape(categorical(repmat(l_us,3,1)'),[],1),reshape(Results_of_alltheta.Theta_US,[],1),20,'k','filled');
ylabel('Avg. Theta');set(gca,'fontsize',14);
subplot(3,1,3),
bar(categorical(l_us),Results_of_alltheta.avg_env_dur_us);hold on;
scatter(reshape(categorical(repmat(l_us,3,1)'),[],1),reshape(Results_of_alltheta.env_dur_us,[],1),20,'k','filled');
ylabel('Avg. UR Dur(s)');set(gca,'fontsize',14);
%beeswarm(reshape(categorical(repmat(l_us,3,1)'),[],1),reshape(Results_of_alltheta.Theta_US,[],1),'dot_size',.5,'overlay_style','sd');

%save
path_summary='H:\1.Test US\2.Tail free����Data from 117\summary.mat';
load(path_summary);
for ii=2:length(Path)
    path=Path(ii);
    ind=find(strcmp(Path,path));
    if ~isempty(ind)
        Path(ind,:)=path;
    else
        Path(end+1,:)=path;
    end
    ind=find(strcmp(Path,path));
    load(fullfile(Path(ii),'Results_of_alltheta.mat'));
    name=fieldnames(Results_of_alltheta);
    for ii=1:length(name)
        if ~isfield(Data,name{ii})
            Data=setfield(Data,name{ii},[]);
        end
        if min(size(getfield(Results_of_alltheta,name{ii})))<=1
            b=getfield(Data,name{ii})';b(ind,1:length(getfield(Results_of_alltheta,name{ii})),:)=getfield(Results_of_alltheta,name{ii});
        elseif min(size(getfield(Results_of_alltheta,name{ii})))>1 && ~strcmp(name{ii},'spon_env_loc')
            b=getfield(Data,name{ii})';
            for jj=1:trial.acq_block_trial
                bb=getfield(Results_of_alltheta,name{ii});
                b{ind,jj}=bb(:,jj,:);
            end
        elseif min(size(getfield(Results_of_alltheta,name{ii})))>1 && strcmp(name{ii},'spon_env_loc')
            b=getfield(Data,name{ii})';
            for jj=1:size(getfield(Results_of_alltheta,name{ii}),1)
                bb=getfield(Results_of_alltheta,name{ii});
                b{ind,jj}=bb(jj,:,:);
            end
        end
        Data=setfield(Data,name{ii},b');
    end
end
save(path_summary,'Path','Data');
