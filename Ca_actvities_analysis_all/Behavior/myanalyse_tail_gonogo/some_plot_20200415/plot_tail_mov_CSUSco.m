%updata from F:\DUlab\FC
%analyse\Ca_actvities_analyze\Behavior\myanalyse_tail_gonogo\some_plot_20200415\plot_tail_mov_US_2_1.m
%updata: add CS,for CS_US pair
%zyq,20200521

clc;clear all;close all;
set(0,'defaultfigurecolor','w');

[name path]=uigetfile('I:\1.Test US\2.Tail free——Data from 117\*.mat','MultiSelect', 'on','alltailpos,movement detection');
if ~iscell(name)==1
    load([path name]);
else
    for ii=1:length(name)
        load([path name{ii}]);
    end
end
disp(path)
%% alltheta
controlpoint=8;basepoint=1;%%%%！！！！计算尾动角度的点
nframes=length(movlabel);
basevec = [selx(controlpoint) sely(controlpoint)] - [selx(basepoint) sely(basepoint)];
movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(basepoint); sely(basepoint)], 1, size(alltailpos,3));
% distance of control points
dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));
% theta, anti-clockwise as positive
theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
theta(theta > 180) = theta(theta > 180) - 360;
theta = -theta;
theta00 = 0;
alltheta_r = theta00 * ones(nframes, 1);
alltheta_r(movind) = theta;

alltheta=alltheta_r;%alltheta(end+1:244800)=0;alltheta(244801:end)=[];
alltheta(frameb.per_cycle*trial.total+1:end)=[];
%%%%%%%%%%如果错帧，跑到这里
cut_thr=[0.5 5];%%%%！！！！矫正参数
aa=abs(diff([alltheta(1);alltheta]));alltheta(find(aa<cut_thr(1) & abs(alltheta)<cut_thr(2)))=0;

alltheta_n=normalize(alltheta,'zscore');
d_alltheta=diff([alltheta(1);alltheta]);
d_alltheta_n=diff([alltheta_n(1);alltheta_n]);
%% setpara
[fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=0;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num)
%plot_behavior_onset(a.delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
%spon
T_spon={};T_US={};l_spon={};l_US={};T_non_spon=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.spon_bef(2)-1)*frameb.per_cycle+1:(trial.spon_bef(3))*frameb.per_cycle;
            l_spon{ii}='Bef.';
        case mat2cell([2:trial.acq_block_num],1,ones(1,trial.acq_block_num-1))
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)]; 
        case trial.acq_block_num+1
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
            %l_spon{ii}='Test';
            l_spon{ii}=['Acq. Interv.' num2str(ii-1)]; 
        case trial.acq_block_num+2
            t=(trial.test(3))*frameb.per_cycle+1:(trial.total)*frameb.per_cycle;
            l_spon{ii}='Aft.';
    end
    T_spon{ii}=t;
end
% us
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frameb.per_cycle+1:(trial.hab(3))*frameb.per_cycle;
            l_US{ii}='Hab';
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle ;
            
            l_US{ii}=['Acq.' num2str(ii-1)];
        case trial.acq_block_num+2
            t=(trial.test(2)-1)*frameb.per_cycle+1:(trial.test(3))*frameb.per_cycle;
            l_US{ii}='Tst';
    end
     T_US{ii}=t;
     T_non_spon=[T_non_spon,t];
end
stimUS=zeros(1,frameb.per_cycle*trial.total);ind_US=[];kk=1;
for ss=1:trial.acq_block_num
    for tt=1:trial.acq_block_trial
        ind=(trial.hab(3)+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frameb.per_cycle+frameb.us_start;
        stimUS(ind:ind+max(round(frameb.us_dur(ss)),1)-1)=1;%23
        ind_US(kk,1)=ind;kk=kk+1;
    end
end
trial_ind_CS=[trial.hab(2):trial.hab(3)];
for ii=1:trial.acq_block_num
    trial_ind_CS=[trial_ind_CS,[trial.hab(3)+1:trial.hab(3)+trial.acq_block_trial]+(ii-1)*(trial.acq_block_interval+trial.acq_block_trial)];
end
trial_ind_CS=[trial_ind_CS,trial.test(2):trial.test(3)];
stimCS=zeros(1,frameb.per_cycle*trial.total);ind_CS=[];kk=1;
for tt=trial_ind_CS
    stimCS((tt-1)*frameb.per_cycle+frameb.cs_start:(tt-1)*frameb.per_cycle+frameb.cs_end)=1;%23
    ind_CS(kk,1)=(tt-1)*frameb.per_cycle+frameb.cs_start;
    ind_CS(kk,2)=(tt-1)*frameb.per_cycle+frameb.cs_end;kk=kk+1;
end
%figure,plot(stimCS,'r');hold on;plot(stimUS,'k');
%% plot
intermovinterv_spon=3*60;
intermovinterv_US=0.5*60;
merge_ind_US=0.5*60;

a=ceil(max(alltheta)/5)*5*[-1 1];T=[1:trial.total*frameb.per_cycle]/60/60;
figure('position',[50,250,1700,450]),plot(T,stimUS*180,'k','linewidth',1.5);hold on;plot(T,stimCS*180,'r','linewidth',1.5);hold on;
plot(T,alltheta,'linewidth',1);xlim([min(T) max(T)]);ylim(a);xlabel('Time(min)','fontsize',16);ylabel('Tail Angle');
set(gca,'fontsize',20,'linewidth',1.5);
%spon
Results_of_alltheta=struct;
a=50;max(alltheta);
a=a*[-1 1];
Results_of_alltheta.Theta_avg_spon=[];Results_of_alltheta.Theta_sum_spon=[];
Results_of_alltheta.Theta_avg_spon_n=[];Results_of_alltheta.Theta_sum_spon_n=[];
Results_of_alltheta.avg_env_dur_spon=[];Results_of_alltheta.avg_inter_env_interv_spon=[];Results_of_alltheta.M_num_spon=[];
Results_of_alltheta.spon_env_loc=nan(trial.acq_block_num+1,100,2);
figure('position',[50,1,60*trial.hab(1),120*trial.acq_block_num+2]);
for ii=1:trial.acq_block_num+2
     t=T_spon{ii};
    subaxis(trial.acq_block_num+2,1,ii,'SpacingVertical',0.04,'SpacingHorizontal',0.01);
    plot(T(t),stimUS(t)*180,'linewidth',1.5);hold on;plot(T(t),alltheta(t),'b','linewidth',1);
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
%us
% figure('position',[50,1,200*trial.acq_block_trial,120*trial.acq_block_num+2]);
% for ii=1:trial.acq_block_num+2
%     t=T_US{ii};
%     subaxis(trial.acq_block_num+2,1,ii,'SpacingVertical',0.05,'SpacingHorizontal',0.01);
%     %subplot(trial.acq_block_num,1,ii);%subplot(ceil(trial.acq_block_num/2),2,ii);
%     plot(T(t),stimCS(t)*180,'k','linewidth',1);hold on;plot(T(t),stimUS(t)*180,'r-','linewidth',1);hold on;
%     plot(T(t),alltheta(t),'b','linewidth',1);hold on;    %plot(T(t),abs([d_alltheta(t)]));hold on;
%     xlim([min(T(t)) max(T(t))]);ylim(a);title(l_US{ii},'fontsize',15);box on;
%     xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',14,'linewidth',1.2);
% end
Results_of_alltheta.Theta_avg_US=[];Results_of_alltheta.Theta_sum_US=[];
Results_of_alltheta.Theta_avg_US_n=[];Results_of_alltheta.Theta_sum_US_n=[];
Results_of_alltheta.avg_env_dur_us=[];Results_of_alltheta.avg_inter_env_interv_us=[];Results_of_alltheta.M_num_US=[];
Results_of_alltheta.Theta_US=nan(trial.acq_block_num+2,trial.acq_block_trial);Results_of_alltheta.Theta_US_n=nan(trial.acq_block_num+2,trial.acq_block_trial);Results_of_alltheta.env_dur_us=nan(trial.acq_block_num+2,trial.acq_block_trial);
Results_of_alltheta.env_loc_us=nan(trial.acq_block_num+2,trial.acq_block_trial,2);
Results_of_alltheta.M_num_US_spon=[];
figure('position',[50,1,200*trial.acq_block_trial,122*trial.acq_block_num+2]);
for ii=1:trial.acq_block_num+2
    t=T_US{ii};
    subaxis(trial.acq_block_num+2,1,ii,'SpacingVertical',0.05,'SpacingHorizontal',0.01);
    %subplot(trial.acq_block_num,1,ii);%subplot(ceil(trial.acq_block_num/2),2,ii);
    plot(T(t),stimUS(t)*180,'r-','linewidth',1.5);hold on;
    plot(T(t),stimCS(t)*180,'k','linewidth',1.5);hold on;
    plot(T(t),alltheta(t),'b','linewidth',1);hold on;    %plot(T(t),abs([d_alltheta(t)]));hold on;
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
        %cs
        aa=(repmat(t(env_loc),length(ind_CS(:,1)),1)-repmat(ind_CS(:,1),1,length(env_loc)));
        ind=find(aa>=-merge_ind_US & abs(aa)<=(frameb.cs_end-frameb.cs_start)+1);[~,ind]=ind2sub([size(aa,1),size(aa,2)],ind);        
        env_loc_cs=env_loc(ind);env_num=length(ind);
        scatter(T(t(env_loc_cs)),alltheta(t(env_loc_cs)),12,'r','filled');hold on;
        %us
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
        Results_of_alltheta.env_loc_us(ii,idx,1)=t(env_loc); Results_of_alltheta.env_loc_us(ii,idx,2)=t(env_end);
        else
        Results_of_alltheta.Theta_US(ii,:)=nan(1,trial.acq_block_trial);
        Results_of_alltheta.Theta_US_n(ii,:)=nan(1,trial.acq_block_trial);
        Results_of_alltheta.env_dur_us(ii,:)=nan(1,trial.acq_block_trial);
        Results_of_alltheta.env_loc_us(ii,:,1)=nan(1,trial.acq_block_trial);
        Results_of_alltheta.env_loc_us(ii,:,2)=nan(1,trial.acq_block_trial);
        end
    else
        Results_of_alltheta.Theta_avg_US(ii)=nan;
        Results_of_alltheta.Theta_avg_US_n(ii)=nan;
        Results_of_alltheta.avg_env_dur_us(ii)=nan;
        Results_of_alltheta.M_num_US(ii)=0;
        Results_of_alltheta.Theta_US(ii,:)=nan;
        Results_of_alltheta.Theta_US_n(ii,:)=nan;
         Results_of_alltheta.env_loc_us(ii,:,:)=nan;
    end
    %figure,plot(alltheta(tt));hold on;plot(abs(d_alltheta(tt)));
    xlim([min(T(t)) max(T(t))]);ylim(a);title(l_US{ii},'fontsize',15);box on;
    xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',14,'linewidth',1.2);
end

h=figure;
[env_num,env.env_loc,env.env_end,~]=env_detect(alltheta(T_non_spon),intermovinterv_US,1);
plot_tailmovement(alltheta(T_non_spon),fs,frameb,trial,env.env_loc);
set(axes,'position',[0.695,0.4,0.12,0.20]);box on;
bar(categorical({l_US{2:end-1}}),Results_of_alltheta.M_num_US(2:end-1));
title('UR. Num');set(gca,'fontsize',14);ylim([0 trial.acq_block_trial]);
set(axes,'position',[0.85,0.4,0.12,0.20]);box on;
bar(categorical({l_US{2:end-1}}),Results_of_alltheta.Theta_avg_US(2:end-1));hold on;
scatter(reshape(categorical(repmat({l_US{2:end-1}},3,1)'),[],1),reshape(Results_of_alltheta.Theta_US(2:end-1,:),[],1),20,'k','filled');
title('Avg. Theta');set(gca,'fontsize',14);
set(axes,'position',[0.695,0.08,0.12,0.20]);box on;
bar(categorical({l_US{2:end-1}}),Results_of_alltheta.avg_env_dur_us(2:end-1));hold on;
scatter(reshape(categorical(repmat({l_US{2:end-1}},3,1)'),[],1),reshape(Results_of_alltheta.env_dur_us(2:end-1,:),[],1),20,'k','filled');
title('Avg. UR Dur(s)');set(gca,'fontsize',14);
%beeswarm(reshape(categorical(repmat(l_us,3,1)'),[],1),reshape(Results_of_alltheta.Theta_US,[],1),'dot_size',.5,'overlay_style','sd');
set(axes,'position',[0.9,0.08,0.1,0.02]);box on;text(0,0,path(end-14:end),'fontsize',14);
set(gca,'visible','off','yticklabel',[],'xticklabel',[])

%statistic
figure('position',[50 1 260*2 245*3]),
subplot(3,2,1),
bar(categorical(l_spon),Results_of_alltheta.M_num_spon,'XDataMode','manual','xdata',[1:length(l_spon)]);
title('Mov. Rate/Min');set(gca,'xticklabels',l_spon,'XTickLabelRotation',45,'fontsize',12);
subplot(3,2,3),
bar(categorical(l_spon), Results_of_alltheta.avg_env_dur_spon,'XDataMode','manual','xdata',[1:length(l_spon)]);
title('Avg. Mov Dur(s)');set(gca,'xticklabels',l_spon,'XTickLabelRotation',45,'fontsize',12);
subplot(3,2,5),
bar(categorical(l_spon),Results_of_alltheta.avg_inter_env_interv_spon,'XDataMode','manual','xdata',[1:length(l_spon)]);
title('Avg. Inter Mov Interv(s)');set(gca,'xticklabels',l_spon,'XTickLabelRotation',45,'fontsize',12);
subplot(3,2,2),
bar(categorical(l_spon),Results_of_alltheta.Theta_avg_spon,'XDataMode','manual','xdata',[1:length(l_spon)]);
title('Avg. Theta');set(gca,'xticklabels',l_spon,'XTickLabelRotation',45,'fontsize',12);
subplot(3,2,4),
bar(categorical(l_spon),Results_of_alltheta.Theta_sum_spon,'XDataMode','manual','xdata',[1:length(l_spon)]);
title('Sum Theta');set(gca,'xticklabels',l_spon,'XTickLabelRotation',45,'fontsize',12);


save([path 'Results_of_alltheta.mat'],'Results_of_alltheta','env','controlpoint','alltheta','intermovinterv_spon','intermovinterv_US','merge_ind_US','cut_thr');


