clc;clear all;close all;
set(0,'defaultfigurecolor','w');

[name path]=uigetfile('H:\1.Test US\2.Tail free¡ª¡ªData from 117\*.mat','MultiSelect', 'on','alltailpos,movement detection;US response');
if ~iscell(name)==1
    load([path name]);
else
    for ii=1:length(name)
        load([path name{ii}]);
    end
end
disp(path)
%% US response
% [imgname, dirpath] = uigetfile('/mnt/W/data/fear conditioning_ZYQ/Test US/*.seq','*.seq');
% imgpath = [dirpath imgname]
% [i,~]=readtailseq(imgpath, 1);
% US_I_m=zeros(size(i,1),size(i,2),length(ind_US));
% H=figure;
% n=30;
% for ii=1:length(ind_US)
%     I=zeros(size(i,1),size(i,2),n);
%     for jj=1:n
%     [I(:,:,jj), ~] = readtailseq(imgpath, ind_US(ii)+jj-1);
%     end
%     US_I_m(:,:,ii)=max(I,[],3);
%     subplot(trial.acq_block_trial,trial.acq_block_number,ii),
%     imshow(US_I_m(:,:,ii),[min(min(US_I_m(:,:,ii))) max(max(US_I_m(:,:,ii)))]);
% end
% savepath = checkpath([dirpath '/' imgname '_find_tail/']);
% savefig(H, [savepath '/US RESPONSE','.fig']);
% save([savepath '/US RESPONSE.mat'],'US_I_m');
%% alltheta
nframes=length(movlabel);
controlpoint=8;
basevec = [selx(controlpoint) sely(controlpoint)] - [selx(1) sely(1)];
movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(1); sely(1)], 1, size(alltailpos,3));
% distance of control points
dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));
% theta, anti-clockwise as positive
theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
theta(theta > 180) = theta(theta > 180) - 360;
theta = -theta;
theta00 = 0;
alltheta = theta00 * ones(nframes, 1);
alltheta(movind) = theta;
%alltheta(52120:53770)=0;alltheta(56170:57600)=0;alltheta(59830:60350)=0;alltheta(98610:98780)=0;alltheta(98900:99020)=0;%20200411fish3
%alltheta(40803:43248)=0;alltheta(65590:65819)=0;alltheta(65880:66284)=0;alltheta(40667:40780)=0;alltheta(2723:2828)=0;%20200410fish2
%alltheta(7083:7256)=0;%20200411fish6
alltheta(find(abs(alltheta)<3))=0;
figure; plot(alltheta);
alltheta_n=normalize(alltheta,'zscore');
d_alltheta=diff([alltheta(1);alltheta]);
d_alltheta_n=diff([alltheta_n(1);alltheta_n]);
%% setpara
frame=[];
frame.per_cycle=20*60;
frame.us_start=round(13.2*60+1);
frame.us_dur=[0.001 0.001 0.01 0.01 0.1]*60;
trial.hab=15;
trial.acq_block_trial=6;
trial.acq_block_number=5;
trial.acq_block_interval=15;
trial.test=15;
trial.total=trial.hab+trial.test+trial.acq_block_trial*trial.acq_block_number+trial.acq_block_interval*(trial.acq_block_number-1);
stimUS=zeros(1,frame.per_cycle*trial.total);ind_US=[];kk=1;
for ss=1:trial.acq_block_number
    for tt=1:trial.acq_block_trial
        ind=(trial.hab+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frame.per_cycle+frame.us_start;
        stimUS(ind:ind+max(round(frame.us_dur(ss)),1)-1)=1;%23
        ind_US(kk,1)=ind;kk=kk+1;
    end
end
find( stimUS==1);%figure,plot( stimUS);
% figure('position',[50,50,1750,950]),
% for ii=1:length(ind_US)
%     subaxis(trial.acq_block_trial,trial.acq_block_number,ii,'SpacingVertical',0.01,'SpacingHorizontal',0.01);
%     imshow(US_I_m(:,:,ii),[min(min(US_I_m(:,:,ii))) max(max(US_I_m(:,:,ii)))]);
% end
%% plot
a=20*[-1 1];
T=[1:trial.total*frame.per_cycle];
figure('position',[50,250,1700,450]),plot(T,stimUS*180,'r','linewidth',1.5);hold on;
plot(T,alltheta_n);xlim([min(T) max(T)]);ylim(a);
T=[1:trial.total*frame.per_cycle]/60/60;
figure('position',[50,250,1700,450]),plot(T,stimUS*180,'r','linewidth',1.5);hold on;
plot(T,alltheta,'linewidth',1);xlim([min(T) max(T)]);ylim(a);xlabel('Time(min)','fontsize',16);ylabel('Tail Angle');
set(gca,'fontsize',20,'linewidth',1.5);

% avg angle
%spon
a=50*[-1 1];
Theta_avg_spon=[];Theta_sum_spon=[];M_num_spon=[];s={};
Theta_avg_spon_n=[];Theta_sum_spon_n=[];
figure('position',[50,50,1500,950]);
for ii=1:trial.acq_block_number+1
    switch ii
        case 1
            t=1:trial.hab*frame.per_cycle;
            s{ii}='Base';
        case {2,3,4,5,6}
            t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
            s{ii}=['Session I.' num2str(ii-1)];
        case 6
            t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
            s{ii}='Test';
    end
    %subplot(trial.acq_block_number+1,1,ii);
    subaxis(trial.acq_block_number+1,1,ii,'SpacingVertical',0.05,'SpacingHorizontal',0.01);
    plot(T(t),stimUS(t)*180,'r','linewidth',1.5);hold on;plot(T(t),alltheta(t));
    xlim([floor(min(T(t))) max(T(t))]);ylim(a);title(s{ii},'fontsize',16);box on;
    xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',16,'linewidth',1.5);
    [env_num,env_loc,~]=env_detect(alltheta(t),3*60,1);scatter(T(t(env_loc)),alltheta(t(env_loc)));
    
    Theta_avg_spon_n(ii)=sum(abs(d_alltheta_n(t)),'omitnan')/env_num;
    Theta_sum_spon_n(ii)=sum(abs(d_alltheta_n(t)),'omitnan')/(length(t)/60/60);
    Theta_avg_spon(ii)=sum(abs(d_alltheta(t)),'omitnan')/env_num;
    Theta_sum_spon(ii)=sum(abs(d_alltheta(t)),'omitnan')/(length(t)/60/60);
    M_num_spon(ii)=max(env_num)/(length(t)/60/60);
end
% us
Theta_avg_US=[];Theta_sum_US=[];M_num_US=[];s_us={};
Theta_avg_US_n=[];Theta_sum_US_n=[];
figure('position',[50,50,1200,800]);
for ii=1:trial.acq_block_number
    t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
        (trial.hab+ii*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
    
    subaxis(trial.acq_block_number,1,ii,'SpacingVertical',0.05,'SpacingHorizontal',0.01);
    %subplot(trial.acq_block_number,1,ii);%subplot(ceil(trial.acq_block_number/2),2,ii);
    plot(T(t),stimUS(t)*180,'r','linewidth',1.5);hold on;
    plot(T(t),alltheta(t));hold on;    plot(T(t),abs([d_alltheta(t)]));hold on;
    xlim([floor(min(T(t))) max(T(t))]);ylim(a);title(['Session ' num2str(ii)],'fontsize',16);box on;
    xticks([floor(min(T(t))) max(T(t))]);set(gca,'fontsize',16,'linewidth',1.5);
    s_us{ii}=['Session ' num2str(ii)];
    [env_num,env_loc,env_end]=env_detect(alltheta(t),3*60,1);scatter(T(t(env_loc)),alltheta(t(env_loc)),12,'g','filled');hold on;
    M_num_US(ii)=env_num;
    Theta_sum_US(ii)=sum(abs(d_alltheta(t)),'omitnan')/(length(t)/60/60);
    Theta_sum_US_n(ii)=sum(abs(d_alltheta_n(t)),'omitnan')/(length(t)/60/60);
    %     
    tt=[];
    if ~isempty(env_loc)
        scatter(T(t(env_end)),alltheta(t(env_end)),12,'k','filled');hold on;
        aa=min(abs(repmat(t(env_loc),length(ind_US),1)-repmat(ind_US,1,length(env_loc))),[],1);
        ind=find(aa<=200);env_loc(ind);
        env_loc=env_loc(ind);env_num=length(ind);env_end=env_end(ind);
        for zz=1:length(ind)
            tt=[tt,env_loc(zz):env_end(zz)];
        end
        tt=t(tt);
        scatter(T(t(env_loc)),alltheta(t(env_loc)),12,'m','filled');
    end
    %
    Theta_avg_US(ii)=sum(abs(d_alltheta(tt)),'omitnan')/env_num;
    Theta_avg_US_n(ii)=sum(abs(d_alltheta_n(tt)),'omitnan')/env_num;
    %figure,plot(alltheta(tt));hold on;plot(abs(d_alltheta(tt)));
end
%statistic
figure('position',[50 50 275 775]),
subplot(3,1,1),
bar(categorical(s),M_num_spon);colormap('parula')
ylabel('Mov. Rate/Min');set(gca,'fontsize',14);
subplot(3,1,2),
bar(categorical(s),Theta_sum_spon);colormap('parula')
ylabel('Sum Theta/Min');set(gca,'fontsize',14);
subplot(3,1,3),
bar(categorical(s),Theta_avg_spon);colormap('parula')
ylabel('Avg. Theta');set(gca,'fontsize',14);

figure('position',[50 50 275 775]),
subplot(3,1,1),
bar(categorical(s_us),M_num_US);
ylabel('Mov. Rate/Min');set(gca,'fontsize',14);
subplot(3,1,2),
bar(categorical(s_us),Theta_sum_US);
ylabel('Sum Theta/Min');set(gca,'fontsize',14);
subplot(3,1,3),
bar(categorical(s_us),Theta_avg_US);
ylabel('Avg. Theta');set(gca,'fontsize',14);

%save 
if ~exist('H:\1.Test US\summary.mat')
    Path=string;
    Data.Theta_avg_US=[];Data.Theta_sum_US=[];Data.M_num_US=[];
    Data.Theta_avg_US_n=[];Data.Theta_sum_US_n=[];
    Data.Theta_avg_spon=[];Data.Theta_sum_spon=[];Data.M_num_spon=[];
    Data.Theta_avg_spon_n=[];Data.Theta_sum_spon_n=[];
    save('H:\1.Test US\summary.mat','Path','Data')
else
    load('H:\1.Test US\summary.mat');
    ind=find(strcmp(Path,path));
    if ~isempty(ind)
        Path(ind,:)=path;
    else
        Path(end+1,:)=path;
    end
    ind=find(strcmp(Path,path));
    Data.Theta_avg_US(ind,:)=Theta_avg_US;Data.Theta_sum_US(ind,:)=Theta_sum_US;Data.M_num_US(ind,:)=M_num_US;
    Data.Theta_avg_US_n(ind,:)=Theta_avg_US_n;Data.Theta_sum_US_n(ind,:)=Theta_sum_US_n;
    Data.Theta_avg_spon(ind,:)=Theta_avg_spon;Data.Theta_sum_spon(ind,:)=Theta_sum_spon;Data.M_num_spon(ind,:)=M_num_spon; 
    Data.Theta_avg_spon_n(ind,:)=Theta_avg_spon_n;Data.Theta_sum_spon_n(ind,:)=Theta_sum_spon_n;
    save('H:\1.Test US\summary.mat','Path','Data');
end
%% movement rate
%sliding
% win=1*60;
% theta_s=nan(size(center_fish_adj,1),1);
% for ii=1:size(center_fish_adj,1)-win
%     theta_s(ii)=mean(alltheta(ii:ii+win-1,1),'omitnan');
% end
% figure('position',[50,250,1700,450]),
% plot(T,stimUS.*max(theta_s),'r--','linewidth',1.5);hold on;
% plot(T,theta_s,'b','linewidth',1.5);ylim([min(theta_s) max(theta_s)]);ylabel('Speed (pixels/s)');hold on;
% xlim([min(T) max(T)]);xlabel('Time(min)');
% title(['Slide window ' num2str(win/60) 's']);set(gca,'fontsize',16);

%env_num=env_detect(alltheta(t),0.1*60,1);

