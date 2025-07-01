v%% load
clc;clear all;close all;
[name path]=uigetfile('I:\1.Test US\2.Tail free¡ª¡ªData from 117\*.csv','Results');
A=xlsread([path name]);
[n p]=uigetfile([path '\*.mat'],'Results_of_alltheta.mat');load([p n]);

n=5;
heartbeat=A(:,3:(size(A,2)-1)/n:size(A,2));heartbeat=heartbeat(:,1:n);
%% setpata
[fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=10;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frameb.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num)
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
stimCS=zeros(1,frameb.per_cycle*trial.total);ind_CS=[];ind_preCS=[];kk=1;
for tt=trial_ind_CS
    stimCS((tt-1)*frameb.per_cycle+frameb.cs_start:(tt-1)*frameb.per_cycle+frameb.cs_end)=1;%23
    ind_preCS(kk,1)=(tt-1)*frameb.per_cycle+1;
    ind_preCS(kk,2)=(tt-1)*frameb.per_cycle+frameb.cs_start-1;
    ind_CS(kk,1)=(tt-1)*frameb.per_cycle+frameb.cs_start;
    ind_CS(kk,2)=(tt-1)*frameb.per_cycle+frameb.cs_end;kk=kk+1;
end
%figure,plot(stimCS,'r');hold on;plot(stimUS,'k');
%% process
heartbeat_s=smoothdata(heartbeat,'gaussian',5);
[dfdf, baselines]=compute_dfdf(heartbeat',101);
figure,plot(heartbeat);hold on;plot(baselines')

dfdf_s=smoothdata(dfdf,2,'movmean',10);
figure,plot(dfdf');hold on;plot(dfdf_s');
noise=dfdf_s(1,:)-dfdf_s(2,:);figure,plot(noise)

heartbeat_aft_pro=dfdf_s;[dfdf_s-repmat(noise,size(dfdf_s,1),1)]';
%figure,plot(heartbeat_aft_pro);

%% findpeaks MinPeakHeight:15%-85%
ECG_Data =heartbeat_aft_pro(1,:)-heartbeat_aft_pro(2,:);T=(1:trial.total*frameb.per_cycle);
[pk, locs] = findpeaks(ECG_Data,'MinPeakProminence',0.01,'MinPeakDistance',60/3); %, 'MinPeakHeight',0.05);
figure;plot(T,ECG_Data,locs,pk,'o');

%% statistic
peakInterval = diff(locs);
figure;subplot(1,2,1);boxplot(peakInterval)
%locs(find(abs(peakInterval>mean(peakInterval)+3*std(peakInterval)))+1)=[];
peakInterval(find(abs(peakInterval>mean(peakInterval)+3*std(peakInterval))))=[];
subplot(1,2,2);boxplot(peakInterval)

xbins = min(peakInterval):max(peakInterval);
figure;hist(peakInterval, xbins)

% heartrate stats
meanIBI = mean(peakInterval); sdIBI=std(peakInterval);
% 10hz, adjust divisor to your imaging rate.
heartRate = (1./(meanIBI/fs.behavior)).*60;
averageInterval = strcat('Average interval between heartbeats is: ',num2str(meanIBI/fs.behavior*1000),'ms');
averageHeartRate = strcat('Average heartRate is: ',num2str(heartRate),'beats per minute (bpm)');
disp(averageInterval);
disp(averageHeartRate);

meanIBI_cs=[];meanIBI_precs=[];kk=1;
for zz=1:trial.acq_block_num+2
    t=T_US{zz};
    for ii=1:length(t)/frameb.per_cycle
        [~,loc_CSz,~]=intersect(locs,[ind_CS(kk,1):ind_CS(kk,2)]);
        meanIBI_cs(zz,ii,:)=mean(diff(locs(loc_CSz)));
        [~,loc_preCSz,~]=intersect(locs,[ind_preCS(kk,1):ind_preCS(kk,2)]);
        meanIBI_precs(zz,ii,:)=mean(diff(locs(loc_preCSz)));
        kk=kk+1;
    end
end
meanIBI_spon=[];
for zz=1:trial.acq_block_num+1
    t=T_spon{zz};
    [~,loc_CSz,~]=intersect(locs,t);
    meanIBI_spon(zz,:)=mean(diff(locs(loc_CSz)));
end
%% plot
y=1:frameb.per_cycle;
colorCSbar=[0.2 0.2 0.2];
locs_low=locs(find(diff(locs)>meanIBI+2*sdIBI));
for zz=1:trial.acq_block_num+2
    t=T_US{zz};[~,locz,~]=intersect(t,locs);[~,loc_lowz,~]=intersect(t,locs_low);
    [~,tailmovlocz,~]=intersect(t,T_non_spon(env.env_loc(:,1)));
    x=ECG_Data(t);
    xr=reshape(x,frameb.per_cycle,[]);
    figure('position',[415,226,1250,250*length(x)/frameb.per_cycle/3]),title(l_US{zz});
    for ii=1:length(x)/frameb.per_cycle
        subplot(length(x)/frameb.per_cycle/3,3,ii),
        yy=(y+(ii-1)*frameb.per_cycle)/fs.behavior;
        patch([min(yy)+(frameb.cs_start-1)/fs.behavior,min(yy)+(frameb.cs_end-1)/fs.behavior,min(yy)+(frameb.cs_end-1)/fs.behavior,min(yy)+(frameb.cs_start-1)/fs.behavior]',...
            repmat([min(x) min(x) max(x) max(x)],1,1)',...
            colorCSbar,'edgecolor',colorCSbar,'facecolor',colorCSbar,'edgealpha',0.2,'facealpha',0.25);hold on;
        plot(yy,xr(:,ii),'b','linewidth',1.5);hold on;
        plot(locz/fs.behavior,x(locz),'k.','markersize',10,'markerfacecolor',[0 0 0]);hold on;
        plot(loc_lowz/fs.behavior,x(loc_lowz),'r.','markersize',16,'markerfacecolor',[0 0 0]);hold on;
        aa=[min(x) max(x)];
        a=(aa(2)-aa(1))/10+aa(1);
         if ~isempty(tailmovlocz)
        plot(tailmovlocz/fs.behavior,a,'k^','markersize',4,'markerfacecolor',[0 0 0]);
         end
        xlim([min(yy)+5 max(yy)-5]);ylim(aa);
        set(gca,'fontsize',14);
    end
end
figure('position',[3,3,1914,120*trial.acq_block_num+1]);
for zz=1:trial.acq_block_num+1
    subplot(trial.acq_block_num+1,1,zz)
    t=T_spon{zz};[~,locz,~]=intersect(t,locs);[~,loc_lowz,~]=intersect(t,locs_low);
    [~,env_loc,env_end,~]=env_detect(alltheta(t),3*60,1);
    [~,tailmovlocz,~]=intersect(t,t(env_loc));
    x=ECG_Data(t);
    yy=t/fs.behavior;
    plot(yy,x,'b','linewidth',1.5);hold on;
    plot(t(locz)/fs.behavior,x(locz),'k.','markersize',10,'markerfacecolor',[0 0 0]);hold on;
    plot(t(loc_lowz)/fs.behavior,x(loc_lowz),'r.','markersize',16,'markerfacecolor',[0 0 0]);hold on;
    a=(max(x)-min(x))/10+min(x);
    if ~isempty(tailmovlocz)
    plot(t(tailmovlocz)/fs.behavior,a,'k^','markersize',4,'markerfacecolor',[0 0 0]);
    end
    xlim([min(yy)+5 max(yy)-5]);ylim([min(x) max(x)]);title(l_spon{zz});
    set(gca,'fontsize',14);
end
save([path 'hearbeat.mat'],'heartbeat','heartbeat_aft_pro','ECG_Data','locs','peakInterval','meanIBI','sdIBI','heartRate','meanIBI_cs','meanIBI_precs','meanIBI_spon');

close all

