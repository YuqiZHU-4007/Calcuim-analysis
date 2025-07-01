clear all; close all; clc;
%%
%parameters
load('heartbeat.mat');
fps = 60;
totframe_percycle = 1200;
CSstart_frame = 601;
CSend_frame = 889;
US_frame = 841;
trials = fix(size(heartbeat,1)/totframe_percycle);
aheartbeat = heartbeat(1:totframe_percycle*trials,7);
trial_hab = 15;
trial_acq = 30;
trial_test = 15;
MinPeakHeight = 2;
MinPeakDistance = 0.3;

for i = 1 : trials
    base1 = mean(aheartbeat(1+(i-1)*totframe_percycle:(CSstart_frame-1)+(i-1)*totframe_percycle,1));
    baseline1 = base1 * ones((CSstart_frame-1),1);
    mheartbeat(1+(i-1)*totframe_percycle:(CSstart_frame-1)+(i-1)*totframe_percycle,1) = aheartbeat(1+(i-1)*totframe_percycle:(CSstart_frame-1)+(i-1)*totframe_percycle,1)-baseline1;
end
for i = 1 : trials
    base2 = mean(aheartbeat(CSstart_frame+(i-1)*totframe_percycle:(CSend_frame-1)+(i-1)*totframe_percycle,1));
    baseline2 = base2 * ones((CSend_frame-CSstart_frame),1);
    mheartbeat(CSstart_frame+(i-1)*totframe_percycle:(CSend_frame-1)+(i-1)*totframe_percycle,1) = aheartbeat(CSstart_frame+(i-1)*totframe_percycle:(CSend_frame-1)+(i-1)*totframe_percycle,1)-baseline2;
end
for i = 1 : trials
    base3 = mean(aheartbeat(CSend_frame+(i-1)*totframe_percycle:totframe_percycle+(i-1)*totframe_percycle,1));
    baseline3 = base3 * ones((totframe_percycle-CSend_frame+1),1);
    mheartbeat(CSend_frame+(i-1)*totframe_percycle:totframe_percycle+(i-1)*totframe_percycle,1) = aheartbeat(CSend_frame+(i-1)*totframe_percycle:totframe_percycle+(i-1)*totframe_percycle,1)-baseline3;
end


% %bin
% bin = 4;
% aftbin_heartbeat = [];
% for i = [1 : bin : (size(mheartbeat,1)-(bin-1))]
%     j = mean(mheartbeat(i:(i+bin-1),1));
%     aftbin_heartbeat(((i+(bin-1))/bin),1) = j;
% end

% %smooth
% mheartbeat_sm200=smooth(aftbin_heartbeat,200);
% figure;plot(aftbin_heartbeat-mheartbeat_sm200);


%habutation,acquisation,test ????
heartbeat_hab = mheartbeat(1:(totframe_percycle)*trial_hab,1);
heartbeat_acq = mheartbeat((totframe_percycle)*trial_hab+1:(totframe_percycle)*(trial_hab+trial_acq),1);
heartbeat_test = mheartbeat((totframe_percycle)*(trial_hab+trial_acq)+1:(totframe_percycle)*trials,1);

resh_heartbeat = reshape(mheartbeat, (totframe_percycle), trials)';
resh_heartbeathab = reshape(heartbeat_hab, (totframe_percycle), trial_hab)';
resh_heartbeatacq = reshape(heartbeat_acq, (totframe_percycle), trial_acq)';
resh_heartbeattest = reshape(heartbeat_test, (totframe_percycle), (trials-trial_hab-trial_acq))';


aftcon_heartbeathab =resh_heartbeathab(:,(CSstart_frame:totframe_percycle));
aftcon_heartbeatacq =resh_heartbeatacq(:,(CSstart_frame:totframe_percycle));
aftcon_heartbeattest =resh_heartbeattest(:,(CSstart_frame:totframe_percycle));

%%
%对trial之间的心律进行比较
for i=1:size(resh_heartbeathab,1)
    [resh_pksHB_hab,resh_locsHB_hab] = findpeaks(resh_heartbeathab(i,:),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum_hab(i,1)=length(resh_pksHB_hab);
end
for i=1:size(resh_heartbeatacq,1)
    [resh_pksHB_acq_bfus,resh_locsHB_acq_bfus] = findpeaks(resh_heartbeatacq(i,1:floor(US_frame)),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum_acqbfus(i,1)=length(resh_pksHB_acq_bfus);
    [resh_pksHB_acq_afus,resh_locsHB_acq_afus] = findpeaks(resh_heartbeatacq(i,floor((US_frame))+1:end),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum_acqafus(i,1)=length(resh_pksHB_acq_afus);
end
for i=1:size(resh_heartbeattest,1)
    [resh_pksHB_test,resh_locsHB_test] = findpeaks(resh_heartbeattest(i,:),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum_test(i,1)=length(resh_pksHB_test);
end

for i=1:size(resh_heartbeat)
    [resh_pksHB,resh_locsHB] = findpeaks(resh_heartbeat(i,1:floor(CSstart_frame)),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum(i,1)=length(resh_pksHB);
    [resh_pksHB,resh_locsHB] = findpeaks(resh_heartbeat(i,floor(CSstart_frame)+1:floor(US_frame)),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum(i,2)=length(resh_pksHB);
    [resh_pksHB,resh_locsHB] = findpeaks(resh_heartbeat(i,floor(US_frame)+1:floor(CSend_frame)),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum(i,3)=length(resh_pksHB);
    [resh_pksHB,resh_locsHB] = findpeaks(resh_heartbeat(i,floor(CSend_frame)+1:end),(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    resh_pksnum(i,4)=length(resh_pksHB);
end
mean_resh_pksnum(1,:)=mean(resh_pksnum(1:trial_hab,:));
mean_resh_pksnum(2,:)=mean(resh_pksnum((trial_hab+1):(trial_hab+10),:));
mean_resh_pksnum(3,:)=mean(resh_pksnum((trial_hab+trial_acq-9):(trial_hab+trial_acq),:));
mean_resh_pksnum(4,:)=mean(resh_pksnum((trial_hab+trial_acq+1):trials,:));

% figure,
% for i=1:4
% subplot(2,2,i),scatter(1:45,resh_pksnum(:,i))/(CSstart_frame/fps);title('before CS');xlim([1 45]);ylim([0 max(resh_pksnum(:,i))+1])
% hold on;line([15 15],[0 max(resh_pksnum(:,i))+1],'Color','red');
% hold on;line([35 35],[0 max(resh_pksnum(:,i))+1],'Color','red');
% hold on;line([25 25],[0 max(resh_pksnum(:,i))+1],'Color','red','LineStyle','--');
% a=(max(resh_pksnum(:,1))+1)/10;b=a*2;
% text(5,a,['Var:' num2str(var(resh_pksnum(1:15,i)))]);text(5,b,['Mean:' num2str(mean(resh_pksnum(1:15,i)))]);
% text(20,a,['Var:' num2str(var(resh_pksnum(16:25,i)))]);text(20,b,['Mean:' num2str(mean(resh_pksnum(16:25,i)))]);
% text(30,a,['Var:' num2str(var(resh_pksnum(26:35,i)))]);text(30,b,['Mean:' num2str(mean(resh_pksnum(26:35,i)))]);
% text(40,a,['Var:' num2str(var(resh_pksnum(36:45,i)))]);text(40,b,['Mean:' num2str(mean(resh_pksnum(36:45,i)))]);    
% end

figure;input=resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps);xlabel('trials','FontSize',20);ylabel('heartbeats/min','FontSize',20);
scatter(1:trials,resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps));title('During CS (before US)','FontSize',20);
ylim([0 max(resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps))+1]);xlim([1 trials])
hold on;line([(trial_hab+0.5) (trial_hab+0.5)],[0 max(resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps))+1],'Color','red');
hold on;line([(trial_hab+trial_acq+0.5) (trial_hab+trial_acq+0.5)],[0 max(resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps))+1],'Color','red');
% hold on;line([(trial_hab+10+0.5) (trial_hab+10+0.5)],[0 max(resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps))+1],'Color','red','LineStyle','--');
a=(max(resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps))+1)/10;b=a*2;
text(5,a,['Var:' num2str(var(input(1:trial_hab,1)))]);text(5,b,['Mean:' num2str(mean(input(1:trial_hab,1)))]);
text(20,a,['Var:' num2str(var(input((trial_hab+1):(trial_hab+10),1)))]);text(20,b,['Mean:' num2str(mean(input((trial_hab+1):(trial_hab+10),1)))]);
text(30,a,['Var:' num2str(var(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);text(30,b,['Mean:' num2str(mean(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);
text(40,a,['Var:' num2str(var(input((trial_hab+trial_acq+1):trials,1)))]);text(40,b,['Mean:' num2str(mean(input((trial_hab+trial_acq+1):trials,1)))]);

figure;
subplot(3,1,1),input=resh_pksnum(:,1)/(CSstart_frame/fps);xlabel('trials','FontSize',20);ylabel('heartbeats/min','FontSize',20);
scatter(1:trials,resh_pksnum(:,1)/(CSstart_frame/fps));title('Before CS','FontSize',20);xlim([1 trials]);ylim([0 max(resh_pksnum(:,1)/(CSstart_frame/fps))+1])
hold on;line([(trial_hab+0.5) (trial_hab+0.5)],[0 max(resh_pksnum(:,1)/(CSstart_frame/fps))+1],'Color','red');
hold on;line([(trial_hab+trial_acq+0.5) (trial_hab+trial_acq+0.5)],[0 max(resh_pksnum(:,1)/(CSstart_frame/fps))+1],'Color','red');
% hold on;line([(trial_hab+10+0.5) (trial_hab+10+0.5)],[0 max(resh_pksnum(:,1)/(CSstart_frame/fps))+1],'Color','red','LineStyle','--');
a=(max(resh_pksnum(:,1)/(CSstart_frame/fps))+1)/(trial_acq/2);b=a*2;
text(5,a,['Var:' num2str(var(input(1:trial_hab,1)))]);text(5,b,['Mean:' num2str(mean(input(1:trial_hab,1)))]);
text(20,a,['Var:' num2str(var(input((trial_hab+1):(trial_hab+10),1)))]);text(20,b,['Mean:' num2str(mean(input((trial_hab+1):(trial_hab+10),1)))]);
text(30,a,['Var:' num2str(var(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);text(30,b,['Mean:' num2str(mean(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);
text(40,a,['Var:' num2str(var(input((trial_hab+trial_acq+1):trials,1)))]);text(40,b,['Mean:' num2str(mean(input((trial_hab+trial_acq+1):trials,1)))]);

subplot(3,1,2),input=resh_pksnum(:,3)/((CSend_frame-US_frame)/fps);xlabel('trials','FontSize',20);ylabel('heartbeats/min','FontSize',20);
scatter(1:trials,resh_pksnum(:,3)/((CSend_frame-US_frame)/fps));title('During CS (after US)','FontSize',20);
ylim([0 max(input)+1]);xlim([1 trials])
hold on;line([(trial_hab+0.5) (trial_hab+0.5)],[0 max(input)+1],'Color','red');
hold on;line([(trial_hab+trial_acq+0.5) (trial_hab+trial_acq+0.5)],[0 max(input)+1],'Color','red');
% hold on;line([(trial_hab+10+0.5) (trial_hab+10+0.5)],[0 max(input)+1],'Color','red','LineStyle','--');
a=(max(input)+1)/10;b=a*2;
text(5,a,['Var:' num2str(var(input(1:trial_hab,1)))]);text(5,b,['Mean:' num2str(mean(input(1:trial_hab,1)))]);
text(20,a,['Var:' num2str(var(input((trial_hab+1):(trial_hab+10),1)))]);text(20,b,['Mean:' num2str(mean(input((trial_hab+1):(trial_hab+10),1)))]);
text(30,a,['Var:' num2str(var(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);text(30,b,['Mean:' num2str(mean(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);
text(40,a,['Var:' num2str(var(input((trial_hab+trial_acq+1):trials,1)))]);text(40,b,['Mean:' num2str(mean(input((trial_hab+trial_acq+1):trials,1)))]);

subplot(3,1,3),input=resh_pksnum(:,4)/((totframe_percycle-CSend_frame)/fps);xlabel('trials','FontSize',20);ylabel('heartbeats/min','FontSize',20);
scatter(1:trials,resh_pksnum(:,4)/((totframe_percycle-CSend_frame)/fps));title('After CS','FontSize',20);
ylim([0 max(input)+1]);xlim([1 trials])
hold on;line([(trial_hab+0.5) (trial_hab+0.5)],[0 max(input)+1],'Color','red');
hold on;line([(trial_hab+trial_acq+0.5) (trial_hab+trial_acq+0.5)],[0 max(input)+1],'Color','red');
% hold on;line([(trial_hab+10+0.5) (trial_hab+10+0.5)],[0 max(input)+1],'Color','red','LineStyle','--');
a=(max(input)+1)/10;b=a*2;
text(5,a,['Var:' num2str(var(input(1:trial_hab,1)))]);text(5,b,['Mean:' num2str(mean(input(1:trial_hab,1)))]);
text(20,a,['Var:' num2str(var(input((trial_hab+1):(trial_hab+10),1)))]);text(20,b,['Mean:' num2str(mean(input((trial_hab+1):(trial_hab+10),1)))]);
text(30,a,['Var:' num2str(var(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);text(30,b,['Mean:' num2str(mean(input((trial_hab+trial_acq-9):(trial_hab+trial_acq),1)))]);
text(40,a,['Var:' num2str(var(input((trial_hab+trial_acq+1):trials,1)))]);text(40,b,['Mean:' num2str(mean(input((trial_hab+trial_acq+1):trials,1)))]);

resh_pksnum_befCS_mean = resh_pksnum(:,1)/(CSstart_frame/fps);
resh_pksnum_durCS1_mean = resh_pksnum(:,2)/((US_frame-CSstart_frame)/fps);
resh_pksnum_durCS2_mean = resh_pksnum(:,3)/((CSend_frame-US_frame)/fps);
resh_pksnum_durCS_mean = (resh_pksnum(:,2) + resh_pksnum(:,3))/((CSend_frame-CSstart_frame)/fps);
resh_pksnum_aftCS_mean = resh_pksnum(:,4)/((totframe_percycle-CSend_frame)/fps);

relat_HB_frequency_durCS1 = resh_pksnum_durCS1_mean./resh_pksnum_befCS_mean;
relat_HB_frequency_durCS2 = resh_pksnum_durCS2_mean./resh_pksnum_befCS_mean;
relat_HB_frequency_durCS = resh_pksnum_durCS_mean./resh_pksnum_befCS_mean;
relat_HB_frequency_aftCS = resh_pksnum_aftCS_mean./resh_pksnum_befCS_mean;

save('heartbeat_freq_crosstrial.mat','resh_pksnum_befCS_mean','resh_pksnum_durCS1_mean',...
    'resh_pksnum_durCS2_mean','resh_pksnum_durCS_mean','resh_pksnum_aftCS_mean','relat_HB_frequency_durCS1',...
    'relat_HB_frequency_durCS2','relat_HB_frequency_durCS','relat_HB_frequency_aftCS');

%%
%????????
[pks_HB,locs_HB] = findpeaks(mheartbeat,(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
%find the hab,acq and test session
[pksHB_hab,locsHB_hab] = findpeaks(heartbeat_hab,(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
[pksHB_acq,locsHB_acq] = findpeaks(heartbeat_acq,(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
[pksHB_test,locsHB_test] = findpeaks(heartbeat_test,(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
period_HB = diff(locs_HB);
figure;plot(period_HB)
freq_HB = 1./period_HB;
Locs_HB= ([locs_HB;NaN(1)]+[NaN(1);locs_HB])/2; %correct time point
Locs_HB = Locs_HB(2:(size(Locs_HB,1)-1),:);
LocsHB_hab= ([locsHB_hab;NaN(1)]+[NaN(1);locsHB_hab])/2; %correct time point
LocsHB_hab = LocsHB_hab(2:(size(LocsHB_hab,1)-1),:);
LocsHB_acq= ([locsHB_acq;NaN(1)]+[NaN(1);locsHB_acq])/2; %correct time point
LocsHB_acq = LocsHB_acq(2:(size(LocsHB_acq,1)-1),:);
LocsHB_test= ([locsHB_test;NaN(1)]+[NaN(1);locsHB_test])/2; %correct time point
LocsHB_test = LocsHB_test(2:(size(LocsHB_test,1)-1),:);
figure;plot(Locs_HB,freq_HB)   %plot heartbeat frequency
    xlabel('Time(s)','FontSize',20)
    ylabel('Frequency(Hz)','FontSize',20)
hold on
plot([1 locsHB_hab(size(locsHB_hab,1))],[max(freq_HB)+1 max(freq_HB)+1],...
    'linestyle','-',...
    'color','c',...
    'LineWidth',20);
plot([locsHB_hab(size(locsHB_hab,1))+locsHB_acq(1,1) locsHB_hab(size(locsHB_hab,1))+locsHB_acq(size(locsHB_acq,1))],...
    [max(freq_HB)+1 max(freq_HB)+1],...
    'linestyle','-',...
    'color','m',...
    'LineWidth',20);
plot([locsHB_hab(size(locsHB_hab,1))+locsHB_acq(size(locsHB_acq,1))+locsHB_test(1,1) locsHB_hab(size(locsHB_hab,1))+locsHB_acq(size(locsHB_acq,1))+locsHB_test(size(locsHB_test,1))],...
    [max(freq_HB)+1 max(freq_HB)+1],...
    'linestyle','-',...
    'color','g',...
    'LineWidth',20);
hold off

[pks_DiffHB,locs_DiffHB] = findpeaks(diff(locs_HB),'MinPeakHeight',0.6);
figure;plot((1:size(mheartbeat,1))/(fps),mheartbeat);
hold on %?????peaks
plot(locs_HB,pks_HB,'o','color','g');
    xlabel('Time(s)','FontSize',20)
    ylabel('Amplitude(a.u.)','FontSize',20)
hold off

figure;plot((1:size(mheartbeat,1))/(fps),mheartbeat);
hold on %??????????peak???
plot(locs_HB(locs_DiffHB),pks_HB(locs_DiffHB),'o',...
    'MarkerSize',6,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor','m');
%??CS???
for i = 1:trials
plot([(1+totframe_percycle*(i-1))/fps ((CSstart_frame-1)+totframe_percycle*(i-1))/fps],[max(mheartbeat)+1 max(mheartbeat)+1],...
    'linestyle','-',...
    'color','r',...
    'LineWidth',20);
plot([(CSstart_frame+totframe_percycle*(i-1))/fps ((CSend_frame-1)+totframe_percycle*(i-1))/fps],[max(mheartbeat)+1 max(mheartbeat)+1],...
    'linestyle','-',...
    'color','k',...
    'LineWidth',20);
plot([(3553+totframe_percycle*(i-1))/fps totframe_percycle*i/fps],[max(mheartbeat)+1 max(mheartbeat)+1],...
    'linestyle','-',...
    'color','r',...
    'LineWidth',20);
end
%??????
for i = (trial_hab+1):(trial_hab+trial_acq)
plot([(US_frame+totframe_percycle*(i-1))/fps (US_frame+totframe_percycle*(i-1))/fps],[-(max(mheartbeat)) max(mheartbeat)],...
    'linestyle','--',...
    'color','k',...
    'LineWidth',2)
end
for i = 1 : trials
plot([totframe_percycle*i/fps totframe_percycle*i/fps],[max(mheartbeat) max(mheartbeat)+2],...
    'linestyle','-',...
    'color','w',...
    'LineWidth',2)
end
    xlabel('Time(s)','FontSize',20)
    ylabel('Amplitude(a.u.)','FontSize',20)
hold off


%%
%frequency in every trial
FreqHB_PerTrial =  zeros(50,trials);
LocsHB_PerTrial =  zeros(50,trials);
for i = 1:trials
    [pksHB_PerTrial,locsHB_PerTrial] = findpeaks(resh_heartbeat(i,:)',fps,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    LocsHB_PerTrial(1:size(locsHB_PerTrial,1),i) = locsHB_PerTrial;
    PeriodHB_PerTrial = diff(locsHB_PerTrial);
    freqHB_PerTrial = 1./PeriodHB_PerTrial;
    FreqHB_PerTrial(1:size(freqHB_PerTrial,1),i) = freqHB_PerTrial;
end 
%save locus of frequency
LocsHB_PerTrial(all(LocsHB_PerTrial==0,2),:) = [];  %delete 0=[]
LocsHB_PerTrial(LocsHB_PerTrial==0) = NaN; 
LocsHB_PerTrial=([LocsHB_PerTrial;NaN(1,trials)]+[NaN(1,trials);LocsHB_PerTrial])/2; %correct time point
LocsHB_PerTrial = LocsHB_PerTrial(2:(size(LocsHB_PerTrial,1)-1),:);
LocsHB_HabPerTrial = LocsHB_PerTrial(:,1:trial_hab);
LocsHB_AcqPerTrial = LocsHB_PerTrial(:,(trial_hab+1):(trial_acq+trial_hab));
LocsHB_TestPerTrial = LocsHB_PerTrial(:,(trial_hab+trial_acq+1):trials);
%save frequency
FreqHB_PerTrial(all(FreqHB_PerTrial==0,2),:) = [];  %delete 0=[]
FreqHB_PerTrial(FreqHB_PerTrial==0) = NaN; 
% 
% for i = 1:size(LocsHB_PerTrial,2)
%     for j = 0.5:0.5:20 
%     a = find(LocsHB_PerTrial(:,i)>=(j-0.5) & LocsHB_PerTrial(:,i)< j);
%     b = size(a)*2;
%     HB_PeaksPerBin(2*j,i) = mean(b);
%     end
% end

for i = 1:size(LocsHB_PerTrial,2)
    for j = 0.5:0.5:20 
    a = find(LocsHB_PerTrial(:,i)>=(j-0.5) & LocsHB_PerTrial(:,i)< j);
    b = FreqHB_PerTrial(a);
    FreqHB_PerTrialStand(2*j,i) = mean(b);
    end
end

%calculate the mean heartbeat in hab, acq and test
FreqHB_HabPerTrial = FreqHB_PerTrial(:,1:trial_hab);
for i = 1:size(LocsHB_HabPerTrial,2)
    for j = 0.5:0.5:20 
    a_hab = find(LocsHB_HabPerTrial(:,i)>=(j-0.5) & LocsHB_HabPerTrial(:,i)< j);
    b_acq = FreqHB_HabPerTrial(a_hab);
    FreqHB_HabPerTrialStand(2*j,i) = mean(b_acq);
    end
end
MeanHB_Hab = mean(FreqHB_HabPerTrialStand,2,'omitnan');
figure;subplot(3,1,1)
plot((1:size(MeanHB_Hab,1))/2-10,MeanHB_Hab,'LineWidth',2)
    title('Trial1-15','FontSize',20)
%     xlabel('Time(s)','FontSize',20)
    ylabel('Frequency(Hz)','FontSize',20)
hold on
plot([CSstart_frame/fps-(CSstart_frame-1)/fps (CSend_frame-1)/fps-(CSstart_frame-1)/fps],[max(MeanHB_Hab)+0.1 max(MeanHB_Hab)+0.1],...
    'linestyle','-',...
    'color','k',...
    'LineWidth',20);
hold off

FreqHB_AcqPerTrial = FreqHB_PerTrial(:,(trial_hab+1):(trial_acq+trial_hab));
for i = 1:size(LocsHB_AcqPerTrial,2)
    for j = 0.5:0.5:20 
    a_acq = find(LocsHB_AcqPerTrial(:,i)>=(j-0.5) & LocsHB_AcqPerTrial(:,i)< j);
    b_acq = FreqHB_AcqPerTrial(a_acq);
    FreqHB_AcqPerTrialStand(2*j,i) = mean(b_acq);
    end
end
MeanHB_Acq = mean(FreqHB_AcqPerTrialStand,2,'omitnan');
subplot(3,1,2)
plot((1:size(MeanHB_Acq,1))/2-10,MeanHB_Acq,'LineWidth',2')
    title('Trial16-35','FontSize',20)
%     xlabel('Time(s)','FontSize',20)
    ylabel('Frequency(Hz)','FontSize',20)
hold on
plot([CSstart_frame/fps-(CSstart_frame-1)/fps (CSend_frame-1)/fps-(CSstart_frame-1)/fps],[max(MeanHB_Acq)+0.1 max(MeanHB_Acq)+0.1],...
    'linestyle','-',...
    'color','k',...
    'LineWidth',20);
% plot([US_frame/fps-(CSstart_frame-1)/fps US_frame/fps-(CSstart_frame-1)/fps],[min(MeanHB_Acq) max(MeanHB_Acq)],...
%     'linestyle','--',...
%     'color','r',...
%     'LineWidth',2)
hold off

FreqHB_TestPerTrial = FreqHB_PerTrial(:,(trial_hab+trial_acq+1):trials);
for i = 1:size(LocsHB_TestPerTrial,2)
    for j = 0.5:0.5:20 
    a_test = find(LocsHB_TestPerTrial(:,i)>=(j-0.5) & LocsHB_TestPerTrial(:,i)< j);
    b_test = FreqHB_TestPerTrial(a_test);
    FreqHB_TestPerTrialStand(2*j,i) = mean(b_test);
    end
end
MeanHB_Test= mean(FreqHB_TestPerTrialStand,2,'omitnan');
subplot(3,1,3)
plot((1:size(MeanHB_Test,1))/2-10,MeanHB_Test,'LineWidth',2)
    title('Trial36-45','FontSize',20)
    xlabel('Time(s)','FontSize',20)
    ylabel('Frequency(Hz)','FontSize',20)
hold on
plot([CSstart_frame/fps-(CSstart_frame-1)/fps (CSend_frame-1)/fps-(CSstart_frame-1)/fps],[max(MeanHB_Test)+0.1 max(MeanHB_Test)+0.1],...
    'linestyle','-',...
    'color','k',...
    'LineWidth',20);
hold off

save('heartbeat_freq_in_trial.mat','FreqHB_PerTrialStand','FreqHB_HabPerTrialStand','FreqHB_AcqPerTrialStand','FreqHB_TestPerTrialStand');

% MeanHB = [MeanHB_Hab(1:39,1);MeanHB_Acq(1:39,1);MeanHB_Test(1:39,1)];
% figure;plot((1:size(MeanHB,1))/2,MeanHB,'LineWidth',2)
%     xlabel('Time(s)','FontSize',20)
%     ylabel('Frequency(Hz)','FontSize',20)
% hold on
% for i = 1 : 3
% plot([CSstart_frame/fps+(i-1)*size(MeanHB,1)/6 (CSend_frame-1)/fps+(i-1)*size(MeanHB,1)/6],[max(MeanHB)+0.1 max(MeanHB)+0.1],...
%     'linestyle','-',...
%     'color','k',...
%     'LineWidth',20);
% end
% plot([(US_frame/fps + size(MeanHB,1)/6) (US_frame/fps + size(MeanHB,1)/6)],[1.5 2.5],...
%     'linestyle','--',...
%     'color','r',...
%     'LineWidth',2)
% plot([1 size(MeanHB,1)/6],[max(MeanHB)+0.2 max(MeanHB)+0.2],...
%     'linestyle','-',...
%     'color','c',...
%     'LineWidth',20);
% plot([1+size(MeanHB,1)/6 2*size(MeanHB,1)/6],[max(MeanHB)+0.2 max(MeanHB)+0.2],...
%     'linestyle','-',...
%     'color','m',...
%     'LineWidth',20);
% plot([1+2*size(MeanHB,1)/6 3*size(MeanHB,1)/6],[max(MeanHB)+0.2 max(MeanHB)+0.2],...  
%     'linestyle','-',...
%     'color','g',...
%     'LineWidth',20);
% hold off

% HB_HabPerTrial = zeros(50,trial_hab);
% for i = 1:trial_hab
%     [pksHB_HabPerTrial,locsHB_HabPerTrial] = findpeaks(resh_heartbeathab(i,:)',(fps),'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
%     periodHB_HabPerTrial = diff(locsHB_HabPerTrial);
%     freqHB_HabPerTrial = 1./periodHB_HabPerTrial;
%     HB_HabPerTrial(1:size(freqHB_HabPerTrial,1),i) = freqHB_HabPerTrial;
% end
% HB_HabPerTrial(HB_HabPerTrial==0) = NaN; 


%??trial????????????
t = (1:size(resh_heartbeat,2))/(fps);
figure;imagesc(t-(CSstart_frame-1)/fps,1:trials,resh_heartbeat(1:end,:));
 colorbar;
    xlabel('Time(s)','FontSize',20)
    ylabel('Trials','FontSize',20)
   % colormap hot
hold on 
 dark_xbar = 2.4;
 dark_ybar = trials+1;
 bar_1 = bar(dark_xbar,dark_ybar,4.8, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
 set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);   
plot([US_frame/fps-(CSstart_frame-1)/fps US_frame/fps-(CSstart_frame-1)/fps],[(trial_hab+1) (trial_hab+trial_acq)],...   
    'linestyle','--',...
    'color','w',...
    'LineWidth',3)
hold off

% %%
% %plot all of the trials
% for i = 1 : trial_hab
%     figure;plot((1:size(resh_heartbeat,2))/fps-(CSstart_frame-1)/fps,resh_heartbeat(i,:));
%         xlabel('Time(s)','FontSize',20)
%         ylabel('Amplitude(a.u.)','FontSize',20)
%     hold on
%     x = [0 (CSend_frame-CSstart_frame)/fps (CSend_frame-CSstart_frame)/fps 0];
%     y = [-max(resh_heartbeat(:))*0.8 -max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8];
%     patch_1 = patch(x,y,'black');
%     set(patch_1,'Facealpha',0.3,'Edgealpha',0.3);
% %     dark_xbar = 2.4;
% %     dark_ybar = max(resh_heartbeat(:))*0.9;
% %     bar_1 = bar(dark_xbar,dark_ybar,4.8, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
% %     set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);
% %     dark_xbar = 2.4;
% %     dark_ybar = -max(resh_heartbeat(:))*0.9;
% %     bar_1 = bar(dark_xbar,dark_ybar,4.8, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
% %     set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);   
% end
% 
% for i = (trial_hab+1) : (trial_hab + trial_acq)
%     figure;x = [0 (CSend_frame-CSstart_frame)/fps (CSend_frame-CSstart_frame)/fps 0];
%     y = [-max(resh_heartbeat(:))*0.8 -max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8];
%     patch_1 = patch(x,y,'black');
%     set(patch_1,'Facealpha',0.3,'Edgealpha',0.3);
%     hold on
%     plot((1:size(resh_heartbeat,2))/fps-(CSstart_frame-1)/fps,resh_heartbeat(i,:));
%         xlabel('Time(s)','FontSize',20)
%         ylabel('Amplitude(a.u.)','FontSize',20)
% %     plot([CSstart_frame/fps-(CSstart_frame-1)/fps (CSend_frame-1)/fps-(CSstart_frame-1)/fps],[max(resh_heartbeat(:))+1 max(resh_heartbeat(:))+1],...
% %     'linestyle','-',...
% %     'color','k',...
% %     'LineWidth',20);
%     plot([US_frame/fps-(CSstart_frame-1)/fps US_frame/fps-(CSstart_frame-1)/fps],[-max(resh_heartbeat(:))*0.9 max(resh_heartbeat(:))*0.9],...
%     'linestyle','--',...
%     'color','r',...
%     'LineWidth',2)
%     hold off 
% end
% 
% for i = (trial_hab + trial_acq+1) : trials
%     figure;plot((1:size(resh_heartbeat,2))/fps-(CSstart_frame-1)/fps,resh_heartbeat(i,:));
%         xlabel('Time(s)','FontSize',20)
%         ylabel('Amplitude(a.u.)','FontSize',20)
%     hold on
%     x = [0 (CSend_frame-CSstart_frame)/fps (CSend_frame-CSstart_frame)/fps 0];
%     y = [-max(resh_heartbeat(:))*0.8 -max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8 max(resh_heartbeat(:))*0.8];
%     patch_1 = patch(x,y,'black');
%     set(patch_1,'Facealpha',0.3,'Edgealpha',0.3);
% end
   
   
%%
%所有的trial都plot在一张图
% y_total_incre = zeros;
% for i = 1 : trials;
%     y_per_incre = 0.8*max(abs(resh_heartbeat(i,:)));
%     y_total_incre(i,:) =  y_per_incre;
% end
% y_final_incre = ((mean(y_total_incre) + std(y_total_incre)));
% %habitation
% figure;
% for i =  5 : trial_hab
%     hold on
%     plot((1:size(resh_heartbeat,2))/fps-3.6,resh_heartbeat(i,:) + y_final_incre*(i-4));
%     %y_per_incre = max(abs(resh_alltheta(i,:)));
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + y);
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + 60*(i-1));
% end
% hold on
% %  ylim([-120 y_final_incre*(trials+1)]);
%  title('Habitation','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Trials','FontSize',20)  
% hold on 
%  dark_xbar = 2;
%  dark_ybar = 90;
%  bar_1 = bar(dark_xbar,dark_ybar,4, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
%  set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);      
% %acqisition (trial_hab+1)-25
%  figure;
% for i =  (trial_hab+1): 25
%     hold on
%     plot((1:size(resh_heartbeat,2))/fps-3.6,resh_heartbeat(i,:) + y_final_incre*(i-trial_hab));
%     %y_per_incre = max(abs(resh_alltheta(i,:)));
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + y);
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + 60*(i-1));
% end
% hold on 
%  dark_xbar = 2;
%  dark_ybar = 80;
%  bar_1 = bar(dark_xbar,dark_ybar,4, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
%  set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);   
% hold on
% plot([1.6 1.6],[0 80],...
%     'linestyle','--',...
%     'color','w',...
%     'LineWidth',3)
% plot([2 2],[0 80],...
%        'linestyle','--',...
%        'Color','w',...
%        'LineWidth',3)
%    title('Acquisition(16-25)','FontSize',20)
%    xlabel('Time(s)','FontSize',20)
%    ylabel('Trials','FontSize',20) 
%  %acqisition 26-(trial_hab+trial_acq)
%  figure;
% for i =  26: (trial_hab+trial_acq)
%     hold on
%     plot((1:size(resh_heartbeat,2))/fps-3.6,resh_heartbeat(i,:) + y_final_incre*(i-25));
%     %y_per_incre = max(abs(resh_alltheta(i,:)));
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + y);
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + 60*(i-1));
% end
% hold on 
%  dark_xbar = 2;
%  dark_ybar = 80;
%  bar_1 = bar(dark_xbar,dark_ybar,4, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
%  set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);   
% hold on
% plot([1.6 1.6],[0 80],...
%     'linestyle','--',...
%     'color','w',...
%     'LineWidth',3)
% plot([2 2],[0 80],...
%        'linestyle','--',...
%        'Color','w',...
%        'LineWidth',3)
%    title('Acquisition(26-35)','FontSize',20)
%    xlabel('Time(s)','FontSize',20)
%    ylabel('Trials','FontSize',20)  
% % test   
%  figure;
% for i =  (trial_hab+trial_acq+1) : trials
%     hold on
%     plot((1:size(resh_heartbeat,2))/fps-3.6,resh_heartbeat(i,:) + 0.8*y_final_incre*(i-(trial_hab+trial_acq)));
%     %y_per_incre = max(abs(resh_alltheta(i,:)));
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + y);
%     %plot((1:length(resh_alltheta))/fps-30,resh_alltheta(i,:) + 60*(i-1));
% end
% hold on 
%  dark_xbar = 2;
%  dark_ybar = 80;
%  bar_1 = bar(dark_xbar,dark_ybar,4, 'FaceColor', [0 0 0], 'EdgeColor', [0 0 0]);
%  set(bar_1,'Facealpha',0.3,'Edgealpha',0.3);   
%    title('Test','FontSize',20)
%    xlabel('Time(s)','FontSize',20)
%    ylabel('Trials','FontSize',20)  

% %小波变换
% %%habutation 
% %average
% f_hab = [];
% wt_hab = [];
% for i= 1 : trial_hab
%     [wt,f,coi] = cwt(aftcon_heartbeathab(i,:),fps);
%     f_hab = [f_hab,f];
%     wt_hab(:,:,i) = wt;
% end
% f_hab = f_hab';
% mean_wthab = mean(wt_hab,3);
% t = (1 : size(abs(mean_wthab),2))/fps;
% figure;imagesc(t,f,abs(mean_wthab),[0.01 0.6]);
% hold on
%  title('Habitation','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20)  
% figure;imagesc(t,f,abs(mean_wthab),[0.01 0.6]);
% hold on
%  ylim([0 4])
%  title('Habitation','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% 
% %acquisation
% %average
% f_acq = [];
% wt_acq = [];
% for i= 1 : trial_acq
%     [wt,f,coi] = cwt(aftcon_heartbeatacq(i,:),fps);
%     f_acq = [f_acq,f];
%     wt_acq(:,:,i) = wt;
% end
% f_acq = f_acq';
% mean_wtacq = mean(wt_acq,3);
% t = (1 : size(abs(mean_wtacq),2))/fps;
% figure;imagesc(t,f,abs(mean_wtacq),[0.01 0.6]);
% hold on
%  title('Acquisition','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% figure;imagesc(t,f,abs(mean_wtacq),[0.01 0.6]);
% hold on
%  ylim([0 4])
%  title('Acquisition','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% 
% %test
% %average 1-5trial
% f_test = [];
% wt_test = [];
% for i= 1 : 5
%     [wt,f,coi] = cwt(aftcon_heartbeattest(i,:),fps);
%     f_test = [f_test,f];
%     wt_test(:,:,i) = wt;
% end
% f_test = f_test';
% mean_wttest = mean(wt_test,3);
% t = (1 : size(abs(mean_wttest),2))/fps;
% figure;imagesc(t,f,abs(mean_wttest),[0.01 0.6]);
% hold on
%  title('Test(tiral 1-5)','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% figure;imagesc(t,f,abs(mean_wttest),[0.01 0.6]);
% hold on
%  ylim([0 4])
%  title('Test(tiral 1-5)','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
%  
% %average 6-10th trial
% f_test = [];
% wt_test = [];
% for i= 1 : 5
%     [wt,f,coi] = cwt(aftcon_heartbeattest(i,:),fps);
%     f_test = [f_test,f];
%     wt_test(:,:,i) = wt;
% end
% f_test = f_test';
% mean_wttest = mean(wt_test,3);
% t = (1 : size(abs(mean_wttest),2))/fps;
% figure;imagesc(t,f,abs(mean_wttest),[0.01 0.6]);
% hold on
%  title('Test(tiral 6-10)','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% figure;imagesc(t,f,abs(mean_wttest),[0.01 0.6]);
% hold on
%  ylim([0 4])
%  title('Test(tiral 6-10)','FontSize',20)
%  xlabel('Time(s)','FontSize',20)
%  ylabel('Frequency(Hz)','FontSize',20) 
% 
% 
% 
% %进行小波变换
% % figure;w = cwt(heartbeat_hab,fps);
% % [wt,f,coi] = cwt(resh_heartbeathab(1,:),fps);
% % ylim([0.001 2]);
% % %smooth
% % heartbeat_habsm20=smooth(heartbeat_hab,20);
% % heartbeat_habsm300=smooth(heartbeat_hab,300);
% % figure;plot(heartbeat_habsm20-heartbeat_habsm300);
% %  ylim([0.01 1])
