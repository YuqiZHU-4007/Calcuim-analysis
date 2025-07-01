%updata from F:\DUlab\FC
%analyse\Ca_actvities_analyze\Behavior\myanalyse_tail_gonogo\some_plot_20200415\plot_tail_mov_US_2_1.m
%updata: setpara_testUS & getUSstim
%zyq,20200430

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

%% setpara
[frameb,framec,trial,fs]=setpara_testUS(4);

% [fs,time,framec,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
% frameb.us_dur=4.8*ones(1,trial.acq_block_num)*fs.behavior;framec.us_dur=4.8*ones(1,trial.acq_block_num)/fs.ca;
[stimCS_b, ind_CS_b,T_spon,T_US,l_spon,l_us]=getUSstim(trial,frameb);%figure,plot( stimUS);
[stimCS_c, ind_CS_c,~,~,~,~]=getUSstim(trial,framec);%figure,plot( stimUS);
 trial.spon=5*60/30;trial.total=6;
%% alltheta
controlpoint=8;
nframes=length(movlabel);
basevec = [selx(controlpoint) sely(controlpoint)] - [selx(1) sely(1)];
movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(1); sely(1)], 1, size(alltailpos,3));
% distance of control points
dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));
% theta, anti-clockwise as positive
theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
theta(theta > 180) = theta(theta > 180) - 360;
theta = -theta;
theta00 = 0;
alltheta_r = theta00 * ones(nframes, 1);
alltheta_r(movind) = theta;

alltheta=alltheta_r;

cut_thr=[0.5 5];
aa=abs(diff([alltheta(1);alltheta]));alltheta(find(aa<cut_thr(1) & abs(alltheta)<cut_thr(2)))=0;
%alltheta=alltheta((trial.spon)*frameb.per_cycle+1:(trial.spon+trial.acq_block_trial)*frameb.per_cycle);
alltheta=alltheta(1:trial.total*frameb.per_cycle);
alltheta_n=normalize(alltheta,'zscore');
d_alltheta=diff([alltheta(1);alltheta]);
d_alltheta_n=diff([alltheta_n(1);alltheta_n]);

%% plot
T_b=[1:trial.total*frameb.per_cycle]/fs.b/60;T_c=[1:trial.total*framec.per_cycle];
figure('position',[50,250,1700,450]),plot(T_b,stimCS_b*max(alltheta_n(:)),'r','linewidth',1.5);hold on;
plot(T_b,alltheta_n);xlim([min(T_b) max(T_b)]);xlabel('Time(min)','fontsize',16);ylabel('Tail Angle');set(gca,'fontsize',20,'linewidth',1.5);
%save([path 'Results_of_alltheta.mat'],'Results_of_alltheta','controlpoint','alltheta','intermovinterv_spon','intermovinterv_US','merge_ind_US','cut_thr');
title(path(end-14:end))
% n=1:24;
% c=hsv(length(n));
% [h, incre] = sepplot((1:15611)*0.02, aa(:,n), c,0.8);
% xlim([1 15611]*0.02);set(gca,'fontsize',14)
