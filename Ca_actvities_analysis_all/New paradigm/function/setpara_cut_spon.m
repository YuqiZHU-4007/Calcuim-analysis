function [fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara_cut_spon(behaviorpath)
ishuc=true;
fs.ca=0.4;%s0.28
fs.behavior=60;%1/s;1275/21

time.per_cycle=30;%s20
time.cs_start=12;%601;10
time.cs_end=16.8;%888;14.8
time.us_start=16.7;%793;13.2

frame.per_cycle = round(time.per_cycle/fs.ca);%20/s
frame.cs_start = round(time.cs_start/fs.ca+1);
frame.cs_end = round(time.cs_end/fs.ca+1);
frame.us_start = round(time.us_start/fs.ca+1);%34:793;

frameb.per_cycle =round(time.per_cycle*fs.behavior);
frameb.cs_start = round(time.cs_start*fs.behavior+1);
frameb.cs_end = round(time.cs_end*fs.behavior+1);
frameb.us_start =  round(time.us_start*fs.behavior+1);

% fs.ca=0.4;%s
% fs.behavior=60;%1/s
% 
% time.per_cycle=20;%s
% time.cs_start=10;%601
% time.cs_end=14.8;%888
% time.us_start=13.2;%793
% 
% frame.per_cycle = round(time.per_cycle/fs.ca);%20/s
% frame.cs_start = round(time.cs_start/fs.ca+1);
% frame.cs_end = round(time.cs_end/fs.ca+1);
% frame.us_start = round(time.us_start/fs.ca+1);%34:793;
% 
% frameb.per_cycle =round(time.per_cycle*fs.behavior);
% frameb.cs_start = round(time.cs_start*fs.behavior+1);
% frameb.cs_end = round(time.cs_end*fs.behavior+1);
% frameb.us_start =  round(time.us_start*fs.behavior+1);
%6:×Ô·¢
%15£ºhab
%6*4£ºacq
%6:test
%6:spon
%1:trial num;2:start trial;3:end trial
a=0;trial.spon_bef=[a min(1,a) a];
a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
trial.acq_block_num=8;
trial.acq_block_trial=3;
trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);

if isempty(behaviorpath)
    re_startpoint=[];startpoint=[];y_3sd=[]; 
else
    a=load(([behaviorpath]));
    if isfield(a,'re_startpoint_sd')
    re_startpoint=a.re_startpoint_sd;
    else
       re_startpoint=a.re_start_end_point(:,1:2); 
    end
    startpoint=a.startpoint_sd;
    y_3sd=a.y_3sd;%y_t!!!!!!!!!!
end
% [behaviorname,behaviorpath]=uigetfile(['G:\data-vamt-nonlearner\*' '.mat'],'behavior');
% [fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([behaviorpath behaviorname]);
% save([behaviorpath '\para'],'fs','time','frame','frameb','trial','re_startpoint','startpoint','y_3sd','ishuc','-v7.3');


