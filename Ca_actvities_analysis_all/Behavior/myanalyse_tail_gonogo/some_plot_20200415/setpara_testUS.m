%% setpara
function [frameb,framec,trial,fs]=setpara_testUS(exp_id)
if nargin==0
    exp_id=1;
end
frameb=struct;
framec=struct;
switch exp_id
    case 1  %5*6,ITI:20
        trial.hab=15;
        trial.acq_block_trial=6;
        trial.acq_block_num=5;
        trial.acq_block_interval=15;
        trial.test=15;
        frameb.per_cycle=20*60;
        frameb.us_start=round(13.2*60+1);
        frameb.us_dur=[0.001 0.001 0.01 0.01 0.1]*60;
    case 2 %8*3,ITI:20
        trial.hab=15;
        trial.acq_block_trial=3;
        trial.acq_block_num=8;
        trial.acq_block_interval=15;
        trial.test=15;
        frameb.per_cycle=20*60;
        frameb.us_start=round(13.2*60+1);
        frameb.us_dur=0.001*ones(1,trial.acq_block_num)*60;
    case 3 %8*3,ITI:30
        trial.hab=10;
        trial.acq_block_trial=3;
        trial.acq_block_num=8;
        trial.acq_block_interval=10;
        trial.test=10;
        frameb.per_cycle=30*60;
        frameb.us_start=round(13.2*60+1);
        frameb.us_dur=0.001*ones(1,trial.acq_block_num)*60;
    case 4 %only cs
        trial.hab=0;
        trial.acq_block_trial=6;
        trial.acq_block_num=1;
        trial.acq_block_interval=0;
        trial.test=0;
        fs.b=60;fs.c=2.5;
        frameb.per_cycle=100*fs.b;
        frameb.us_start=round(10*fs.b+1);
        frameb.us_dur=4.8*ones(1,trial.acq_block_num)*fs.b;
        framec.per_cycle=100*fs.c;
        framec.us_start=round(10*fs.c+1);
        framec.us_dur=4.8*ones(1,trial.acq_block_num)*fs.c;
end
trial.total=trial.hab+trial.test+trial.acq_block_trial*trial.acq_block_num+trial.acq_block_interval*(trial.acq_block_num-1);