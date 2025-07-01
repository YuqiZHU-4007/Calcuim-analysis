function [fs,time,frame,frameb,trial,T_non_spon_be,T_spon_be,T_non_spon_ca,T_spon_ca]=setpara_spon()
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
%frame.us_dur=0.01*ones(1,trial.acq_block_num)*fs.behavior;

%1:trial num;2:start trial;3:end trial
a=10;trial.spon_bef=[a min(1,a) a];
a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
trial.acq_block_num=8;
trial.acq_block_trial=3;
trial.acq_block_interval=10;
trial.acq = [trial.acq_block_trial*trial.acq_block_num ...
    trial.hab(3)+1 ...
    trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval];
a=6;trial.test =[a trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1 ...
    trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+a];
a=10;trial.spon_aft=[a ,trial.test(3)+1,trial.test(3)+a];

trial.total=trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1)+trial.acq_block_interval*(trial.acq_block_num);

T_non_spon_be=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frameb.per_cycle+1:(trial.hab(3))*frameb.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(2)-1)*frameb.per_cycle+1:(trial.test(3))*frameb.per_cycle;
    end
    T_non_spon_be=[T_non_spon_be,t];
end
T_spon_be={};
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.spon_bef(2)-1)*frameb.per_cycle+1:(trial.spon_bef(3))*frameb.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frameb.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frameb.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(3))*frameb.per_cycle+1:(trial.total)*frameb.per_cycle;
    end
    T_spon_be{ii}=t;
end

T_non_spon_ca=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3))*frame.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle;
    end
    T_non_spon_ca=[T_non_spon_ca,t];
end
T_spon_ca={};
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.spon_bef(2)-1)*frame.per_cycle+1:(trial.spon_bef(3))*frame.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
%         case trial.acq_block_num+1
%             t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
%                 (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
        case trial.acq_block_num+2
            t=(trial.test(3))*frame.per_cycle+1:(trial.total)*frame.per_cycle;
    end
    T_spon_ca{ii}=t;
end
