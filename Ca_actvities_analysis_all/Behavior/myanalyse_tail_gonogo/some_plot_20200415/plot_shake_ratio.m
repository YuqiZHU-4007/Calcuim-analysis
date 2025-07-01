clc;clear all;close all;
global fs
global frameb
global trial

load('I:\1.Test US\2.Tail free¡ª¡ªData from 117\20210101\fish1\Results_of_alltheta.mat');
for ii=1
    [fs,time,~,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=10;
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
end

ii=33;Name(ii)
b=V_trial_r_preCS{3}(1:6,ii);bb=V_trial_r_preCS{3}(31:end,ii);
c=V_trial_r_CS{3}(1:6,ii);cc=V_trial_r_CS{3}(31:end,ii);
pre=[b,c]';aft=[bb,cc]';
 P_ranksum_Vcscs(ii)
P_ranksum_Vprecscs_post(ii)
P_ranksum_Vsponspon(ii)

startpoint=env.env_loc;
re_startpoint(:,1)=fix(startpoint/frameb.per_cycle)+1;
re_startpoint(:,2)=mod(startpoint,frameb.per_cycle);
a=0;trial.spon_bef=[a min(1,a) a];
a=6;trial.hab = [a trial.spon_bef(3)+1 trial.spon_bef(3)+a];
trial.acq_block_num=8;
trial.acq_block_trial=3;
trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
a=6;trial.test =[a trial.acq(3)+1 trial.acq(3)+a];
a=0;trial.spon_aft=[a trial.test(3)+1 trial.test(3)+a];
trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);

%% 
isbinary=1;
[p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(alltheta(T_non_spon)',re_startpoint,[],[],[],isbinary);
P_bino_movratio_cscs_b=p_value.CS_15_6;
P_bino_movratio_precscs_post_b=p_value.aftcond_befCS_CS_6;
P_bino_movratio_sponspon_b=p_value.states_6;
%%
isbinary=0;
[p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(alltheta(T_non_spon)',re_startpoint,[],[],[],isbinary);
pre_con=[bef_cond.befCS_shake/2 bef_cond.CS_shake];
aft_con=[aft_cond.befCS_shake/2 aft_cond.CS_shake];
P_bino_movratio_cscs=p_value.CS_15_6
P_bino_movratio_precscs_post=p_value.aftcond_befCS_CS_6
P_bino_movratio_sponspon=p_value.states_6

