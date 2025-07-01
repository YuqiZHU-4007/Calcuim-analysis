function [P_value, dur_cond] = binom_test(re_startpoint,trial,frameb,bef_cond,dur_cond,aft_cond)
%%
%比较学习前后CS阶段的行为
%比较学习后CS和spontaneous的行为
%比较学习前后状态变化
%比较学习过程中CS阶段的行为
%比较学习过程中CS和spontaneous的行为
%%

%%%%%%%%%%%%%%hab all trials, test all trials%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!
if ~isempty(aft_cond.CS_shake)
    bef_cond_CS_15 = bef_cond.CS_shake(1:trial.hab(1),1)';
    aft_cond_CS_6 = aft_cond.CS_shake(1:trial.test(1),1)';
    prob_shake_bef_cond_CS_15 = mean(bef_cond_CS_15,2);%baseline ratio
    if prob_shake_bef_cond_CS_15 == 0
        prob_shake_bef_cond_CS_15 = 0.05;
    end
    event_shake_aft_cond_CS_6 = sum(aft_cond_CS_6,2);%tested ratio
    num_aft_cond_CS_6 = size(aft_cond_CS_6,2);
    P_value.CS_15_6 = 1-binocdf(event_shake_aft_cond_CS_6,num_aft_cond_CS_6,prob_shake_bef_cond_CS_15);
            else
      P_value.CS_15_6 =nan;
end
%compare CS with spontaneous!!!!!!!!!!!!!!!!!!!!!
if ~isempty(aft_cond.befCS_shake)
    aft_cond_befCS_6 = aft_cond.befCS_shake(1:trial.test(1),1)';
    aft_cond_CS_6 = aft_cond.CS_shake(1:trial.test(1),1)';
    prob_shake_aft_cond_befCS_6 = mean(aft_cond_befCS_6,2);
    if prob_shake_aft_cond_befCS_6==0
        prob_shake_aft_cond_befCS_6=0.05;
    end
    event_shake_aft_cond_CS_6 = sum(aft_cond_CS_6,2);
    num_aft_cond_CS_6 = size(aft_cond_CS_6,2);
    P_value.aftcond_befCS_CS_6 = 1-binocdf(event_shake_aft_cond_CS_6,num_aft_cond_CS_6,prob_shake_aft_cond_befCS_6);
    %[P_value.aftcond_befCS_CS_6_signrank,~,~] = signrank(aft_cond.befCS_shake(1:trial.test(1),1),aft_cond.CS_shake(1:trial.test(1),1),'alpha',0.05,'tail','left');
            else
      P_value.aftcond_befCS_CS_6 =nan;
end
%%%%compare whether states changed
if ~isempty(aft_cond.befCS_shake)
    aft_cond_befCS_6 = aft_cond.befCS_shake(1:trial.test(1),1)';
    bef_cond_befCS_15 = bef_cond.befCS_shake(1:trial.hab(1),1)';
    prob_shake_bef_cond_befCS_15 = mean(bef_cond_befCS_15,2);
    if prob_shake_bef_cond_befCS_15==0
        prob_shake_bef_cond_befCS_15=0.05;
    end
    event_shake_bef_cond_befCS_15 = sum(aft_cond_befCS_6,2);
    num_bef_cond_befCS_15 = size(aft_cond_befCS_6,2);
    P_value.states_6 = 1-binocdf(event_shake_bef_cond_befCS_15,num_bef_cond_befCS_15,prob_shake_bef_cond_befCS_15);
            else
       P_value.states_6 =nan;
end
%%
%compare acq(6 trials) with hab(all trials)!!!!!!!!!!!!!!!!!!!!!!!!!!
%compare CS before and during learning
P_value.dur_learnning_CS=[];c=1;slid_win=6;
for ii=[6 trial.hab(1)]
    b=0;
    bef_cond_CS192 = bef_cond.CS192_shake((trial.hab(1)-ii+1):trial.hab(1),1)';
    prob_shake_bef_cond_CS192 = mean(bef_cond_CS192,2);%baseline ratio
    if prob_shake_bef_cond_CS192==0
        prob_shake_bef_cond_CS192=0.05;
    end
    for jj=1:trial.acq(1)-slid_win+1
        dur_cond_CS = dur_cond.CS192_shake(jj:jj+slid_win-1,1)';
        event_shake_dur_cond_CS = sum(dur_cond_CS,2);
        num_dur_cond_CS = size(dur_cond_CS,2);
        b=b+1;
        P_value.dur_learnning_CS(jj,c) = 1-binocdf(event_shake_dur_cond_CS,num_dur_cond_CS,prob_shake_bef_cond_CS192);
    end
    c=c+1;
end
if ~isempty( P_value.dur_learnning_CS)
    P_value.dur_learnning_CS_6 = P_value.dur_learnning_CS(:,1);
    P_value.dur_learnning_CS_15 = P_value.dur_learnning_CS(:,2);
else
    P_value.dur_learnning_CS_6 = [];
    P_value.dur_learnning_CS_15 = [];
end

%%%%compare CS with spontaneous during learning
P_value.durcond_befCS_CS=[];
for ii = 1:trial.acq(1)-slid_win+1
    dur_cond_befCS = dur_cond.befCS_shake(ii:ii+slid_win-1,1)';
    dur_cond_CS = dur_cond.CS192_shake(ii:ii+slid_win-1,1)';
    prob_shake_dur_cond_befCS = mean(dur_cond_befCS,2);
    if prob_shake_dur_cond_befCS == 0
        prob_shake_dur_cond_befCS = 0.05;
    end
    event_shake_dur_cond_CS = sum(dur_cond_CS,2);
    num_dur_cond_CS = size(dur_cond_CS,2);
    P_value.durcond_befCS_CS(ii,1) = 1-binocdf(event_shake_dur_cond_CS,num_dur_cond_CS,prob_shake_dur_cond_befCS);
end
%%%%compare spontaneous before learning with spontaneous during learning
if ~isempty(dur_cond.befCS_shake)
    for ii = 1:trial.acq(1)-slid_win+1
        bef_cond_befCS_15 = bef_cond.befCS192_shake(1:trial.hab(1),1)';
        dur_cond_befCS = dur_cond.befCS_shake(ii:ii+slid_win-1,1)';
        prob_shake_bef_cond_befCS_15 = mean(bef_cond_befCS_15,2);
   
    if prob_shake_bef_cond_befCS_15==0
        prob_shake_bef_cond_befCS_15=0.05;
    end
    event_shake_bef_cond_befCS_15 = sum(dur_cond_befCS,2);
    num_bef_cond_befCS_15 = size(dur_cond_befCS,2);
    P_value.durcond_states_6(ii,1) = 1-binocdf(event_shake_bef_cond_befCS_15,num_bef_cond_befCS_15,prob_shake_bef_cond_befCS_15); 
    end
else
    P_value.durcond_states_6 =nan;
end

%%
%%%%%%%%%%%%%hab 15 trials, test 3 trials%%%%%%%%%%%
%compare hab with test during CS
if ~isempty(aft_cond.CS_shake)
    bef_cond_CS_15 = bef_cond.CS_shake(1:(trial.hab(1)),1)';
    aft_cond_CS_3 = aft_cond.CS_shake(1:3,1)';
    prob_shake_bef_cond_CS_15 = mean(bef_cond_CS_15,2);
    if prob_shake_bef_cond_CS_15 == 0
        prob_shake_bef_cond_CS_15 = 0.05;
    end
    event_shake_aft_cond_CS_3 = sum(aft_cond_CS_3,2);
    num_aft_cond_CS_3 = size(aft_cond_CS_3,2);
    P_value.CS_15_3 = 1-binocdf(event_shake_aft_cond_CS_3,num_aft_cond_CS_3,prob_shake_bef_cond_CS_15);
else
    P_value.CS_15_3=nan;
end
%compare CS with spontaneous
if ~isempty(aft_cond.befCS_shake)
    aft_cond_befCS_3 = aft_cond.befCS_shake(1:3,1)';
    prob_shake_aft_cond_befCS_3 = mean(aft_cond_befCS_3,2)/2;
    if prob_shake_aft_cond_befCS_3==0
        prob_shake_aft_cond_befCS_3=0.05;
    end
    event_shake_aft_cond_CS_3 = sum(aft_cond_CS_3,2);
    num_aft_cond_CS_3 = size(aft_cond_CS_3,2);
    P_value.aftcond_befCS_CS_3 = 1-binocdf(event_shake_aft_cond_CS_3,num_aft_cond_CS_3,prob_shake_aft_cond_befCS_3);
    else
     P_value.aftcond_befCS_CS_3=nan;
end
%%%%compare whether states changed
if ~isempty(aft_cond.befCS_shake)
    aft_cond_befCS_3 = aft_cond.befCS_shake(1:3,1)';
    bef_cond_befCS_15 = bef_cond.befCS_shake(1:trial.hab(1),1)';
    prob_shake_bef_cond_befCS_15 = mean(bef_cond_befCS_15,2);
    if prob_shake_bef_cond_befCS_15==0
        prob_shake_bef_cond_befCS_15=0.05;
    end
    event_shake_bef_cond_befCS_15 = sum(bef_cond_befCS_15,2);
    num_bef_cond_befCS_15 = size(bef_cond_befCS_15,2);
    P_value.states_3 = 1-binocdf(event_shake_bef_cond_befCS_15,num_bef_cond_befCS_15,prob_shake_bef_cond_befCS_15);
        else
      P_value.states_3=nan;
end
%%%compare CS with spontaneous during learning
%compare CS with before 192 points
P_value.durcond_befCS192_CS=[];
for ii = 1:trial.acq_block_trial
    dur_cond_befCS192 = dur_cond.befCS192_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    dur_cond_CS = dur_cond.CS192_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    prob_shake_dur_cond_befCS192 = mean(dur_cond_befCS192,2);
    if prob_shake_dur_cond_befCS192 == 0
        prob_shake_dur_cond_befCS192 = 0.05;
    end
    event_shake_dur_cond_CS = sum(dur_cond_CS,2);
    num_dur_cond_CS = size(dur_cond_CS,2);
    PValue_dur_befCS192_CS = 1-binocdf(event_shake_dur_cond_CS,num_dur_cond_CS,prob_shake_dur_cond_befCS192);
    P_value.durcond_befCS192_CS = [P_value.durcond_befCS192_CS;PValue_dur_befCS192_CS];
end

%compare CS with before 384 points
P_value.durcond_befCS384_CS=[];
for ii = 1:trial.acq_block_trial
    dur_cond_befCS384 = dur_cond.befCS384_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    dur_cond_CS = dur_cond.CS192_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    prob_shake_dur_cond_befCS384 = mean(dur_cond_befCS384,2);
    if prob_shake_dur_cond_befCS384 == 0
        prob_shake_dur_cond_befCS384 = 0.05;
    end
    event_shake_dur_cond_CS = sum(dur_cond_CS,2);
    num_dur_cond_CS = size(dur_cond_CS,2);
    PValue_dur_befCS384_CS = 1-binocdf(event_shake_dur_cond_CS,num_dur_cond_CS,prob_shake_dur_cond_befCS384);
    P_value.durcond_befCS384_CS = [P_value.durcond_befCS384_CS;PValue_dur_befCS384_CS];
end

%compare CS with before 576 points
P_value.durcond_befCS576_CS=[];
for ii = 1:trial.acq_block_trial
    dur_cond_befCS576 = dur_cond.befCS576_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    dur_cond_CS = dur_cond.CS192_shake((1+trial.acq_block_num*(ii-1)):(trial.acq_block_num+trial.acq_block_num*(ii-1)),1)';
    prob_shake_dur_cond_befCS576 = mean(dur_cond_befCS576,2);
    if prob_shake_dur_cond_befCS576 == 0
        prob_shake_dur_cond_befCS576 = 0.05;
    end
    event_shake_dur_cond_CS = sum(dur_cond_CS,2);
    num_dur_cond_CS = size(dur_cond_CS,2);
    PValue_dur_befCS576_CS = 1-binocdf(event_shake_dur_cond_CS,num_dur_cond_CS,prob_shake_dur_cond_befCS576);
    P_value.durcond_befCS576_CS = [P_value.durcond_befCS576_CS;PValue_dur_befCS576_CS];
end


end


