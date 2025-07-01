function [stimUS, ind_US,T_spon,T_US,l_spon,l_US]=getUSstim(trial,frame)
stimUS=zeros(1,frame.per_cycle*trial.total);ind_US=[];kk=1;
for ss=1:trial.acq_block_num
    for tt=1:trial.acq_block_trial
        ind=(trial.hab+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frame.per_cycle+frame.us_start;
        stimUS(ind:ind+max(round(frame.us_dur(ss)),1)-1)=1;%23
        ind_US(kk,1)=ind;kk=kk+1;
    end
end
find( stimUS==1);%figure,plot( stimUS);
%spon
T_spon={};T_US={};l_spon={};l_US={};
for ii=1:trial.acq_block_num+1
    switch ii
        case 1
            t=1:trial.hab*frame.per_cycle;
            l_spon{ii}='Base.';
        case mat2cell([2:trial.acq_block_num],1,ones(1,trial.acq_block_num-1))
            t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
            l_spon{ii}=['S Interv.' num2str(ii-1)]; 
        case trial.acq_block_num+1
            t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
            %l_spon{ii}='Test';
            l_spon{ii}=['S Interv.' num2str(ii-1)]; 
    end
    T_spon{ii}=t;
end
% us
for ii=1:trial.acq_block_num
    t=(trial.hab+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
        (trial.hab+ii*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
    T_US{ii}=t;
    l_US{ii}=['Session ' num2str(ii)];
end

