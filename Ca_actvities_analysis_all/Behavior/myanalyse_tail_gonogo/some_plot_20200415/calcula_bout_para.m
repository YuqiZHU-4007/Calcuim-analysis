%Bout tail max amplitude
%Bout mean tail angle
%Max tail change
%Mean tail change
%Sum tail change
%Bout duration
%Beat number
%Beat number2
function [P,alltheta_out]=calcula_bout_para(alltheta,env_loc,env_end)
alltheta=normalize(alltheta,'zscore');
alltheta_out=nan(length(env_loc),500);
P=struct;
for ii=1:length(env_loc)
    if ~isnan(env_loc(ii))
        x=alltheta(env_loc(ii):env_end(ii));
        P(ii).max_amp=max(x);
        P(ii).mean_amp=mean(abs(x));
        P(ii).max_tail_change=max(diff(x));
        P(ii).mean_tail_change=nanmean(abs(diff(x)));
        P(ii).sum_tail_change=sum(abs(diff(x)));
        P(ii).dur=env_end(ii)-env_loc(ii);
        a=length(find(x>0));b=length(find(x<0)); 
        P(ii).beat_num=a-b;
        P(ii).beat_num_2=length(find(abs(x)<=0.05));
        alltheta_out(ii,1:length(x))=x;
    else
        P(ii).max_amp=nan;
        P(ii).mean_amp=nan;
        P(ii).max_tail_change=nan;
        P(ii).mean_tail_change=nan;
        P(ii).sum_tail_change=nan;
        P(ii).dur=nan;
        P(ii).beat_num=nan;
        P(ii).beat_num_2=nan;
    end
end


