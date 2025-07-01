function [isevent]=isevent_20190508(xx,ref_win,ind,trial)

smoothwin=3;
rasing_width=[1 100];
decling_width=[4 100];

indc=ind(1:end-1);
%xxs=smoothdata(xx([ref_win indc]),'movmean',smoothwin);
xxs=xx([ref_win indc]);
dxxs= detrend(xxs(ref_win-ref_win(1)+1));
thres=mean(dxxs)+3*std(dxxs);
% base=getbaseline_cut3sd_strict_for_new_paradigm(xx',ref_win,trial);
% thres=base.rep_m+3*base.rep_sd;
[m,indm]=max(xxs(length(ref_win)+1:end));
k=[m-xxs(length(ref_win)+1:end)]'./abs([1:length(indc)]-indm);%m-xxs([ind(1):indm+ind(1)-1]-ref_win(1)+1)./(indm+ind(1)-1-[ind(1):indm+ind(1)-1])';
k(isinf(k))=[];
if ~isempty(k)
    [sm,indr]=max(k);
    [ssm,indd]=min(k);
    rasing_dur=indm-indr+1;
    decline_dur=indd-indm+1;
    %                 figure,plot(xx);hold on;plot(thres*ones(size(xs)));hold on;
    %                 scatter([inds+ind(1)-1 indm+ind(1)-1],[xs(inds+ind(1)-1) xs(indm+ind(1)-1)]);hold on;
    %                 plot([1:50],thres/ length(ind)*[1:50]+xs(indm+ind(1)-1)-(inds+ind(1)-1)*thres/ length(ind))    
    if  m> thres  && sm>= (thres-min(xxs(indc-ref_win(1)+1)))/length(ind) && ~isinf(sm) ...
            && ssm>= -(thres-min(xxs(indc-ref_win(1)+1)))/length(ind) && ~isinf(ssm) ...
             &&  rasing_dur>=rasing_width(1) && rasing_dur<=rasing_width(2) ...
              &&  decline_dur>=decling_width(1) && decline_dur<=decling_width(2) 
        isevent=true;
    else
       isevent=false;
    end
else
    isevent=false;
end

