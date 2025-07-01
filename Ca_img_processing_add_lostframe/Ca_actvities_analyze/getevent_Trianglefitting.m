function event=getevent_Trianglefitting(a,baseline,maxthres,rasing_width,decling_width,frame)
%rasing_width:从起始点到最高点的间隔（包括最高点，包括起始点）
%decling_width：从最高点到结束点（包括最高点，包括结束点）
%baseline：认为用反应的thres
%maxthres：反应最高点的thres
%一般认为起始点>1sd（baseline）且最高点>2sd（maxthres）,是反应的基本条件，再看反应的width（上升和下降）排除噪声
%20180925修改版：考虑波动导致峰内下降，错误判断最高点，使得快上升慢下降的条件不满足（求最高点时smoothdata）
%a为列向量

% ii=2;jj=16; a=activities_preCS_dfdf_aftcorrect{ii}(:,jj);
% aa=reshape(a,frame.per_cycle,trial.total)';
% base=getbaseline_cut3sd(aa,1:frame.cs_start-1);
% m=base.rep_m;sd=base.rep_sd;
% % m=mean(aa(:,1:frame.cs_start-1),2);m=kron(m,ones(frame.per_cycle,1));
% % sd=std(aa(:,1:frame.cs_start-1),[],2);sd=kron(sd,ones(frame.per_cycle,1));
% baseline=m+1*sd;maxthres=m+3*sd;maxthres2=m+2*sd;
% rasing_width=1;rasing_width2=2;
% decling_width=4;decling_width2=5;%decling_width+rasing_width-1<8

event=[];
ind=find(a>baseline);
width_max=ceil(ind/frame.per_cycle)*frame.per_cycle-ind;
ind(find(width_max==0))=[];width_max(find(width_max==0))=[];
xx=1;event.end_ind=[];event.end_ind(xx)=0;event.rep_ind=zeros(size(a));
sm_a=smoothdata(a,'movmean',3);
for kk=1:length(ind)
    if ind(kk)>=event.end_ind(end) && a(ind(kk))>=a(max(ind(kk)-1,1))
        for zz=width_max(kk):-1:1
            if sum(sm_a(ind(kk):ind(kk)+zz-1)>=baseline(ind(kk):ind(kk)+zz-1))>=zz-1
                [ma,max_ind]=max(a(ind(kk):ind(kk)+zz-1));
                rasing_dur=(max_ind+ind(kk)-1)-ind(kk)+1;
                decline_dur=zz-(max_ind)+1;
                %                 if sum(smoothdata(a(ind(kk):ind(kk)+zz-1),2,'movmean',3)>=baseline(ind(kk):ind(kk)+zz-1))>=zz-3 ...%考虑acq US前面一段
                %                 else
                %                 end
                if max(sm_a(ind(kk):ind(kk)+zz-1))>=maxthres(ind(kk)) && rasing_dur>=rasing_width &&  decline_dur>=decling_width && mod(ind(kk),frame.per_cycle)~=1%不是trial的第一个点
                    if rasing_dur==1
                        event.ind(xx)=ind(kk)-1;event.end_ind(xx)=ind(kk)+zz-1;event.dur(xx)=zz+1;
                        event.rasing_dur(xx)=rasing_dur; event.decline_dur(xx)=decline_dur;
                        event.max_ind(xx)=max_ind+ind(kk)-1;
                        event.rep_ind(ind(kk)-1)=1;
                        xx=xx+1;break;
                    else
                        event.ind(xx)=ind(kk);event.end_ind(xx)=ind(kk)+zz-1;event.dur(xx)=zz;
                        event.rasing_dur(xx)=rasing_dur; event.decline_dur(xx)=decline_dur;
                        event.max_ind(xx)=max_ind+ind(kk)-1;
                        event.rep_ind(ind(kk))=1;
                        xx=xx+1;break;
                    end
                    %                 elseif ma>=maxthres2(ind(kk)) && rasing_dur>=rasing_width2 &&  decline_dur>=decling_width2
                    %                     event.ind(xx)=ind(kk);event.end_ind(xx)=ind(kk)+zz-1;event.dur(xx)=zz;
                    %                     event.rasing_dur(xx)=rasing_dur; event.decline_dur(xx)=decline_dur;
                    %                     event.max_ind(xx)=max_ind+ind(kk)-1;
                    %                     xx=xx+1;break;261
                end
            end
        end
    end
end

% %figure,plot(activities{ii}(:,jj));
% figure,plot(a);hold on;plot(baseline);hold on;plot(maxthres);hold on;
% scatter(event.ind,a(event.ind),'r');hold on;
% scatter(event.max_ind,a(event.max_ind),'y');hold on;
% scatter(event.end_ind,a(event.end_ind),'o');hold on;
% plot(smoothdata(a,'movmean',3),'k');