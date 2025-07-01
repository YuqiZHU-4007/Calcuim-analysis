function event=getevent_Trianglefitting_strict_for_new_paradigm(a,baseline,maxthres,rasing_width,decling_width,frame,trial,base)
%rasing_width:从起始点到最高点的间隔（包括最高点，包括起始点）;区间【最小 最大】
%decling_width：从最高点到结束点（包括最高点，包括结束点）;区间【最小 最大】
%baseline：认为用反应的thres
%maxthres：反应最高点的thres
%一般认为起始点>1sd（baseline）且最高点>2sd（maxthres）,是反应的基本条件，再看反应的width（上升和下降）排除噪声
%20180925修改版：考虑波动导致峰内下降，错误判断最高点，使得快上升慢下降的条件不满足（求最高点时smoothdata）
%a为列向量
%20181016修改：因为trial的格式变化，把用到trial的地方改成新的格式

% ii=4;jj=1; a=activities_preCS_dfdf_aftcorrect{ii}(:,jj);
% aa=reshape(a,frame.per_cycle,trial.total)';
% base=getbaseline_cut3sd(aa,1:frame.cs_start-1,trial);
% m=base.rep_m;sd=base.rep_sd;
% % m=mean(aa(:,1:frame.cs_start-1),2);m=kron(m,ones(frame.per_cycle,1));
% % sd=std(aa(:,1:frame.cs_start-1),[],2);sd=kron(sd,ones(frame.per_cycle,1));
% baseline=m+2*sd;maxthres=m+3*sd;maxthres2=m+2*sd;
% rasing_width=[1 100];rasing_width2=2;
% decling_width=[4 100];decling_width2=5;%decling_width+rasing_width-1<8

event=[];
ind=find(a>baseline);
width_max=ceil(ind/frame.per_cycle)*frame.per_cycle-ind;
ind(find(width_max==0))=[];width_max(find(width_max==0))=[];
xx=1;event.end_ind=[];event.end_ind(xx)=0;event.rep_ind=zeros(size(a));event.ind=[];
sm_a=smoothdata(a,'movmean',3);
for kk=1:length(ind)
    if ind(kk)>event.end_ind(end) && a(ind(kk))>=a(max(ind(kk)-1,1))
        for zz=width_max(kk):-1:1
            if sum(a(ind(kk):ind(kk)+zz-1)>=baseline(ind(kk):ind(kk)+zz-1))>=zz-1
                [ma,max_ind]=max(a(ind(kk):ind(kk)+zz-1));
                rasing_dur=(max_ind+ind(kk)-1)-ind(kk)+1;
                decline_dur=zz-(max_ind)+1;
                if max(sm_a(ind(kk):ind(kk)+zz-1))>=maxthres(ind(kk)) && rasing_dur>=rasing_width(1) && rasing_dur<=rasing_width(2)...
                        &&  decline_dur>=decling_width(1) && decline_dur<=decling_width(2) ...
                        && mod(ind(kk),frame.per_cycle)~=1%不是trial的第一个点
                    if rasing_dur==1 && sum(a(max_ind+ind(kk)-1:max_ind+ind(kk)-1+1)>=maxthres(ind(kk)))==2
                        event.ind(xx)=ind(kk);event.end_ind(xx)=ind(kk)+zz-1;event.dur(xx)=zz+1;
                        event.rasing_dur(xx)=rasing_dur; event.decline_dur(xx)=decline_dur;
                        event.max_ind(xx)=max_ind+ind(kk)-1;
                        event.rep_ind(ind(kk)-1)=1;
                        xx=xx+1;break;
                    elseif rasing_dur~=1
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
% %做onset精细判断
%event.ind1=event.ind;
if ~isempty(event.ind)
%     for ii=1:length(event.ind)
%         onset=event.ind(ii);
%         if a(onset)>baseline(onset)
%             for zz=min(8,mod(onset,frame.per_cycle)-1):-1:1 %找前面一段是上升相
%                 m=mean(diff(base.output_cut3sd(ceil(onset/frame.per_cycle),1:frame.cs_start-1)));%波动情况
%                 sd=std(diff(base.output_cut3sd(ceil(onset/frame.per_cycle),1:frame.cs_start-1)));
%                 if sum(diff(sm_a(onset-zz:onset))>0)>=zz-1 && diff(a(onset-zz:onset-zz+1))>m+2*sd...
%                         && sum(diff(a(onset:onset+event.rasing_dur(ii)-1))>0)==event.rasing_dur(ii)-1 ... %本身属于持续上升相
%                         && ceil((onset-zz)/frame.per_cycle)==ceil(onset/frame.per_cycle)
%                     
%                     event.ind(ii)=onset-zz;
%                     %event.ind1(ii)=onset-zz;
%                     %onset-zz;
%                     break;
%                 end
%             end
%         end
%     end
%     
    %acq阶段trial的前4个点有反应不算，只适用于连续成像时候，反应还没下去的情况
    event.end_ind(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
    event.dur(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
    event.rasing_dur(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
    event.decline_dur(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
    event.max_ind(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
    event.rep_ind(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=0;
    event.ind(find(ceil(event.ind/frame.per_cycle)>trial.hab(3) & ceil(event.ind/frame.per_cycle)<=trial.acq(3) & mod(event.ind,frame.per_cycle)<=8))=[];
end
% for ii=1:length(event.ind)
%     onset=event.ind(ii);
%     [mi,min_ind]=min(a(onset-8:event.max_ind(ii)));
%     peakheight=a(event.max_ind(ii))-mi;
%     ind=find(sm_a(onset-8:event.max_ind(ii))>0.3*peakheight);
%     if ~isempty(ind)
%         event.ind1(ii)=ind(1)+onset-8-1-1;
%     end
% end
%
% figure,plot(a);hold on;
% scatter(event.ind,a(event.ind),'k');hold on;
% scatter(event.ind1,a(event.ind1),'r');hold on;
% %
% %figure,plot(activities{ii}(:,jj));
% figure,plot(a);hold on;plot(baseline);hold on;plot(maxthres);hold on;
% scatter(event.ind,a(event.ind),'r');hold on;
% scatter(event.max_ind,a(event.max_ind),'y');hold on;
% scatter(event.end_ind,a(event.end_ind),'o');hold on;
% plot(smoothdata(a,'movmean',3),'k');