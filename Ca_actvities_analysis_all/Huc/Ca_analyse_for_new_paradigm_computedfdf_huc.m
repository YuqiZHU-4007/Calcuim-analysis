function [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf_huc(activities_new)
global fs
global time
global frame
global frameb
global trial

set(0,'defaultfigurecolor','w');

%set para.
%[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);

activities=[];
if ~isempty(activities_new(:))& size(activities_new,1)<trial.total*frame.per_cycle
    activities=[activities_new(1,:);activities_new;activities_new(size(activities_new,1),:);];
elseif size(activities_new,1)==trial.total*frame.per_cycle
    activities=activities_new;
elseif isempty(activities_new(:))  
    error('activities_new is empty');
end

if size(activities,2)<trial.total*frame.per_cycle
    error('frame number is wrong');
end

%detrend
% for ii = 1:size(activities,1)
%     for jj=1:size(activities{ii},2)
%     activities_detrend{ii}(:,jj)=detrend(activities{ii}(:,jj));
%     end
% end
% %smooth and 补齐到3000
% if size(activities,1)<trial.total*frame.per_cycle
%     %activities{ii,1}=[activities{ii,1}(1,:);activities{ii,1};];
%     a=ones(frame.per_cycle*trial.total,size(activities,2)).*repmat(activities(end,:),frame.per_cycle*trial.total,1);
%     a(1:size(activities,1),:)=activities;
%     activities_total=a;
% else
%     activities_total=activities;
% end
%figure,plot_rawtrace_trials(activities(1:2100,2),[],fs,frame,trial,[],1);title('aft process')
%dfdf,用每个trial前段时间的基线
activities_preCS_baseline=[];activities_preCS_dfdf=[];
a=activities; 
for jj=1:size(a,2)
    aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
    base=getbaseline_cut3sd_strict_for_new_paradigm(aa,frame.cs_start-12:frame.cs_start-1,trial);
    activities_preCS_baseline(:,jj)=base.rep_m;
    activities_preCS_dfdf(:,jj)=(a(:,jj)-base.rep_m)./base.rep_m;
end


%去下降点
activities_preCS_dfdf_aftcorrect=[];
a=activities_preCS_dfdf;
for jj=1:size(a,2)
    aa=a(:,jj); sort_aa=sort(aa); sort_aa= sort_aa(1:ceil(length(sort_aa)*0.975));
    m=mean(sort_aa);sd=std(sort_aa);
    
    ind1=find((aa(2:end-2)-aa(1:end-3))<-(m+2*sd))+1;
    ind2=find((aa(2:end-2)-aa(4:end))<-(m+sd))+1;
    ind_cut1_1=intersect(ind1,ind2);
    ind1=find((aa(3:end-1)-aa(1:end-3))<-(m+2*sd))+2;
    ind2=find((aa(2:end-1)-aa(3:end))<-(m+sd))+1;
    ind_cut1_2=intersect(ind1,ind2);
    ind_cut2=intersect(ind_cut1_1,ind_cut1_2-1);ind_cut2=union(ind_cut2,ind_cut2+1);
    
    ind1=find((aa(2:end-2)-aa(1:end-3))<-(m+sd))+1;
    ind2=find((aa(2:end-2)-aa(4:end))<-(m+2*sd))+1;
    ind_cut1_1=intersect(ind1,ind2);
    ind1=find((aa(3:end-1)-aa(1:end-3))<-(m+sd))+2;
    ind2=find((aa(2:end-1)-aa(3:end))<-(m+2*sd))+1;
    ind_cut1_2=intersect(ind1,ind2);
    ind_cut3=intersect(ind_cut1_1,ind_cut1_2-1);ind_cut3=union(ind_cut3,ind_cut3+1);
    ind1=find((aa(2:end-1)-aa(1:end-2))<-(m+sd))+1;
    ind2=find((aa(2:end-1)-aa(3:end))<-(m+2*sd))+1;
    ind_cut1_1=intersect(ind1,ind2);
    ind1=find((aa(2:end-1)-aa(1:end-2))<-(m+2*sd))+1;
    ind2=find((aa(2:end-1)-aa(3:end))<-(m+sd))+1;
    ind_cut1_2=intersect(ind1,ind2);
    ind_cut1=union(ind_cut1_1,ind_cut1_2);
    
    ind_cut=union(ind_cut1,ind_cut2);ind_cut=union(ind_cut,ind_cut3);
    
    aa_cut=myinterp(1:length(aa),aa',ind_cut,'max');
    activities_preCS_dfdf_aftcorrect(:,jj)=aa_cut;
    %         figure,plot(aa,'b');hold on;plot(aa_cut,'r');hold on;
    %         scatter(ind_cut,aa(ind_cut),[],'k');
end

%计算onset
activities_event_preCS=[];
a=activities_preCS_dfdf_aftcorrect;
for jj=1:size(a,2)
    event=getevent(a(:,jj),trial,frame);
    activities_event_preCS{jj,1}=event;
end


% %计算曲线下积分
% if trial.hab(1)<=trial.test(1)
%     indcut1=0;
%     indcut2=trial.test(1)-trial.hab(1);
% else
%     indcut1=trial.hab(1)-trial.test(1);
%     indcut2=0;
% end
% a=activities_preCS_dfdf_aftcorrect;
% area.CS_hab_tst=zeros(2,size(a,2));
% area.CS_acq_block=zeros(trial.acq_block_num,size(a,2));
% area.CS_acq_1st_block=zeros(trial.acq_block_num,size(a,2));
% for jj=1:size(a,2)
%     rawacti=[];
%     aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
%     rawacti(:,1)=mean(aa(trial.hab(2):trial.hab(3)-indcut1,:),1);
%     rawacti(:,2)=mean(aa(trial.test(2):trial.test(3)-indcut2,:),1);
%     area.CS_hab_tst(:,jj) = sum(rawacti(frame.cs_start:(frame.cs_end-1),:));%trapz(frame.cs_start:(frame.cs_end-1),rawacti(frame.cs_start:(frame.cs_end-1),:));
%     areaa=[];
%     for ii=1:trial.acq_block_num
%         ind1=trial.acq_block_trial*(ii-1)+1+trial.hab(3);
%         ind2=trial.acq_block_trial*(ii-1)+1+trial.hab(3)+trial.acq_block_trial-1;
%         aaa=mean(aa(ind1:ind2,:));
%         areaa = [areaa,sum(aaa(frame.cs_start:(frame.us_start-1)))];%trapz(aaa(frame.cs_start:(frame.cs_end-1)))
%     end
%     area.CS_acq_block(:,jj) = areaa;
%    area.CS_acq_1st_block(:,jj)=sum(aa(1+trial.hab(3):trial.acq_block_trial:trial.acq(3),frame.cs_start:(frame.us_start-1)),2);%trapz(aa(1+trial.hab(3):trial.acq_block_trial:trial.acq(3),frame.cs_start:(frame.cs_end-1)));
% end
%save([inputpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf');
% %test
% a=activities_new{11,1}(:,1);%b=activities_total_cut{1,1}(:,1);
% figure,plot(a,'b');hold on; %plot(b,'r')
