function [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf(activities_new)
global fs
global time
global frame
global frameb
global trial
global re_startpoint


set(0,'defaultfigurecolor','w');

%set para.
%[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);

%补齐前后两帧
activities={};
for ii = 1:size(activities_new,1)
    if ~isempty(activities_new{ii,1}(:))
        activities{ii,1}=[activities_new{ii,1}(1,:);activities_new{ii,1};activities_new{ii,1}(size(activities_new{ii,1},1),:);];
    else
        activities{ii,1}=[];
    end
    %activities{ii,1}(trial.acq(3)*frame.per_cycle+1:(trial.acq(3)+2)*frame.per_cycle,:)=[];
end

% % %插少帧点，没有少帧则注释掉
% for jj=trial.acq(2)+2*trial.acq_block_trial:trial.acq(2)+2*trial.acq_block_trial+3-1
%     for ii = 1:size(activities_new,1)
%         activities{ii,1}((jj-1)*frame.per_cycle+1:end+1,:)=[activities{ii,1}((jj-1)*frame.per_cycle+1,:);...
%             activities{ii,1}((jj-1)*frame.per_cycle+1:end,:)];
%     end
% end
% for ii = 1:size(activities_new,1)
%     activities_new{ii,1}=activities{ii,1}(2:end-1,:);
% end

%报告少帧
if size(activities{ii},1)<trial.total*frame.per_cycle
    error('frame number is wrong');
end
%detrend
% for ii = 1:size(activities,1)
%     for jj=1:size(activities{ii},2)
%     activities_detrend{ii}(:,jj)=detrend(activities{ii}(:,jj));
%     end
% end
% %smooth and 补齐到3000，
% for ii = 1:size(activities,1)
%     if ~isempty(activities_new{ii,1}(:))
%         if size(activities{ii},1)<trial.total*frame.per_cycle
%             %activities{ii,1}=[activities{ii,1}(1,:);activities{ii,1};];
%             a=ones(frame.per_cycle*trial.total,size(activities{ii},2)).*repmat(activities{ii}(end,:),frame.per_cycle*trial.total,1);
%             a(1:size(activities{ii},1),:)=activities{ii};
%             activities_total{ii,1}=a;
%         else
%             activities_total{ii,1}=activities{ii};
%         end
%     else
%        activities_total{ii,1}=[];
%     end
% end
% %测试smooth和去离群值，效果都不好
% aa=reshape(activities{ii}(:,jj),frame.per_cycle,trial.total)';
% test_smooth(aa',fs,frame,trial)
% test_isoutlier(aa',fs,frame,trial)
        
%dfdf,用每个trial前段时间的基线
activities_preCS_baseline={};activities_preCS_dfdf={};
baseline_win=1:frame.cs_start-1;%ceil(4.8/fs.ca+1)
for ii=1:max(size(activities))
    a=activities{ii};
    for jj=1:size(a,2)  
    aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
    %base=getbaseline_cut3sd(aa,1:frame.cs_start-1);
    base=getbaseline_cut3sd_strict_for_new_paradigm(aa,baseline_win,trial);
    activities_preCS_baseline{ii,1}(:,jj)=base.rep_m;
    activities_preCS_baseline_cutpersent{ii,1}(:,jj)=base.baseline_cut_percent_results;
    activities_preCS_dfdf{ii,1}(:,jj)=(a(:,jj)-base.rep_m)./base.rep_m;
    %activities_preCS_percent{ii,1}{jj,1}=base.baseline_cut_percent;
    
%     figure,plot(a(:,jj));hold on;plot(reshape(base.output_cut3sd',1,[])','r');hold on;
%     plot(base.rep_m);
    end
end

%去下降点
activities_preCS_dfdf_aftcorrect={};
for ii=1:max(size(activities))
    a=activities_preCS_dfdf{ii};
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
        activities_preCS_dfdf_aftcorrect{ii,1}(:,jj)=aa_cut;
%         figure,plot(aa,'b');hold on;plot(aa_cut,'r');hold on;
%         scatter(ind_cut,aa(ind_cut),[],'k');
    end
end
% figure,
% plot_rawtrace_trials( activities{1,1}(:,2),event,fs,frame,trial,startpoint,1)
% figure,
% plot_rawtrace_trials( activities_preCS_dfdf{1,1}(:,2),event,fs,frame,trial,startpoint,1)
% figure,
% plot_rawtrace_trials(activities_preCS_dfdf_aftcorrect{1,1}(:,2),event,fs,frame,trial,startpoint,1)

% figure,plot_rawtrace_trials(activities_preCS_dfdf_aftcorrect{1,1}(:,1),[],fs,frame,trial,[],1);title('aft process')
% figure,plot_rawtrace_trials(activities{1,1}(:,1),[],fs,frame,trial,[],1);title('raw trace')

%计算onset
activities_event_preCS={};
for ii=1:max(size(activities))
    a=activities_preCS_dfdf_aftcorrect{ii};
    for jj=1:size(a,2) 
        event=getevent(a(:,jj),trial,frame);
        activities_event_preCS{ii,1}{jj,1}=event;
    end
end


%save([inputpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf');

