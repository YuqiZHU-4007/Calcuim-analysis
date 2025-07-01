set(0,'defaultfigurecolor','w');
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 36;

%6:自发
%15：hab
%6*4：acq
%6:test
%6:spon

%1:trial num;2:start trial;3:end trial
trial.spon_bef=[6*3 1 6*3];
trial.test = [6 trial.spon_bef(3)+1 trial.spon_bef(3)+6];
trial.total =trial.spon_bef(1)+trial.test(1);

fs_ca=0.4;
fs_behavior=60;

[inputname,inputpath]=uigetfile('G:\data\.mat','activities');
activities_new=importdata([inputpath,inputname]);
%save('I:\fish2_copy\analysis_manually\activities_new.mat','activities_new');
activities={};
for ii = 1:size(activities_new,1)
    activities{ii,1}=[activities_new{ii,1}(1,:);activities_new{ii,1};activities_new{ii,1}(size(activities_new{ii,1},1),:);];
end

%detrend
% for ii = 1:size(activities,1)
%     for jj=1:size(activities{ii},2)
%     activities_detrend{ii}(:,jj)=detrend(activities{ii}(:,jj));
%     end
% end
%smooth and 补齐到3000
for ii = 1:size(activities,1)
    if size(activities{ii},1)<trial.total*frame.per_cycle
        %activities{ii,1}=[activities{ii,1}(1,:);activities{ii,1};];
        a=ones(frame.per_cycle*trial.total,size(activities{ii},2)).*repmat(activities{ii}(end,:),frame.per_cycle*trial.total,1);
        a(1:size(activities{ii},1),:)=activities{ii};
        activities_total{ii,1}=a;
    else
        activities_total{ii,1}=activities{ii};   
    end
end


%dfdf,用每个trial前段时间的基线
for ii=1:size(activities,1)
    a=activities_total{ii};
    for jj=1:size(a,2)  
    aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
    base=getbaseline_cut3sd(aa,1:frame.cs_start-1);
    activities_preCS_baseline{ii,1}(:,jj)=base.rep_m;
    activities_preCS_dfdf{ii,1}(:,jj)=(a(:,jj)-base.rep_m)./base.rep_m;
    %activities_preCS_percent{ii,1}{jj,1}=base.baseline_cut_percent;
%     figure,plot(a(:,jj));hold on;plot(reshape(base.output_cut3sd',1,[])','r');hold on;
%     plot(base.rep_m);
    end
end

%去下降点
for ii=1:size(activities,1)
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

%计算onset
for ii=1:size(activities,1)
    a=activities_preCS_dfdf_aftcorrect{ii};
    for jj=1:size(a,2) 
        aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
        base=getbaseline_cut3sd(aa,1:frame.cs_start-1);
        m=base.rep_m;sd=base.rep_sd;
        % m=mean(aa(:,1:frame.cs_start-1),2);m=kron(m,ones(frame.per_cycle,1));
        % sd=std(aa(:,1:frame.cs_start-1),[],2);sd=kron(sd,ones(frame.per_cycle,1));
        baseline=m+2*sd;maxthres=m+3*sd;
        rasing_width=[1 100];
        decling_width=[4 100];%decling_width+rasing_width-1<8
        event=getevent_Trianglefitting_strict_for_new_paradigm_day2(a(:,jj),baseline,maxthres,rasing_width,decling_width,frame,base);
        activities_event_preCS{ii,1}{jj,1}=event;
    end
end

save([inputpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf');
% %test
% a=activities_new{11,1}(:,1);%b=activities_total_cut{1,1}(:,1);
% figure,plot(a,'b');hold on; %plot(b,'r')
