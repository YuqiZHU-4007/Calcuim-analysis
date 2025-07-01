clc;clear all;
set(0,'defaultfigurecolor','w');
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 34;

trial.total = 60;
trial.hab = 15;
trial.acq = 30;
trial.test = 15;
fs_ca=0.4;
fs_behavior=60;

[inputname,inputpath]=uigetfile('G:\data\.mat','activities');
activities_new=importdata([inputpath,inputname]);
%save('I:\fish2_copy\analysis_manually\activities_new.mat','activities_new');
activities={};
for ii = 1:size(activities_new,1)
   activities{ii,1}=[activities_new{ii,1}(1,:);activities_new{ii,1};activities_new{ii,1}(size(activities_new{ii,1},1),:);];
    %activities{ii,1}=activities_new{ii,1};
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
    a=activities_total{ii,1};
    for jj=1:size(a,2)  
    aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
    base=getbaseline_cut3sd_strict(aa,1:frame.cs_start-1,trial);
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
        base=getbaseline_cut3sd_strict(aa,1:frame.cs_start-1,trial);
        m=base.rep_m;sd=base.rep_sd;
        % m=mean(aa(:,1:frame.cs_start-1),2);m=kron(m,ones(frame.per_cycle,1));
        % sd=std(aa(:,1:frame.cs_start-1),[],2);sd=kron(sd,ones(frame.per_cycle,1));
        baseline=m+2*sd;maxthres=m+3*sd;
        rasing_width=[1 100];
        decling_width=[4 100];%decling_width+rasing_width-1<8
        event=getevent_Trianglefitting_strict(a(:,jj),baseline,maxthres,rasing_width,decling_width,frame,trial,base);
        activities_event_preCS{ii,1}{jj,1}=event;
    end
end

save([inputpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf');
% %test
% a=activities_new{11,1}(:,1);%b=activities_total_cut{1,1}(:,1);
% figure,plot(a,'b');hold on; %plot(b,'r')

%画图
outpath='H:\20180623\fish2\summary_rec155212227\rec155212227_calcium_raw\activities\figure\raw activity_dfdf_baseline_preCs_aftcorrect_20180929';
st=1;
input=activities_preCS_dfdf_aftcorrect;
%画热图
checkpath([outpath '\heatmap']);
checkpath([outpath '\raw trace']);
for ii=1:size(activities,1)
    for jj=1:size(activities{ii},2)
        %input=activities{ii}(:,jj);
        lab='dfdf';b=0;%Intensity
        
        a=input{ii}(:,jj);%a(find(a==0),:)=[];
        aa=reshape(a,frame.per_cycle,trial.total)';
        x=([1:frame.per_cycle]-frame.cs_start)*0.4;
        y=1:60;
        h1=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        set(axes,'position',[0.03,0.05,0.95,0.9]);
        %set(gca,'YTick',0:1:trial.total);
        imagesc(x,y,aa);hold on %[min(min_activity_dfdf) max(max_activity_dfdf)];[-1 1]
        title([num2str(ii) 'layer-' num2str(jj)]);xlabel('time(s)');ylabel('trial num.');colorbar;
        line(([frame.cs_start frame.cs_start]-frame.cs_start)*0.4,[1 trial.total],'color',[0 0 0],'linestyle','--');hold on;
        line(([frame.us_start frame.us_start]-frame.cs_start)*0.4,[trial.hab+1-0.5 trial.hab+trial.acq+0.5],'color','r','linestyle','-');hold on;
        line(([frame.cs_end frame.cs_end]-frame.cs_start)*0.4,[1 trial.total],'color','k','linestyle','--');hold on;
        
        %behavior
        onset=[];onset=re_startpoint_sd;
        line(([onset(:,2) onset(:,2)]'-10*fs_behavior)/fs_behavior,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[1 1 0]);hold on;
        event=activities_event_preCS{ii,1}{jj,1};onset=[];
        onset(:,1)=ceil(event.ind/frame.per_cycle);
        onset(:,2)=mod(event.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[1 1 1],'linewidth',1.2);hold on;
        onset=[];onset(:,1)=ceil(event.end_ind/frame.per_cycle);
        onset(:,2)=mod(event.end_ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[0.5 0.5 0.5],'linewidth',1.2);hold on;

        %raw dfdf
        h2=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        subplot(3,1,1);
       plot([1:frame.per_cycle*trial.total]*0.4,a);hold on;xlim([1 trial.hab*frame.per_cycle]*0.4);hold on;

        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        %plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'r','LineStyle','--');hold on
        
        scatter(event.ind*0.4,a(event.ind),9,'k','filled');hold on;
        title([num2str(ii) '-layer' num2str(jj) ' Hab']);xlabel('time(s)');ylabel(lab);
        
        subplot(3,1,2);
        plot([1:frame.per_cycle*trial.total]*0.4,a,'r');hold on; xlim([trial.hab*frame.per_cycle+1 (trial.hab+trial.acq)*frame.per_cycle]*0.4);
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)'*0.4,[min(a(:)) max(a(:))],'k');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'b','LineStyle','--');hold on
               
        scatter(event.ind*0.4,a(event.ind),9,'k','filled');hold on;
        title(['Acq']);xlabel('time(s)');ylabel(lab);
        
        subplot(3,1,3);
        plot([1:frame.per_cycle*trial.total]*0.4,a);hold on;xlim([(trial.hab+trial.acq)*frame.per_cycle+1 (trial.total)*frame.per_cycle]*0.4);
        
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        %plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'r','LineStyle','--');hold on
              
        scatter(event.ind*0.4,a(event.ind),9,'k','filled');hold on;
        title(['Test']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
        
        switch st
            case 1
                saveas(h1,[outpath '\heatmap\' num2str(ii) 'layer-' num2str(jj) '.tif']);
                saveas(h2,[outpath '\raw trace\' num2str(ii) 'layer-' num2str(jj) '.tif']);
               % disp('1')
            case 2
                saveas(h1,[checkpath([outpath '\heatmap\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.tif']);
                savefig(h1,[checkpath([outpath '\heatmap\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)]);
                saveas(h2,[checkpath([outpath '\raw trace\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.tif']);
                savefig(h2,[checkpath([outpath '\raw trace\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)]);
        end
        close(h1);
        close(h2);
    end
end

%算均值
checkpath([outpath '\mean acorss trials']);
for ii=1:size(activities,1)
    for jj=1:size(activities{ii},2)
        %input=activities{ii}(:,jj);
        
        a=reshape(input{ii}(:,jj),frame.per_cycle,trial.total)';
        mean_hab_CSCS=mean(a(1:trial.hab,frame.cs_start:frame.cs_end-1),2);
        mean_acq_CSUS=mean(a(trial.hab+1:trial.acq+trial.hab,frame.cs_start:frame.us_start-1),2);
        mean_acq_USCS=mean(a(trial.hab+1:trial.acq+trial.hab,frame.us_start:frame.cs_end-1),2);
        mean_test_CSCS=mean(a(trial.acq+trial.hab:trial.total,frame.cs_start:frame.cs_end-1),2);
        h3=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        
        subplot(4,1,1);plot(mean_hab_CSCS);hold on;
        scatter(1:length(mean_hab_CSCS),mean_hab_CSCS);ylim([min(a(:)) max(a(:))]);xlim([1 trial.hab]);
        title([num2str(ii) 'layer-' num2str(jj) '---mean of dfdf between CS-on and CS-off(Hab.)']);xlabel('trials');
        
        subplot(4,1,2);plot(mean_acq_CSUS);hold on;
        scatter(1:length(mean_acq_CSUS),mean_acq_CSUS);ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_acq_CSUS)]);
        title('between CS-on and US-on(Acq.)');xlabel('trials');
        
        subplot(4,1,3);plot(mean_acq_USCS);hold on;
        scatter(1:length(mean_acq_USCS),mean_acq_USCS);ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_acq_USCS)]);
        title('between US-on and CS-off(Acq.)');xlabel('trials');
        
        subplot(4,1,4);plot(mean_test_CSCS,'b');hold on;plot(mean_hab_CSCS,'y');hold on;%legend('Test','Hab.')
        scatter(1:length(mean_test_CSCS),mean_test_CSCS,'r');hold on;
        scatter(1:length(mean_hab_CSCS),mean_hab_CSCS,'g');
        ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_test_CSCS)]);
        title('between CS-on and CS-off(Test)');xlabel('trials');
        
        switch st
            case 1
                saveas(h3,[outpath '\mean acorss trials\' num2str(ii) 'layer-' num2str(jj) '.tif']);
            case 2
                saveas(h3,[checkpath([outpath '\mean acorss trials\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.tif']);
                savefig(h3,[checkpath([outpath '\mean acorss trials\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)])
        end
        close(h3);
    end
end

%单独画acq和test（hab）
checkpath([outpath '\dfdf_acq']);
checkpath([outpath '\dfdf_test']);
for ii=1:size(activities,1)
    for jj=1:size(activities{ii},2)
        aa=input{ii}(:,jj); %activities_smooth_dfdf_process{ii}(:,jj);
        a=aa;a(find(a==0))=[];
        
        %figure,plot(activities{ii}(:,jj));hold on;plot(activities_baseline{ii}(:,jj));
        %figure,plot(a,'b');hold on;plot(activities_dfdf_process{ii}(:,jj),'r');legend()
        %h4=figure;
        for zz=1:5
         x=[0:6*frame.per_cycle-1]*0.4/60+8*(zz-1);
        h4=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        plot(x,aa([trial.hab*frame.per_cycle+1:(trial.hab+6)*frame.per_cycle]+6*frame.per_cycle*(zz-1),:),'r');hold on
        line([repmat([frame.cs_start:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.cs_end:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.us_start:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','g','LineStyle','--');hold on;
        line([repmat([0:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','b','LineStyle','-');hold on;
        title([num2str(ii) 'layer-' num2str(jj) '  Acq.-' num2str(zz)]);ylim([min(a(:)) max(a(:))]);%
        xlabel('time(min)');ylabel('dfdf');
        event=activities_event_preCS{ii,1}{jj,1};onset=[];
        onset=event.ind(find(event.ind>=(trial.hab*frame.per_cycle+1)+6*frame.per_cycle*(zz-1) & event.ind<=(trial.hab+6)*frame.per_cycle+6*frame.per_cycle*(zz-1)));
        scatter((onset-6*frame.per_cycle*(zz-1)-trial.hab*frame.per_cycle-1)*0.4/60+8*(zz-1),aa(onset),8,'k','filled');
        
        switch st
            case 1
                saveas(h4,[outpath '\dfdf_acq\' num2str(ii) 'layer-' num2str(jj) '-' num2str(zz) '.tif']);
            case 2
                saveas(h4,[checkpath([outpath '\dfdf_acq\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '-' num2str(zz) '.tif']);
                savefig(h4,[checkpath([outpath '\dfdf_acq\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '-' num2str(zz)]);
        end
        close(h4);
        end
        
        
        
        h5=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        x=[0:trial.hab*frame.per_cycle]*0.4/60;
        plot(aa([1:trial.hab*frame.per_cycle],:));hold on;
        plot(aa([(trial.hab+trial.acq)*frame.per_cycle+1:trial.total*frame.per_cycle],:),'r');hold on
        %legend('Hab.','Test');
        ylim([min(a(:)) max(a(:))]);
        xlabel('frames');ylabel('dfdf');title([num2str(ii) 'layer-' num2str(jj) '  Test(r) and Hab.(b)' num2str(zz)]);
        line([repmat([frame.cs_start:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on;
        line([repmat([frame.cs_end:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color',[1 0.5 0],'LineStyle','--');hold on;
        %line([repmat([frame.us_start:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','g','LineStyle','--');hold on;
        line([repmat([0:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','g','LineStyle','-');
        
        switch st
            case 1
                saveas(h5,[outpath '\dfdf_test\' num2str(ii) 'layer-' num2str(jj) '.png']);
            case 2
                saveas(h5,[checkpath([outpath '\dfdf_test\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.tif']);
                savefig(h5,[checkpath([outpath '\dfdf_test\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)]);
        end
        %saveas(h5,[outpath '\dfdf_test\'  num2str(ii) 'layer-' num2str(jj) '.png']);
        close(h5);
    end
end