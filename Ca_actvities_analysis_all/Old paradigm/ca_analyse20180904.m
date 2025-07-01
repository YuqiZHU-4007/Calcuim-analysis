set(0,'defaultfigurecolor','w')
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 34;

trial.total = 60;
trial.hab = 15;
trial.acq = 30;
trial.test = 15;
%save('I:\fish2_copy\analysis_manually\activities_new.mat','activities_new');
activities=activities_new;
for ii = 1:size(activities,1)
    activities{ii,1}=[activities{ii,1}(1,:);activities{ii,1};activities{ii,1}(size(activities{ii,1},1),:);];
end

%detrend
% for ii = 1:size(activities,1)
%     for jj=1:size(activities{ii},2)
%     activities_detrend{ii}(:,jj)=detrend(activities{ii}(:,jj));
%     end
% end
%smooth and 补齐到3000
if size(activities{1,1},1)<3000
    for ii = 1:size(activities,1)
        %activities{ii,1}=[activities{ii,1}(1,:);activities{ii,1};];
        a=zeros(frame.per_cycle*trial.total,size(activities{ii},2));
        a(1:size(activities{ii},1),:)=activities{ii};
        activities_total{ii,1}=a;
        activities_smooth{ii,1}=smoothdata(activities{ii},1,'gaussian',10,'omitnan');
    end
else
    activities_total=activities;
end
%dfdf(不同层基线波动不同，基线平的用均值或加大winsize)
activities_dfdf = cell(16,1);
activities_baseline=cell(16,1);
for ii = 1:size(activities,1)%16
    [dfdf, baselines] = compute_dfdf(activities{ii}(:,:)',45,0.1);
    activities_dfdf{ii,1} = dfdf';
    activities_baseline{ii,1}=baselines';
    [dfdf, baselines] = compute_dfdf(activities_smooth{ii}(:,:)',45,0.1);
    activities_smooth_dfdf{ii,1} = dfdf';
    activities_smooth_baseline{ii,1}=baselines';
end

%补齐到3000
for ii = 1:size(activities,1)
    x=activities_smooth_dfdf{ii};
    a=zeros(frame.per_cycle*trial.total,size(activities{ii},2));
    a(1:size(x,1),:)=x;
    activities_smooth_dfdf{ii,1}=a;
    
    x=activities_dfdf{ii};
    a=zeros(frame.per_cycle*trial.total,size(activities{ii},2));
    a(1:size(x,1),:)=x;
    activities_dfdf{ii,1}=a;
end

%process,去基线突变
activities_smooth_dfdf_process = activities_smooth_dfdf;
for ii=1:size(activities,1)
    a=activities{ii};
    for jj=1:size(a,2)
        %figure,plot(activities_dfdf{ii}(:,jj));
        %figure,plot(a(:,jj));
        ind=51:frame.per_cycle:size(a,1);
        %判断在每个trial的起始点是不是有基线突变
        %         erase_point=ind(find(abs(a(ind,jj)-a(max(ind-1,1),jj))./activities_baseline{ii}(ind,jj)>...
        %             (mean(activities_dfdf{ii}(:,jj))+2*std(activities_dfdf{ii}(:,jj))))); %dfdf>0.05
        erase_point2=ind(find(abs(a(ind,jj)-a(ind-1,jj))>...
            abs(mean([a(ind-4,jj) a(ind-3,jj) a(ind-2,jj) a(ind-1,jj)]-[a(ind-5,jj) a(ind-4,jj) a(ind-3,jj) a(ind-2,jj)],2)+...
            2*std([a(ind-4,jj) a(ind-3,jj) a(ind-2,jj) a(ind-1,jj)]-[a(ind-5,jj) a(ind-4,jj) a(ind-3,jj) a(ind-2,jj)],[],2))...
            ));% 和前一个点的差大于m（前5个点之间的差异）+2*sd
        if ~isempty(erase_point2)
            for zz=1:length(erase_point2) %基线突变的点的后12个点的dfdf设为均值
                a_gaussian= smoothdata(a(erase_point2(zz):min(erase_point2(zz)+50-1,length(a)),jj),'gaussian',6,'omitnan');
                activities_smooth_dfdf_process{ii}(erase_point2(zz):min(erase_point2(zz)+12-1,length(a)),jj)=(a(erase_point2(zz):min(erase_point2(zz)+12-1,length(a)),jj)-a_gaussian(1:min(1+12-1,length(a_gaussian))))./a_gaussian(1:min(1+12-1,length(a_gaussian)));
                %             figure,plot(a_gaussian,'r');hold on;plot(a(erase_point2(zz):erase_point2(zz)+50-1,jj));legend('gaussian','raw');
                %             figure,plot((a(erase_point2(zz):erase_point2(zz)+50-1,jj)-a_gaussian)./a_gaussian);hold on;
                %             plot(activities_dfdf{ii}(erase_point2(zz):erase_point2(zz)+50-1,jj),'r');legend('gaussian_dfdf','raw_dfdf');
                %mean(activities_dfdf{ii}(erase_point2:erase_point2+frame.per_cycle-1,jj));
            end
        end
        %         figure,plot(activities_dfdf_process{ii}(:,jj),'r');hold on
        %         plot(activities_dfdf{ii}(:,jj),'b');hold on
    end
end

activities_dfdf_process = activities_dfdf;
for ii=1:size(activities,1)
    a=activities{ii};
    for jj=1:size(a,2)
        %figure,plot(activities_dfdf{ii}(:,jj));
        %figure,plot(a(:,jj));
        ind=51:frame.per_cycle:size(a,1);
        %判断在每个trial的起始点是不是有基线突变
        %         erase_point=ind(find(abs(a(ind,jj)-a(max(ind-1,1),jj))./activities_baseline{ii}(ind,jj)>...
        %             (mean(activities_dfdf{ii}(:,jj))+2*std(activities_dfdf{ii}(:,jj))))); %dfdf>0.05
        erase_point2=ind(find(abs(a(ind,jj)-a(ind-1,jj))>...
            abs(mean([a(ind-4,jj) a(ind-3,jj) a(ind-2,jj) a(ind-1,jj)]-[a(ind-5,jj) a(ind-4,jj) a(ind-3,jj) a(ind-2,jj)],2)+...
            2*std([a(ind-4,jj) a(ind-3,jj) a(ind-2,jj) a(ind-1,jj)]-[a(ind-5,jj) a(ind-4,jj) a(ind-3,jj) a(ind-2,jj)],[],2))...
            ));% 和前一个点的差大于m（前5个点之间的差异）+2*sd
        if ~isempty(erase_point2)
            for zz=1:length(erase_point2) %基线突变的点的后12个点的dfdf设为均值
                a_gaussian= smoothdata(a(erase_point2(zz):min(erase_point2(zz)+50-1,length(a)),jj),'gaussian',6,'omitnan');
                activities_dfdf_process{ii}(erase_point2(zz):min(erase_point2(zz)+12-1,length(a)),jj)=(a(erase_point2(zz):min(erase_point2(zz)+12-1,length(a)),jj)-a_gaussian(1:min(1+12-1,length(a_gaussian))))./a_gaussian(1:min(1+12-1,length(a_gaussian)));
                %             figure,plot(a_gaussian,'r');hold on;plot(a(erase_point2(zz):erase_point2(zz)+50-1,jj));legend('gaussian','raw');
                %             figure,plot((a(erase_point2(zz):erase_point2(zz)+50-1,jj)-a_gaussian)./a_gaussian);hold on;
                %             plot(activities_dfdf{ii}(erase_point2(zz):erase_point2(zz)+50-1,jj),'r');legend('gaussian_dfdf','raw_dfdf');
                %mean(activities_dfdf{ii}(erase_point2:erase_point2+frame.per_cycle-1,jj));
            end
        end
        %         figure,plot(activities_dfdf_process{ii}(:,jj),'r');hold on
        %         plot(activities_dfdf{ii}(:,jj),'b');hold on
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
    figure,plot(a(:,jj));hold on;plot(reshape(base.output_cut3sd',1,[])','r');hold on;
    plot(base.rep_m);
    
%     %[~,b]=compute_dfdf(aa,11,0.5);bb=b(:,1:frame.cs_start-1);    
%     b=aa(:,1:frame.cs_start-1);bb=sort(b,2);bb=bb(:,1:fix(frame.cs_start*0.90));
%     
%     %figure,plot(aa(2,:));hold on;plot(bb(2,:));hold on;
%     %line([1 50],[mean(aa(2,:)) mean(aa(2,:))]);hold on; line([1 50],[mean(aa(2,1:25)) mean(aa(2,1:25))],'color','r');
%     
%     b(:,frame.cs_start:frame.per_cycle )= repmat(mean(bb,2),1,length(frame.cs_start:frame.per_cycle))+3*repmat(std(bb,[],2),1,length(frame.cs_start:frame.per_cycle));
%     m=repmat(mean(bb,2),1,size(aa,2));sd=repmat(std(bb,[],2),1,size(aa,2));mm=reshape(m',size(aa,1)*size(aa,2),1);
%     ind=zeros(size(aa));ind=abs(b-m)-3*sd>0; ind(:,frame.cs_start:frame.per_cycle)=0;ind=reshape(ind',size(aa,1)*size(aa,2),1);%bb=reshape(b',size(b,1)*size(b,2),1);
%     if ~isempty(find(ind==1))
%     a(find(ind==1),jj)=mm(find(ind==1));%interp1(1:length(a(:,jj)),a(:,jj),find(ind==1),'linear');
%     %figure,plot(activities_total{ii}(:,jj));hold on;plot(a(:,jj),'r');
%     end
%     aa=reshape(a(:,jj),frame.per_cycle,trial.total)';
%     activities_preCS_baseline{ii,1}(:,jj)=reshape(repmat(mean(aa(:,1:frame.cs_start-1),2),1,size(aa,2))',frame.per_cycle*trial.total,1);
%     activities_preCS_dfdf{ii,1}(:,jj)=(activities_total{ii}(:,jj)-activities_preCS_baseline{ii}(:,jj))./activities_preCS_baseline{ii}(:,jj);
%     %figure,plot(activities_total{ii}(:,jj));hold on;plot(activities_preCS_baseline{ii}(:,jj),'r');
    end
end


%test 用哪个
ii=11;jj=1;
mean(activities_baseline{ii}(:,jj))
figure,plot(activities{ii}(:,jj),'b');hold on;plot(activities_baseline{ii}(:,jj),'y');
plot(activities_smooth{ii}(:,jj),'r');hold on;plot(activities_smooth_baseline{ii}(:,jj),'g');
figure,plot(activities_dfdf{ii}(:,jj),'b');hold on;plot(activities_smooth_dfdf{ii}(:,jj),'r');
   
%smooth dfdf
% activities_dfdf_process_sgolay=cell(size(activities_dfdf));
% for ii=1:size(activities,1)
%    % activities_dfdf_process_guassian{ii}= smoothdata(activities_dfdf_process{ii},1,'gaussian',6,'omitnan'); 
%     a=zeros(frame.per_cycle*trial.total,size(activities_dfdf{ii},2));
%     a(1:size(activities_dfdf{ii},1),:)=smoothdata(activities_dfdf{ii},1,'rloess',6,'omitnan');
%     activities_dfdf_process_sgolay{ii}=a;
%      min_activity_dfdf(ii)=min(activities_dfdf_process_sgolay{ii}(:));
%      max_activity_dfdf(ii)=max(activities_dfdf_process_sgolay{ii}(:));
% %     figure, plot(activities_dfdf_process{ii}(:,6),'b');hold on
% %     %plot(activities_dfdf_process_guassian{ii}(:,6),'r');hold on
% %     plot(activities_dfdf_process_sgolay{ii}(:,6),'y');hold on
% end

%去下降点
for ii=1:size(activities,1)
    a=activities_preCS_dfdf{ii};
    for jj=1:size(a,2) 
        aa=a(:,jj);
        m=mean(aa);sd=std(aa);
        
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
        
%         aa_cut=myinterp(1:length(aa),aa',ind_cut);
%         figure,plot(aa,'b');hold on;plot(aa_cut,'r');hold on;
%         scatter(ind_cut,aa(ind_cut),[],'k');
    end
end

%计算onset


%画图
outpath='H:\20180623_fish1\activity_analyze\figure\raw activity_dfdf_baseline_preCs_aftcorrect';
st=1;
%input=activities_smooth_dfdf_process;
%画热图
checkpath([outpath '\hotmap']);
checkpath([outpath '\raw trace']);
for ii=1:size(activities,1)
    for jj=1:size(activities{ii},2)
        input=activities{ii}(:,jj);
        lab='dfdf';b=0;%Intensity
        
        a=input;a(find(a==0),:)=[];
        aa=reshape(input,frame.per_cycle,trial.total)';
        x=([1:frame.per_cycle]-frame.cs_start)*0.4;
        y=1:60;
        h1=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        set(axes,'position',[0.03,0.05,0.95,0.9]);
        imagesc(x,y,aa);hold on %[min(min_activity_dfdf) max(max_activity_dfdf)];[-1 1]
        title([num2str(ii) 'layer-' num2str(jj)]);xlabel('time(s)');ylabel('trials');colorbar;
        line(([frame.cs_start frame.cs_start]-frame.cs_start)*0.4,[1 trial.total],'color','k','linestyle','--');hold on;
        line(([frame.us_start frame.us_start]-frame.cs_start)*0.4,[trial.hab+1 trial.hab+trial.acq],'color','r','linestyle','-');hold on;
        line(([frame.cs_end frame.cs_end]-frame.cs_start)*0.4,[1 trial.total],'color','k','linestyle','--');hold on;
        event=

        %raw dfdf
        h2=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        set(axes,'position',[0.03,0.05,0.95,0.9]);
        subplot(3,1,1);plot([1:frame.per_cycle*trial.total]*0.4,input);hold on;xlim([1 trial.hab*frame.per_cycle]*0.4);
        title([num2str(ii) '-layer' num2str(jj) ' Hab']);xlabel('time(s)');ylabel(lab);
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        %plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'r','LineStyle','--');hold on
        
        subplot(3,1,2);
        plot([1:frame.per_cycle*trial.total]*0.4,input,'r');hold on;xlim([trial.hab*frame.per_cycle+1 (trial.hab+trial.acq)*frame.per_cycle]*0.4);
        title(['Acq']);xlabel('time(s)');ylabel(lab);
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)'*0.4,[min(a(:)) max(a(:))],'k');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'b','LineStyle','--');hold on
        
        subplot(3,1,3);plot([1:frame.per_cycle*trial.total]*0.4,input);hold on;xlim([(trial.hab+trial.acq)*frame.per_cycle+1 (trial.total)*frame.per_cycle]*0.4);
        title(['Test']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'g');hold on
        %plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
        plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(a(:)) max(a(:))],'y','LineStyle','--');hold on
        plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(a(:)) max(a(:))],'r','LineStyle','--');hold on
        
        switch st
            case 1
                saveas(h1,[outpath '\raw trace\' num2str(ii) 'layer-' num2str(jj) '.png']);
                saveas(h2,[outpath '\raw trace\' num2str(ii) 'layer-' num2str(jj) '.png']);
               % disp('1')
            case 2
                saveas(h1,[checkpath([outpath '\hotmap\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.png']);
                savefig(h1,[checkpath([outpath '\hotmap\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)]);
                saveas(h2,[checkpath([outpath '\raw trace\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.png']);
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
        input=activities{ii}(:,jj);
        
        a=reshape(input,frame.per_cycle,trial.total)';
        mean_hab_CSCS=mean(a(1:trial.hab,frame.cs_start:frame.cs_end-1),2);
        mean_acq_CSUS=mean(a(trial.hab+1:trial.acq+trial.hab,frame.cs_start:frame.us_start-1),2);
        mean_acq_USCS=mean(a(trial.hab+1:trial.acq+trial.hab,frame.us_start:frame.cs_end-1),2);
        mean_test_CSCS=mean(a(trial.acq+trial.hab:trial.total,frame.cs_start:frame.cs_end-1),2);
        h3=figure;title('x') 
        subplot(4,1,1);plot(mean_hab_CSCS);hold on;
        scatter(1:length(mean_hab_CSCS),mean_hab_CSCS);ylim([min(a(:)) max(a(:))]);xlim([1 trial.hab]);
        title([num2str(ii) 'layer-' num2str(jj) '---mean of dfdf between CS-on and CS-off(Hab.)']);xlabel('trials');
        
        subplot(4,1,2);plot(mean_acq_CSUS);hold on;
        scatter(1:length(mean_acq_CSUS),mean_acq_CSUS);ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_acq_CSUS)]);
        title('between CS-on and US-on(Acq.)');xlabel('trials');
        
        subplot(4,1,3);plot(mean_acq_USCS);hold on;
        scatter(1:length(mean_acq_USCS),mean_acq_USCS);ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_acq_USCS)]);
        title('US-on and CS-off(Acq.)');xlabel('trials');
        
        subplot(4,1,4);plot(mean_test_CSCS,'b');hold on;plot(mean_hab_CSCS,'y');hold on;%legend('Test','Hab.')
        scatter(1:length(mean_test_CSCS),mean_test_CSCS,'r');hold on;
        scatter(1:length(mean_hab_CSCS),mean_hab_CSCS,'g');
        ylim([min(a(:)) max(a(:))]);xlim([1 length(mean_test_CSCS)]);
        title('CS-on and CS-off(Test)');xlabel('trials');
        
        switch st
            case 1
                saveas(h3,[outpath '\mean acorss trials\' num2str(ii) 'layer-' num2str(jj) '.png']);
            case 2
                saveas(h3,[checkpath([outpath '\mean acorss trials\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.png']);
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
        aa=activities{ii}(:,jj); %activities_smooth_dfdf_process{ii}(:,jj);
        a=aa;a(find(a==0))=[];
        
        %figure,plot(activities{ii}(:,jj));hold on;plot(activities_baseline{ii}(:,jj));
        %figure,plot(a,'b');hold on;plot(activities_dfdf_process{ii}(:,jj),'r');legend()
        %h4=figure;
        for zz=1:5
         x=[0:6*frame.per_cycle-1]*0.4/60+8*(zz-1);
        h4=figure;plot(x,aa([trial.hab*frame.per_cycle+1:(trial.hab+6)*frame.per_cycle]+6*frame.per_cycle*(zz-1),:),'r');hold on
        line([repmat([frame.cs_start:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.cs_end:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.us_start:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','g','LineStyle','--');hold on;
        line([repmat([0:frame.per_cycle:6*frame.per_cycle]*0.4/60,2,1)']+8*(zz-1),[min(a(:)) max(a(:))],'color','b','LineStyle','-');hold on;
        title([num2str(ii) 'layer-' num2str(jj) '  Acq.-' num2str(zz)]);ylim([min(a(:)) max(a(:))]);%
        xlabel('time(min)');ylabel('dfdf');
        
        switch st
            case 1
                saveas(h4,[outpath '\dfdf_acq\' num2str(ii) 'layer-' num2str(jj) '-' num2str(zz) '.png']);
            case 2
                saveas(h4,[checkpath([outpath '\dfdf_acq\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '-' num2str(zz) '.png']);
                savefig(h4,[checkpath([outpath '\dfdf_acq\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '-' num2str(zz)]);
        end
        close(h4);
        end
        
        
        
        h5=figure;
        x=[0:trial.hab*frame.per_cycle]*0.4/60;
        plot(aa([1:trial.hab*frame.per_cycle],:),'b');hold on;
        plot(aa([(trial.hab+trial.acq)*frame.per_cycle+1:trial.total*frame.per_cycle],:),'r');hold on
        %legend('Hab.','Test');
        ylim([min(a(:)) max(a(:))]);
        xlabel('frames');ylabel('dfdf');title([num2str(ii) 'layer-' num2str(jj) '  Test(r) and Hab.(b)' num2str(zz)]);
        line([repmat([frame.cs_start:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.cs_end:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','m','LineStyle','--');hold on;
        line([repmat([frame.us_start:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','g','LineStyle','--');hold on;
        line([repmat([0:frame.per_cycle:trial.hab*frame.per_cycle],2,1)'],[min(a(:)) max(a(:))],'color','b','LineStyle','-');
        
        switch st
            case 1
                saveas(h5,[outpath '\dfdf_test\' num2str(ii) 'layer-' num2str(jj) '.png']);
            case 2
                saveas(h5,[checkpath([outpath '\dfdf_test\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj) '.png']);
                savefig(h5,[checkpath([outpath '\dfdf_test\' num2str(ii) 'layer\' num2str(jj) '\']) num2str(ii) 'layer-' num2str(jj)]);
        end
        %saveas(h5,[outpath '\dfdf_test\'  num2str(ii) 'layer-' num2str(jj) '.png']);
        close(h5);
    end
end