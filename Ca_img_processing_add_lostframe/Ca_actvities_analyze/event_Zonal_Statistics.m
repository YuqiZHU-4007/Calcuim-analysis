function event_Zonal_Statistics(num,activities,frame,trial)

outpath='H:\20180623_fish1\分区\Olfactory';
%[filename,filepath]=uigetfile('H:\','.*');
num=xlsread(['H:\20180623_fish1\分区\分区_20180623_fish1.xlsx'],'Olfactory');

acti=[];a_event={};a_event_rep=[];
for ii=1:size(num,1)
    acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
    a_event{ii,1}=activities_event_preCS{num(ii,1),1}{num(ii,2),1}';
    a_event_rep=[a_event_rep a_event{ii,1}.ind];
end

%全画
checkpath([outpath '\hotmap']);
checkpath([outpath '\raw trace']);
for ii=1:size(num,1)
    %aa=acti(:,ii);
    lab='dfdf'; b=0;
    aa=reshape(acti(:,ii),frame.per_cycle,trial.total)';
    x=([1:frame.per_cycle]-frame.cs_start)*0.4;
    y=1:60;
    h1=figure;
    set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
    set(axes,'position',[0.03,0.05,0.95,0.9]);
    %set(gca,'YTick',0:1:trial.total);
    imagesc(x,y,aa);hold on %[min(min_activity_dfdf) max(max_activity_dfdf)];[-1 1]
    title([num2str(num(ii,1)) 'layer-' num2str(num(ii,2))]);xlabel('time(s)');ylabel('trial num.');colorbar;
    line(([frame.cs_start frame.cs_start]-frame.cs_start)*0.4,[1 trial.total],'color',[0 0 0],'linestyle','--');hold on;
    line(([frame.us_start frame.us_start]-frame.cs_start)*0.4,[trial.hab+1-0.5 trial.hab+trial.acq+0.5],'color','r','linestyle','-');hold on;
    line(([frame.cs_end frame.cs_end]-frame.cs_start)*0.4,[1 trial.total],'color','k','linestyle','--');hold on;
    
    %behavior
    onset=[];onset=re_startpoint_sd;
    line(([onset(:,2) onset(:,2)]'-10*fs_behavior-1)/fs_behavior,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[1 1 0],'linewidth',1.2);hold on;
    event=activities_event_preCS{num(ii,1),1}{num(ii,2),1};onset=[];
    onset(:,1)=ceil(event.ind/frame.per_cycle);
    onset(:,2)=mod(event.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
    line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[1 1 1],'linewidth',1.2);hold on;
    onset=[];onset(:,1)=ceil(event.end_ind/frame.per_cycle);
    onset(:,2)=mod(event.end_ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
    line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[0.5 0.5 0.5],'linewidth',1.2);hold on;
        
    saveas(h1,[outpath '\hotmap\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
    close(h1);
    
    h2=figure;
    aa=acti(:,ii);
    set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
    subplot(3,1,1);
    plot([1:frame.per_cycle*trial.total]*0.4,aa);hold on;xlim([1 trial.hab*frame.per_cycle]*0.4);hold on;ylim([min(aa(:)) max(aa(:))]);
    plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'g');hold on
    plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'r','LineStyle','--');hold on
    scatter(event.ind*0.4,aa(event.ind),9,'k','filled');hold on;
%     %behavior
%     onset=[];onset=startpoint_sd;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[max(aa(:))-0.1*(max(aa(:))-min(aa(:))) max(aa(:))],'color','k');hold on;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[min(aa(:)) min(aa(:))+0.1*(max(aa(:))-min(aa(:)))],'color','k');hold on;
    title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Hab']);xlabel('time(s)');ylabel(lab);
    
    subplot(3,1,2);
    plot([1:frame.per_cycle*trial.total]*0.4,aa,'r');hold on; xlim([trial.hab*frame.per_cycle+1 (trial.hab+trial.acq)*frame.per_cycle]*0.4);ylim([min(aa(:)) max(aa(:))]);
    plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'g');hold on
    plot(repmat([frame.per_cycle*trial.hab:frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab-1)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'k');hold on
    %plot(repmat([frame.per_cycle*(trial.hab+1):frame.per_cycle*6:frame.per_cycle*(trial.acq+trial.hab)],2,1)',[min(a(:)) max(a(:))],'b');hold on
    plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'b','LineStyle','--');hold on
    scatter(event.ind*0.4,aa(event.ind),9,'k','filled');hold on;
%     %behavior
%     onset=[];onset=startpoint_sd;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[max(aa(:))-0.1*(max(aa(:))-min(aa(:))) max(aa(:))],'color','k');hold on;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[min(aa(:)) min(aa(:))+0.1*(max(aa(:))-min(aa(:)))],'color','k');hold on;
    title(['Acq']);xlabel('time(s)');ylabel(lab);
    
    subplot(3,1,3);
    plot([1:frame.per_cycle*trial.total]*0.4,aa);hold on;xlim([(trial.hab+trial.acq)*frame.per_cycle+1 (trial.total)*frame.per_cycle]*0.4);ylim([min(aa(:)) max(aa(:))]);
    plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'g');hold on
    plot(repmat([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
    plot(repmat([(trial.hab*frame.per_cycle+frame.us_start):frame.per_cycle:frame.per_cycle*(trial.acq+trial.hab)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'r','LineStyle','--');hold on
    scatter(event.ind*0.4,aa(event.ind),9,'k','filled');hold on;
%     %behavior
%     onset=[];onset=startpoint_sd;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[max(aa(:))-0.1*(max(aa(:))-min(aa(:))) max(aa(:))],'color','k');hold on;
%     line(([onset;onset]'-601)/60+frame.cs_start*0.4,[min(aa(:)) min(aa(:))+0.1*(max(aa(:))-min(aa(:)))],'color','k');hold on;
    title(['Test']);xlabel('time(s)');ylabel(lab);%ylim([min(a(:))-b max(a(:))])
    
    saveas(h2,[outpath '\raw trace\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
    close(h2);
end


%总体
%hotmap
figure,
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
set(axes,'position',[0.054,0.08,0.90,0.9]);
imagesc(acti');%colorbar;
hold on;
xlabel('frame');
ylabel('cell number');
%behavior
onset=[];onset=startpoint_sd;
line([onset;onset]'/60/0.4,[1-0.5 size(num,1)+0.5],'color',[0 0 0]);hold on;
set(gca,'YTick',1:size(num,1));

%均值
figure,
y=mean(acti,2);
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
set(axes,'position',[0.054,0.08,0.90,0.6]);
patch([[1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.total/2,1)',...
    [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
    repmat([min(y) min(y) max(y) max(y)],trial.total,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off

line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
    [min(y) max(y)],'color','r','linewidth',1.05);hold on;
line([1:frame.per_cycle:frame.per_cycle*trial.total;1:frame.per_cycle:frame.per_cycle*trial.total],...
    [min(y) max(y)],'color',[0.5 0.5 0.5],'linestyle','--');hold on;

ylim([min(y) max(y)]);xlabel('frames');ylabel('dfdf');
plot(mean(acti,2),'k');hold on;
% line([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total;frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],...
%      [min(y) max(y)],'color','b');hold on;
% line([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total;frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],...
%      [min(y) max(y)],'color','b');hold on;

%onset
figure,
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
set(axes,'position',[0.054,0.08,0.90,0.6]);
xlabel('frames');ylabel('cell number');
patch([[1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
    repmat([0 0 size(num,1) size(num,1)],trial.total/2,1)',...
    [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
    repmat([0 0 size(num,1) size(num,1)],trial.total,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off

line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
    [0 size(num,1)],'color','r','linewidth',1);hold on;
% line([1:frame.per_cycle:frame.per_cycle*trial.total;1:frame.per_cycle:frame.per_cycle*trial.total],...
%     [0 size(num,1)],'color',[0.5 0.5 0.5],'linestyle','--');hold on;

for ii=1:size(num,1)
    line([a_event{ii,1}.ind;a_event{ii,1}.ind],[ii-0.3 ii+0.3],'color',[0 0 0],'linewidth',1);hold on;
    %scatter(a_event{ii,1}.ind,repmat(ii,length(a_event{ii,1}.ind),1),8,'k','filled');hold on;
end
ylim([0 size(num,1)])
set(gca,'YTick',0:1:size(num,1));

%onset统计图
figure,
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
set(axes,'position',[0.054,0.08,0.90,0.6]);
xlabel('frames');ylabel('number of events');
patch([[1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
    [1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
    repmat([0 0 size(num,1) size(num,1)],trial.total/2,1)',...
    [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
patch([[frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total]'...
    [frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total]']',...
    repmat([0 0 size(num,1) size(num,1)],trial.total,1)',...
    'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
    [size(num,1)-1 size(num,1)],'color','r','linewidth',1);hold on;
h=histogram(a_event_rep,frame.per_cycle*trial.total/4,'facecolor','k');hold on
% line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
%     [0 size(num,1)],'color','r','linewidth',1);hold on;
line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
    [0 2],'color','r','linewidth',1);hold on;
xlim([0 frame.per_cycle*trial.total]);ylim([1 size(num,1)])
%y=sum(a_event_rep,2);

%onset随时间变化
a_event_onset1=[];
a_event_onset2=[];
h=figure;
set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
set(axes,'position',[0.054,0.08,0.90,0.6]);
xlabel('trial number');ylabel('time(s)')
x=1:trial.total;
patch([x(1)  x(end) x(end) x(1)],([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*0.4,'g','edgecolor','g','facealpha',0.3,'edgealpha',0.2);hold on;
line([x(1) x(end)],([frame.us_start frame.us_start]-frame.cs_start)*0.4,'color','r');hold on;
xlim([1 trial.total]);ylim(([1 frame.per_cycle]-frame.cs_start)*0.4);
for ii=1:size(num,1)
    onset=[];
    onset(:,1)=ceil(a_event{ii,1}.ind/frame.per_cycle);
    onset(:,2)=mod(a_event{ii,1}.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
    scatter(onset(:,1),(onset(:,2)-frame.cs_start)*0.4,13);hold on;
    %plot(onset(:,1),(onset(:,2)-frame.cs_start)*0.4);hold on;
    a_event_onset1=[a_event_onset1;onset(:,1)];
    a_event_onset2=[a_event_onset2;onset(:,2)];
    
    ind=setdiff(x,intersect(onset(:,1)',x))';
    onset(size(onset,1)+1:size(onset,1)+length(ind),:)=[ind zeros(size(ind))];[~,ind]=sort(onset(:,1));onset=onset(ind,:);

end
onset=[];onset=re_startpoint_sd;%行为
scatter(onset(:,1),(onset(:,2)-10*60-1)/60,12,[0 0 0],'filled');hold on;

a_event_onset=[a_event_onset1 a_event_onset2];m=[];mx=[];sem=[];
for i=trial.hab+1:trial.hab+trial.acq
    ind=find(a_event_onset(:,1)==i & a_event_onset(:,2)<frame.cs_end & a_event_onset(:,2)>=frame.cs_start);
    if ~isempty(ind)
        m(i-trial.hab)=mean((a_event_onset(ind,2)-frame.cs_start)*0.4);
        sem(i-trial.hab)=std((a_event_onset(ind,2)-frame.cs_start)*0.4);%/sqrt(length(ind));
        mx(i-trial.hab)=i;
    else
        m(i-trial.hab)=[]; mx(i-trial.hab)=[];
    end
end
e=errorbar(mx,m,sem,'-*r');hold on

activities_event_acq_serrorbar.Olfactory=[mx;m;sem];
save(['H:\20180623_fish1\分区\' 'activities_event_acq_serrorbar'],'activities_event_acq_serrorbar');







