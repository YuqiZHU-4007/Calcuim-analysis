function event_Zonal_Statistics_for_new_paradigm_20181016
global fs
global time
global frame
global frameb
global trial
global re_startpoint
global startpoint
global y_3sd

set(0,'defaultfigurecolor','w');
[actname,actpath]=uigetfile('G:\data\.mat','activities');
load([actpath,actname]);
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
%load([behaviorpath,behaviorname]);
outputpath=uigetdir('G:\data','outputpath');
[loctname,locpath]=uigetfile('G:\data\.xls','location');

%set para.
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);

[~,sheet,~]=xlsfinfo([locpath,loctname]);sheet

for kk=1:size(sheet,2)
    disp(sheet{kk});
    num=xlsread([locpath,loctname],sheet{kk});
    if isempty(num)
        continue;warning([sheet{kk} 'is empty']);
    end
    outpath=checkpath([outputpath '\' sheet{kk}]);
    acti=[];a_event={};a_event_rep=[];
    for ii=1:size(num,1)
        acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
        a_event{ii,1}=activities_event_preCS{num(ii,1),1}{num(ii,2),1}';
        %a_event_rep=[a_event_rep a_event{ii,1}.ind];
    end
    
    %全画
    checkpath([outpath '\heatmap']);
    checkpath([outpath '\raw trace']);
    checkpath([outpath '\raw trace-sp\']);
    %checkpath([outpath '\raw trace-sp-aft\']);
    
    for ii=1:size(num,1)
        event=activities_event_preCS{num(ii,1),1}{num(ii,2),1};onset=[];
        h1=figure;
        plot_heatmap(acti(:,ii),event,fs,frame,trial,frameb,re_startpoint);
        title([num2str(num(ii,1)) 'layer-' num2str(num(ii,2))]);xlabel('time(s)');ylabel('trial num.');colorbar;
        saveas(h1,[outpath '\heatmap\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h1);
        
        %hab;acq;test的raw trace
        h2=figure;title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Hab']);
        plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,1); 
        saveas(h2,[outpath '\raw trace\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h2);
        
        %spon bef and aft
        h3=figure;title([num2str(num(ii,1)) '-layer' num2str(num(ii,2))]);
        plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,2); 
        saveas(h3,[outpath '\raw trace-sp\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h3);
        
    end
    
    
    %总体
    outpathtotal=checkpath([outpath '\total']);
    %heatmap
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.9]);
    imagesc(acti');%colorbar;
    hold on;
    xlabel('frame');
    ylabel('cell number');
    set(gca,'YTick',1:size(num,1));
    saveas(h,[outpathtotal '\' sheet{kk} '-heatmap' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-heatmap']);
    close(h);
    %均值
    y=mean(acti,2);event=getevent(y,trial,frame);
    %heatmap
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    plot_rawtrace(y,[1 frame.per_cycle*trial.total],[],fs,frame,trial,startpoint);legend('off');
    saveas(h,[outpathtotal '\' sheet{kk} '-mean' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-mean']);
    close(h);
    %all trace
    h=figure;
    plot_rawtrace_trials(acti,event,fs,frame,trial,startpoint,1);
    saveas(h,[outpathtotal '\' sheet{kk} '-all trace' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-all trace']);
    %mean-heatmap
    h1=figure;
    plot_heatmap(y,event,fs,frame,trial,frameb,re_startpoint);
    title([sheet{kk} '-mean']);xlabel('time(s)');ylabel('trial num.');colorbar;
    saveas(h1,[outpathtotal '\' sheet{kk} '-mean-heatmap' '.tif']);
    savefig(h1,[outpathtotal '\' sheet{kk} '-mean-heatmap']);
    close(h1);
    %hab;acq;test的raw trace
    h2=figure;title([sheet{kk} '-mean' ' Hab']);
    plot_rawtrace_trials(y,event,fs,frame,trial,startpoint,1);
    saveas(h2,[outpathtotal '\' sheet{kk} '-mean-raw trace' '.tif']);
    savefig(h2,[outpathtotal '\' sheet{kk} '-mean-raw trace']);
    close(h2);
    
    %onset随时间变化
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    [summ1,summ2]=onset_errorbar(a_event,fs,trial,frame,re_startpoint);
    activities_event_acq_errorbar{kk,1}=[summ1.mx;summ1.m;summ1.sem]';
    activities_event_acq_errorbar{kk,2}=sheet{kk};
    weight=[[trial.acq(1):-1:1]/trial.acq(1)];
    activities_event_acq_errorbar{kk,3}=[summ2.ind;summ2.m;summ2.ind.*weight]';
    saveas(h,[outpathtotal '\' sheet{kk} '-errorbar onset' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-errorbar onset']);
    close(h);    
end
save([outputpath '\activities_event_acq_errorbar'],'activities_event_acq_errorbar');

% for kk=1:9
% activities_event_acq_serrorbar{kk,1}=[activities_event_acq_serrorbar{kk,1}]';
% end

% %onset
% h=figure;
% set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
% set(axes,'position',[0.054,0.08,0.90,0.6]);
% xlabel('frames');ylabel('cell number');
% patch([[1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
%     repmat([0 0 size(num,1) size(num,1)],ceil(trial.total/2),1)',...
%     [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
% patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
%     repmat([0 0 size(num,1) size(num,1)],trial.test(3)-trial.hab(2)+1,1)',...
%     'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
% 
% line([frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3));frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3))],...
%     [0 size(num,1)],'color','r','linewidth',1);hold on;
% % line([1:frame.per_cycle:frame.per_cycle*trial.total;1:frame.per_cycle:frame.per_cycle*trial.total],...
% %     [0 size(num,1)],'color',[0.5 0.5 0.5],'linestyle','--');hold on;
% 
% for ii=1:size(num,1)
%     if ~isempty(a_event{ii,1}.ind)
%         line([a_event{ii,1}.ind;a_event{ii,1}.ind],[ii-0.3 ii+0.3],'color',[0 0 0],'linewidth',1);hold on;
%     end
%     %scatter(a_event{ii,1}.ind,repmat(ii,length(a_event{ii,1}.ind),1),8,'k','filled');hold on;
% end
% ylim([0 size(num,1)]);xlim([1 frame.per_cycle*trial.total]);
% set(gca,'YTick',0:1:size(num,1));
% saveas(h,[outpathtotal '\' sheet{kk} '-onset' '.tif']);
% savefig(h,[outpathtotal '\' sheet{kk} '-onset']);
% close(h);
% 
% %onset统计图
% h=figure;
% set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
% set(axes,'position',[0.054,0.08,0.90,0.6]);
% xlabel('frames');ylabel('number of events');
% patch([[1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
%     [1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
%     repmat([0 0 size(num,1) size(num,1)],ceil(trial.total/2),1)',...
%     [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
% patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
%     [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
%     repmat([0 0 size(num,1) size(num,1)],trial.test(3)-trial.hab(2)+1,1)',...
%     'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
% line([frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3));frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3))],...
%     [size(num,1)-1 size(num,1)],'color','r','linewidth',1);hold on;
% h2=histogram(a_event_rep,frame.per_cycle*trial.total/1,'facecolor','k');hold on
% % line([frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq);frame.us_start+trial.hab*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.hab+trial.acq)],...
% %     [0 size(num,1)],'color','r','linewidth',1);hold on;
% line([frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3));frame.us_start+trial.hab(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*(trial.acq(3))],...
%     [0 2],'color','r','linewidth',1);hold on;
% xlim([0 frame.per_cycle*trial.total]);ylim([1 size(num,1)])
% %y=sum(a_event_rep,2);
% saveas(h,[outpathtotal '\' sheet{kk} '-statistic onset' '.tif']);
% savefig(h,[outpathtotal '\' sheet{kk} '-statistic onset']);
% close(h);

