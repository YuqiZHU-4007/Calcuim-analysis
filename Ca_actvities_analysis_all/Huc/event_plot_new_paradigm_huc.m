function event_plot_new_paradigm_huc(locpath,activities_preCS_dfdf_aftcorrect,activities_event_preCS,outputpath)
global fs
global frame
global frameb
global trial
global re_startpoint
global startpoint

if ~isempty(find(locpath))
    [~,sheet,~]=xlsfinfo(locpath);
else
    sheet={'total'};
end
num=[];activities_event_acq_errorbar=cell(size(sheet,2),1);
for kk=5:6%1:size(sheet,2)
    disp(sheet{kk});
    if ~isempty(find(locpath))
        num=xlsread(locpath,sheet{kk});
    else
        num=1:size(activities_preCS_dfdf_aftcorrect,2);
    end
    if isempty(num)
        continue;warning([sheet{kk} 'is empty']);
    end
    outpath=checkpath([outputpath '\' sheet{kk}]);
    acti=[];a_event=cell(size(num,1),1);%a_event_rep=[];
    acti=activities_preCS_dfdf_aftcorrect(:,num(:,2));
    for ii=1:size(num,1)
    a_event{ii,1}=activities_event_preCS{num(ii,2)};
    end
    meanacti(:,kk)=mean(acti,2);
    acti_name(:,kk)=string(sheet{kk});
    %全画
    checkpath([outpath '\heatmap']);
    checkpath([outpath '\raw trace']);
    %checkpath([outpath '\raw trace-sp-aft\']);
    
%     for ii=1:ceil(size(num,1)/500):size(num,1)
%         event=activities_event_preCS{num(ii,2),1};onset=[];
%         h1=figure;
%         plot_heatmap(acti(:,ii),event,fs,frame,trial,frameb,re_startpoint);
%         title([num2str(num(ii,1)) 'layer-' num2str(num(ii,2))]);xlabel('time(s)');ylabel('trial num.');colorbar;
%         saveas(h1,[outpath '\heatmap\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
%         close(h1);
%         
%         %hab;acq;test的raw trace
%         h2=figure;title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Hab']);
%         plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,1);
%         saveas(h2,[outpath '\raw trace\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
%         close(h2);
%         
%         %spon bef and aft
%         if trial.spon_bef(1)~=0 && trial.spon_aft(1)~=0
%             checkpath([outpath '\raw trace-sp\']);
%             h3=figure;title([num2str(num(ii,1)) '-layer' num2str(num(ii,2))]);
%             plot_rawtrace_trials(acti(:,ii),event,fs,frame,trial,startpoint,2);
%             saveas(h3,[outpath '\raw trace-sp\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
%             close(h3);
%         end
%         
%     end
    
    %总体
    outpathtotal=checkpath([outpath '\total']);
    y=mean(acti,2);event=getevent(y,trial,frame);
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
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    plot_rawtrace(y,[1 frame.per_cycle*trial.total],[],fs,frame,trial,startpoint);legend('off');
    saveas(h,[outpathtotal '\' sheet{kk} '-mean' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-mean']);
    close(h);
%     %all trace
%     h=figure;
%     plot_rawtrace_trials(acti,event,fs,frame,trial,startpoint,1);
%     saveas(h,[outpathtotal '\' sheet{kk} '-all trace' '.tif']);
%     savefig(h,[outpathtotal '\' sheet{kk} '-all trace']);
%     close(h);
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
    %mean_acq
    hh1=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    plot_meantrace_acq(y,fs,trial,frame);
    saveas(hh1,[ outpath '\'  'meantrace_acq' '.tif']);
    savefig(hh1,[ outpath '\'  'meantrace_acq']);
    close(hh1);
    %mean_hab_test
    hh2=figure;
    plot_mean_hab_test(y,fs,trial,frame,frameb,startpoint);
    saveas(hh2,[ outpath '\'  'trace_hab_test' '.tif']);
    savefig(hh2,[ outpath '\'  'trace_hab_test']);
    close(hh2);
    
    %all traces sepplot
    ind=1:size(num,1);
    if ~isempty(ind)
        str=strcat({char(string(num(ind,1)))},{char(char('-')*ones(size(num(ind,1))))},{char(string(num(ind,2)))});
        ind_frames=1:trial.total*frame.per_cycle;
        h=figure;
        set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.88]) ;
        plot_sepplot(ind_frames,acti(:,ind),str,frame,trial);hold on;
        saveas(h,[ outpath '\'  'all traces sepplot' '.tif']);
        savefig(h,[ outpath '\'  'all traces sepplot']);close(h);
        %figure,Correlational_analysis(acti(:,ind));
        h=figure;subplot(1,2,1);Correlational_analysis(acti((trial.hab(2)-1)*frame.per_cycle+1:trial.hab(3)*frame.per_cycle,ind));
        subplot(1,2,2);Correlational_analysis(acti((trial.test(2)-1)*frame.per_cycle+1:trial.test(3)*frame.per_cycle,ind));
        saveas(h,[ outpath '\'  'all traces corr_hab_test' '.tif']);
        savefig(h,[ outpath '\'  'all traces corr_hab_test']);
        close(h);
    end
end
save([outputpath '\activities_event_acq_errorbar'],'activities_event_acq_errorbar','meanacti','acti_name');

% for kk=1:9
% activities_event_acq_serrorbar{kk,1}=[activities_event_acq_serrorbar{kk,1}]';
% end




