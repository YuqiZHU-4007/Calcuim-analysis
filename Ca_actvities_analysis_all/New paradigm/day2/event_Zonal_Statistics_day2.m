function event_Zonal_Statistics_day2

set(0,'defaultfigurecolor','w');
frame.per_cycle = 50;
frame.cs_start = 26;
frame.cs_end = 38;
frame.us_start = 36;
%1:trial num;2:start trial;3:end trial
trial.spon_bef=[6*3 1 6*3];
trial.test = [6 trial.spon_bef(3)+1 trial.spon_bef(3)+6];
trial.total =trial.spon_bef(1)+trial.test(1);
fs_ca=0.4;
fs_behavior=60;
frameb.cs_start=10*fs_behavior+1;

%[filename,filepath]=uigetfile('H:\','.*');
[inputname,inputpath]=uigetfile('G:\data\.xlsx','location');
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
[actname,actpath]=uigetfile('G:\data\.mat','activities_aft_process');
load([behaviorpath,behaviorname]);
load([actpath,actname]);
outputpath=uigetdir('G:\data','outputpath');


[~,sheet,~]=xlsfinfo([inputpath,inputname]);

for kk=1:size(sheet,2)
    disp(sheet{kk});
    outpath=checkpath([outputpath '\' sheet{kk}]);
    num=xlsread([inputpath,inputname],sheet{kk});
    
    acti=[];a_event={};a_event_rep=[];
    for ii=1:size(num,1)
        acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
        a_event{ii,1}=activities_event_preCS{num(ii,1),1}{num(ii,2),1}';
        a_event_rep=[a_event_rep a_event{ii,1}.ind];
    end
    
    %全画
    checkpath([outpath '\heatmap-test']);
    checkpath([outpath '\raw trace-sp']);
    checkpath([outpath '\raw trace-test\']);
    %checkpath([outpath '\raw trace-sp-aft\']);
    for ii=1:size(num,1)
        %figure,plot(activities_new{1}(:,1))
        %aa=acti(:,ii);
        lab='dfdf'; b=0;
        aa=reshape(acti(trial.spon_bef(3)*frame.per_cycle+1:trial.test(3)*frame.per_cycle,ii),frame.per_cycle,trial.test(1))';
        x=([1:frame.per_cycle]-frame.cs_start)*0.4;
        y=1:trial.test(1);
        h1=figure;
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        set(axes,'position',[0.03,0.05,0.95,0.9]);
        %set(gca,'YTick',0:1:trial.total);
        imagesc(x,y,aa);hold on %[min(min_activity_dfdf) max(max_activity_dfdf)];[-1 1]
        title([num2str(num(ii,1)) 'layer-' num2str(num(ii,2))]);xlabel('time(s)');ylabel('trial num.');colorbar;
        line(([frame.cs_start frame.cs_start]-frame.cs_start)*0.4,[1-0.5 trial.test(1)+0.5],'color',[0 0 0],'linestyle','--');hold on;
        line(([frame.cs_end frame.cs_end]-frame.cs_start)*0.4,[1-0.5 trial.test(1)+0.5],'color','k','linestyle','--');hold on;
        
        %behavior
        onset=[];onset=re_startpoint_sd;
        line(([onset(:,2) onset(:,2)]'-frameb.cs_start)/fs_behavior,[onset(:,1)-0.5 onset(:,1)+0.5]'-trial.spon_bef(3),'color',[1 1 0],'linewidth',1.2);hold on;
        %event
        event=activities_event_preCS{num(ii,1),1}{num(ii,2),1};onset=[];
        onset(:,1)=ceil(event.ind/frame.per_cycle);
        onset(:,2)=mod(event.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]'-trial.spon_bef(3),'color',[1 1 1],'linewidth',1.2);hold on;
        onset=[];onset(:,1)=ceil(event.end_ind/frame.per_cycle);
        onset(:,2)=mod(event.end_ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        line(([onset(:,2) onset(:,2)]'-frame.cs_start)*0.4,[onset(:,1)-0.5 onset(:,1)+0.5]'-trial.spon_bef(3),'color',[0.5 0.5 0.5],'linewidth',1.2);hold on;
        
        saveas(h1,[outpath '\heatmap-test\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h1);
        
        h2=figure;%spon
        aa=acti(:,ii);
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        plot([1:frame.per_cycle*trial.total]*0.4,aa);hold on;xlim([1 trial.spon_bef(3)*frame.per_cycle]*0.4);hold on;ylim([min(aa(:)) max(aa(:))]);
        scatter(event.ind*0.4,aa(event.ind),9,'k','filled');hold on;
        %behavior
        onset=[];onset=startpoint_sd;
        scatter(([onset]'-frameb.cs_start)/fs_behavior+frame.cs_start*0.4,ones(length(onset),1)*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),9,'r','filled');hold on;
        title([num2str(num(ii,1)) '-layer' num2str(num(ii,2)) ' Spon.']);xlabel('time(s)');ylabel(lab);
        saveas(h2,[outpath '\raw trace-sp\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h2);
        
        h3=figure;%test
        aa=acti(:,ii);
        set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.85]) ;
        plot([1:frame.per_cycle*trial.total]*0.4,aa);hold on;xlim([trial.spon_bef(3)*frame.per_cycle trial.test(3)*frame.per_cycle]*0.4);hold on;ylim([min(aa(:)) max(aa(:))]);
        plot(repmat([frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.total],2,1)'*0.4,[min(aa(:)) max(aa(:))],'g');hold on
        plot(repmat([frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        plot(repmat([frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)],2,1)'*0.4,[min(aa(:)) max(aa(:))],'color',[1 0.5 0],'LineStyle','--');hold on
        scatter(event.ind*0.4,aa(event.ind),9,'k','filled');hold on; 
        %behavior
        onset=[];onset=startpoint_sd;
        scatter(([onset]'-frameb.cs_start)/fs_behavior+frame.cs_start*0.4,ones(length(onset),1)*(min(aa(:))+0.1*(max(aa(:))-min(aa(:)))),9,'r','filled');hold on;
        title(['Test']);xlabel('time(s)');ylabel(lab);
        saveas(h3,[outpath '\raw trace-test\' num2str(num(ii,1)) 'layer-' num2str(num(ii,2)) '.tif']);
        close(h3)
    end
    
    %总体
    outpathtotal=checkpath([outpath '\total']);
    %hotmap
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
    y=mean(acti,2);
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    patch([[frame.per_cycle*trial.spon_bef(3)+1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
        repmat([min(y) min(y) max(y) max(y)],ceil(trial.test(1)/2),1)',...
        [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
    patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
        repmat([min(y) min(y) max(y) max(y)],trial.test(1),1)',...
        'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
    
    ylim([min(y) max(y)]);xlim([1 frame.per_cycle*trial.total]);xlabel('frames');ylabel('dfdf');
    plot(mean(acti,2),'k');hold on;
    % line([frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total;frame.cs_start:frame.per_cycle:frame.per_cycle*trial.total],...
    %      [min(y) max(y)],'color','b');hold on;
    % line([frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total;frame.cs_end:frame.per_cycle:frame.per_cycle*trial.total],...
    %      [min(y) max(y)],'color','b');hold on;
    saveas(h,[outpathtotal '\' sheet{kk} '-mean' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-mean']);
    close(h);
    
    
    %onset
    h=figure;
    set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    set(axes,'position',[0.054,0.08,0.90,0.6]);
    xlabel('frames');ylabel('cell number');
     patch([[frame.per_cycle*trial.spon_bef(3)+1:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+frame.per_cycle:frame.per_cycle*2:frame.per_cycle*trial.total]'...
        [frame.per_cycle*trial.spon_bef(3)+1:frame.per_cycle*2:frame.per_cycle*trial.total]']',...
        repmat([min(y) min(y) max(y) max(y)],ceil(trial.test(1)/2),1)',...
        [0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facecolor',[0.5 0.5 0.5],'edgealpha',0.1,'facealpha',0.2);hold on %trial
    patch([[frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_end+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]'...
        [frame.cs_start+trial.spon_bef(3)*frame.per_cycle:frame.per_cycle:frame.per_cycle*trial.test(3)]']',...
        repmat([min(y) min(y) max(y) max(y)],trial.test(1),1)',...
        'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
    for ii=1:size(num,1)
        line([a_event{ii,1}.ind;a_event{ii,1}.ind],[ii-0.3 ii+0.3],'color',[0 0 0],'linewidth',1);hold on;
        %scatter(a_event{ii,1}.ind,repmat(ii,length(a_event{ii,1}.ind),1),8,'k','filled');hold on;
    end
    ylim([0 size(num,1)]);xlim([1 frame.per_cycle*trial.total]);
    set(gca,'YTick',0:1:size(num,1));
    saveas(h,[outpathtotal '\' sheet{kk} '-onset' '.tif']);
    savefig(h,[outpathtotal '\' sheet{kk} '-onset']);
    close(h);
    
end





