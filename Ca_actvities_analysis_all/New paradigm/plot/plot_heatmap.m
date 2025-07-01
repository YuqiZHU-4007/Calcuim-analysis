function plot_heatmap(a,event,fs,frame,trial,frameb,re_startpoint)

if nargin<2
    event=[];
end

b=0;
linewidth=3;
aa=reshape(a,frame.per_cycle,trial.total)';
x=([1:frame.per_cycle]-frame.cs_start)*fs.ca;
y=1:trial.total;

set(gcf,'Units','normalized', 'position',[0.01,0.05,0.96,0.9]) ;
set(axes,'position',[0.03,0.05,0.95,0.9]);
%set(gca,'YTick',0:1:trial.total);
imagesc(x,y,aa);hold on %[min(min_activity_dfdf) max(max_activity_dfdf)];[-1 1]
line(([frame.cs_start frame.cs_start]-frame.cs_start)*fs.ca,[trial.hab(2) trial.test(3)]+0.5,'color',[1 0.6 0.78],'linestyle','--','LineWidth',linewidth);hold on;
line(([frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,[trial.hab(2) trial.test(3)]+0.5,'color',[1 0.6 0.78],'linestyle','--','LineWidth',linewidth);hold on;
line(([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,[trial.hab(3)+1-0.5 trial.acq(3)+0.5],'color','r','linestyle','-','LineWidth',linewidth);hold on;

%behavior
linewidthb=1.3;
onset=[];onset=re_startpoint;
line(([onset(:,2) onset(:,2)]'-frameb.cs_start)/fs.behavior,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[0 0 0],'LineWidth',linewidthb);hold on;
%calcuim event
if ~isempty(event)
    onset=[];
    if ~isempty(event.ind)
        onset(:,1)=ceil(event.ind/frame.per_cycle);
        onset(:,2)=mod(event.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        line(([onset(:,2) onset(:,2)]'-frame.cs_start)*fs.ca,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[1 1 1],'LineWidth',linewidthb);hold on;
        %     onset=[];onset(:,1)=ceil(event.end_ind/frame.per_cycle);%event end
        %     onset(:,2)=mod(event.end_ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        %     line(([onset(:,2) onset(:,2)]'-frame.cs_start)*fs_ca,[onset(:,1)-0.5 onset(:,1)+0.5]','color',[0.5 0.5 0.5],'LineWidth',linewidthb);hold on;
    end
end
set(gca,'FontSize',13);