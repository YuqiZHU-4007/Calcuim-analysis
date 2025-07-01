function plot_realtime_trace(trialind)
set(0,'defaultfigurecolor','k');

[actname,actpath]=uigetfile('G:\data\.mat','activities aft process');
load([actpath,actname]);
[loctname,locpath]=uigetfile([actpath '.xls'],'location');
[behaviorname,behaviorpath]=uigetfile([actpath '.mat'],'behavior');
%load([behaviorpath,behaviorname]);
outputpath=uigetdir(actpath,'outputpath');

%set para.
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);
trialind=[trial.acq(2) trial.acq(3)];%%»­

[~,sheet,~]=xlsfinfo([locpath,loctname]);

for kk=1:size(sheet,2)
    disp(sheet{kk});
    num=xlsread([locpath,loctname],sheet{kk});
    if isempty(num)
        continue;warning([sheet{kk} 'is empty']);
    end
    acti=[];
    for ii=1:size(num,1)
        acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
    end
    meantrace(:,kk)=mean(acti,2);
end


imagepath=uigetdir(actpath,'imagepath');
listdir = dir([imagepath '/*.JPG']);
imagebehavior=zeors(401,1031,length(listdir));
if ~isempty(listdir)
    for ii = 1:length(listdir)
        filename = listdir(ii).name;
        imagebehavior(:,:,ii)=imread([imagepath '/' filename]);
    end
end

%real time plot
tracecolor(:,1)=string(["Raphe";"Pr-5HT";"OB";"LC";"Pr-DA";"PT";"Pallium";"Sub pallium";"CG"]);
tracecolor(:,2)=num2str([[0 1 0];[1 0 0];[0 1 1];[0.75 0 0.75];[1 1 1];[0.45 0.26 0.26];[0.15 0.23 0.37];[0 0 0];[1 1 1]]);%Raphe;PT,GC,Pr5HT,PrDA,Pallium,OB,Subpallium
frameind=[(trialind(1)-1)*frame.per_cycle+1:trialind(2)*frame.per_cycle];
bin=6;
frmaeind_beh=[(trialind(1)-1)*frameb.per_cycle+1:1:trialind(2)*frameb.per_cycle]-trial.hab(3)*frameb.per_cycle;
showtrace=meantrace(:,[4 1 8]);
tracename={'Pr','Raphe','OB'};%string(sheet{[1 4 8]});
%tracename(find(showtrace(frame.cs_start,:)==0 & showtrace(frame.us_start,:)==0))=[];
%showtrace(:,find(showtrace(frame.cs_start,:)==0 & showtrace(frame.us_start,:)==0))=[];
edge=[min(min(showtrace(frameind,:))) max(max(showtrace(frameind,:)))]+[-0.05 0.05]*(max(max(showtrace(frameind,:)))-min(min(showtrace(frameind,:))));
checkpath([outputpath '\' 'realtimetrace2']);
t=1;
h=figure;
F(1) = getframe(h);
F(2:length(frmaeind_beh)/bin)=F(1);close(h);
zz=1;
for t = 1:bin:length(frmaeind_beh)
    h=figure;
    set(gcf,'position',[1,1,1920,1500],'color','k');
    %image
    %subplot(4,1,1)
    set(axes,'position',[ 0.06,0.62,0.62,0.34],'color','w');
    I=imshow(imagebehavior(:,:,t));hold on;axis normal
    %CS
    if single(mod(t,frameb.per_cycle))>=frameb.cs_start &&  single(mod(t,frameb.per_cycle))<frameb.cs_end &&  single(mod(t,frameb.per_cycle))~=frameb.us_start
        text(79.7019,49.1861,'CS','color',[1 1 1],'fontsize',50,'HorizontalAlignment','center');
    
    elseif single(mod(t,frameb.per_cycle))==frameb.us_start || single(mod(t,frameb.per_cycle))==frameb.us_start+1 %US 
         text(79.7019,49.1861,'US','color',[1 0 0],'fontsize',50,'HorizontalAlignment','center');
    end
      
    %trace
    frameinds=frameind(1)+[1:ceil(t/(fs.ca*fs.behavior))]-1;
    times =frameinds*fs.ca;%s
    %     zz=0;
    for kk=1:3%size(showtrace,2)
        %         if isempty(find(showtrace(frame.cs_start,kk)~=0 & showtrace(frame.us_start,kk)~=0))
        %             continue;warning([tracename(kk) 'is empty']);
        %         end
        %         zz=zz+1;
        %         subplot(length(find(showtrace(frame.cs_start,:)~=0)),1,zz);
        switch kk
            case 1
                %outerposition=[0.0021,0.4,1.0,0.16];1
                position=[0.06,0.4405,0.85,0.15];
            case 2
                %outerposition=[0.0016,0.2195,1.0000,0.1604];2
                position=[0.06,0.2805,0.85,0.15];
            case 3
                %outerposition=[0.001,0.0098,1.0000,0.201];3
                position=[ 0.06,0.12,0.85,0.15];
        end
        %subplot(4,1,kk+1);
        set(axes,'position',position,'color','k');%'outerposition',outerposition,
        y=showtrace(frameinds,kk);
        realtime_frames(times,y,[frameind(1)*fs.ca frameind(end)*fs.ca],str2num(tracecolor(kk,2)),edge,frame,trial,fs,startpoint);hold on;
        if kk==3
            set(gca,'xTick',[],'xcolor','w','ycolor','w','tickdir','out','ticklength',[0.005,0.001]);set(gca,'fontsize',20,'linewidth',4);
        else
            set(gca,'xTick',[],'xcolor','k','ycolor','w','tickdir','out','ticklength',[0.005,0.001]);set(gca,'fontsize',20,'linewidth',4);
        end
        ylabel(tracename{kk});
    end
    %behavior
    set(axes,'position',[ 0.06,0.07,0.85,0.03],'color','k');
    onset=[];onset=startpoint;onset(find(onset>times(end)*fs.behavior | onset<times(1)*fs.behavior))=[];onset;
    s2=scatter(onset/fs.behavior,ones(size(onset))*(0.0),250,[1 1 0],'^','filled');hold on;
    set(gca,'xtick',[],'ytick',[],'xcolor','k','ycolor','k','color','k');%set(gca,'fontsize',20,'linewidth',4);
    xlim([frameind(1)*fs.ca frameind(end)*fs.ca]);
    xlabel('Time(s)');
    set(gcf,'color','k','InvertHardcopy','off');
    F(zz) = getframe(h); % online show
    zz=zz+1;
    saveas(h,[outputpath '\' 'realtimetrace2' '\' num2str(t) '.jpg']);
    %print(h,[outputpath '\' 'realtimetrace2' '\' num2str(t)],'-djpeg');
    %exportsetupdlg
    close(h);
end
%F(201:end)=[];

obj = VideoWriter([outputpath '\' 'realtimetrace2.avi'],'Uncompressed AVI');
obj.FrameRate=60;%1/fs.ca*60;
open(obj)
writeVideo(obj,F);
close(obj);delete(obj);
