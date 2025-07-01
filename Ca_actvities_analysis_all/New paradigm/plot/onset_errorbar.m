function [summ1,summ2]=onset_errorbar(a_event,fs,trial,frame,re_startpoint)
%onset随时间变化

linewidth=1.2;scattersize=20;
x=trial.hab(2):trial.test(3);
%CS bar
color=[0.5 0.5 0.5];
p1=patch([trial.hab(2)  trial.test(3) trial.test(3) trial.hab(2)],([frame.cs_start frame.cs_start frame.cs_end frame.cs_end]-frame.cs_start)*fs.ca,...
    color,'edgecolor',color,'facealpha',0.2,'edgealpha',0.25);hold on;
%US bar
l1=line([trial.acq(2) trial.acq(3)],([frame.us_start frame.us_start]-frame.cs_start)*fs.ca,'color','r','linewidth',linewidth);hold on;

a_event_onset1=[];
a_event_onset2=[];
for ii=1:size(a_event,1)
    onset=[];
    if ~isempty(a_event{ii,1}.ind)
        onset(:,1)=ceil(a_event{ii,1}.ind/frame.per_cycle);
        onset(:,2)=mod(a_event{ii,1}.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
        s2=scatter(onset(:,1),(onset(:,2)-frame.cs_start)*fs.ca,scattersize);hold on;
        onset(find(onset(:,2)>=frame.cs_end | onset(:,2)<frame.cs_start),:)=[];
        if ~isempty(onset)
            tb=tabulate(onset(:,1));ind=find(tb(:,2)>1);
            if ~isempty(ind)
                for jj=1:length(ind)
                    onset(find(onset(:,1)==tb(ind(jj))),:)=[];
                    onset(end+1,:)=[tb(ind(jj)) min(onset(find(onset(:,1)==tb(ind(jj))),2))];
                end
            end
            a_event_onset1=[a_event_onset1;onset(:,1)];
            a_event_onset2=[a_event_onset2;onset(:,2)];
        end
        %         ind=setdiff(x,intersect(onset(:,1)',x))';
        %         onset(size(onset,1)+1:size(onset,1)+length(ind),:)=[ind zeros(size(ind))];[~,ind]=sort(onset(:,1));onset=onset(ind,:);
    end
end
onset=[];onset=re_startpoint;%行为
s1=scatter(onset(:,1),(onset(:,2)-((frame.cs_start-1)*fs.ca)*fs.behavior-1)/fs.behavior,scattersize,[0 0 0],'^','filled');hold on;

a_event_onset=[a_event_onset1,a_event_onset2];summ1=[];summ2=[];
if ~isempty(a_event_onset)
    for i=trial.hab(3)+1:trial.acq(3)
        ind=find(a_event_onset(:,1)==i);
        if ~isempty(ind)
            summ1.m(i-trial.hab(3))=mean((a_event_onset(ind,2)-frame.cs_start))*fs.ca;
            summ1.sem(i-trial.hab(3))=std((a_event_onset(ind,2)-frame.cs_start))*fs.ca;%/sqrt(length(ind));
            summ1.mx(i-trial.hab(3))=i;
        else
            summ1.m(i-trial.hab(3))=0; summ1.sem(i-trial.hab(3))=0;summ1.mx(i-trial.hab(3))=0;
        end
        if ~isempty(ind) &&  mean(a_event_onset(ind,2))<frame.us_start
            summ2.m(i-trial.hab(3))=mean(a_event_onset(ind,2));
            summ2.ind(i-trial.hab(3))=true;
            summ1.mx(i-trial.hab(3))=i;
        else
            summ2.m(i-trial.hab(3))=0; summ2.ind(i-trial.hab(3))=0;summ2.mx(i-trial.hab(3))=i;
        end
    end
    
    summ1.m(summ1.mx==0)=[];summ1.sem(summ1.mx==0)=[];summ1.mx(summ1.mx==0)=[];
    e=errorbar(summ1.mx,summ1.m,summ1.sem,'-*r','LineWidth',1.5);hold on
    xlim([x(1) x(end)]);ylim(([1 frame.per_cycle+1]-frame.cs_start)*fs.ca);
    xlabel('trial number');ylabel('Time(s)')
    legend([p1,l1,s1,s2],{'CS','US','Behavior event onset','Ca event onset'},'position',[0.793726744588737,0.737768636100846,0.157894733348531,0.117283947346796]);
    set(gca,'FontSize',13);box on;
end