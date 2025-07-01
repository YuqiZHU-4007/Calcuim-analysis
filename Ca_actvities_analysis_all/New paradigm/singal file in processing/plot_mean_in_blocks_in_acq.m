clc;clear all;close all;
set(0,'defaultfigurecolor','w');

[actname,actpath]=uigetfile('G:\data\.mat','activities aft process');
load([actpath,actname]);
[loctname,locpath]=uigetfile([actpath '.xls'],'location');
[behaviorname,behaviorpath]=uigetfile([actpath '.mat'],'behavior');
%load([behaviorpath,behaviorname]);
outputpath=uigetdir(actpath,'outputpath');

%set para.
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd]=setpara([behaviorpath,behaviorname]);

[~,sheet,~]=xlsfinfo([locpath,loctname]);sheet

area_CS_acq_1st_per_block = [];
area_CS_acq_mean_per_block = [];
area_CS_hab_test = [];
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
    end
    
    aa=mean(acti,2);
    meanacti(:,kk)=mean(acti,2);
    acti_name(:,kk)=string(sheet{kk});
    hh1=figure;
    %set(gcf,'Units','normalized', 'position',[0.05,0.05,0.9,0.6]) ;
    [area_CS_acq_block,area_CS_acq_1st_block] = plot_meantrace_acq(aa,fs,trial,frame);
    area_CS_acq_1st_per_block = [area_CS_acq_1st_per_block;area_CS_acq_1st_block];
    area_CS_acq_mean_per_block = [area_CS_acq_mean_per_block;area_CS_acq_block];
    saveas(hh1,[ outpath '\'  'meantrace_acq' '.tif']);
    savefig(hh1,[ outpath '\'  'meantrace_acq']);
    close(hh1);
    
    
    hh2=figure;
    [rawacti, area_CS_hab_tst] = plot_mean_hab_test(aa,fs,trial,frame,frameb,startpoint);
    area_CS_hab_test = [area_CS_hab_test;area_CS_hab_tst];
    saveas(hh2,[ outpath '\'  'trace_hab_test' '.tif']);
    savefig(hh2,[ outpath '\'  'trace_hab_test']);
    close(hh2);
end
% %bb_20180623fish1=bb;
% aa=mean(acti,2);%aa=zscore(aa);
% bb=bb_20180623fish2;%bb=zscore(bb);
% for ii=1:trial.acq_block_num
%     board=[min(min(aa),min(bb)) max(max(aa),max(bb))];
%     figure,
%     line([frame.us_start frame.us_start],[board(1) board(2)],'color',[1 0 0],'linestyle','--');hold on
%     patch([frame.cs_start...
%         frame.cs_end...
%         frame.cs_end'...
%         frame.cs_start],...
%         [board(1) board(1) board(2) board(2)],...
%         'g','edgecolor','g','facecolor','g','edgealpha',0.2,'facealpha',0.3);hold on %CS on-off
%     ind1=trial.acq_block_trial*(ii-1)+1+trial.hab(3);
%     ind2=trial.acq_block_trial*(ii-1)+1+trial.hab(3)+trial.acq_block_trial-1;
%     %subplot(trial.acq_block_num,1,ii),
%      plot(mean(reshape(bb((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)'));hold on
%      plot(mean(reshape(aa((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)'));hold on
%
%      mma=mean(reshape(aa((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)');
%      inda{ii}=find(mma>mean(mma(20:frame.cs_start-1))+1*std(mma(1:frame.cs_start-1)));inda{ii}=min(inda{ii}(find(inda{ii}>=frame.cs_start)));
%      mmb=mean(reshape(bb((ind1-1)*frame.per_cycle+1:ind2*frame.per_cycle),frame.per_cycle,trial.acq_block_trial)');
%      indb{ii}=find(mmb>mean(mmb(20:frame.cs_start-1))+1*std(mmb(1:frame.cs_start-1)));indb{ii}=min(indb{ii}(find(indb{ii}>=frame.cs_start)));
%      maxium(ii,:)=[max(mma(frame.cs_start:frame.us_start-1))...
%          max(mmb(frame.cs_start:frame.us_start-1))];
%     legend('us','cs','0623fish2','1123fish4');
%     ylim([board(1) board(2)]);
% end
