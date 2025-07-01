%求onset加权index；CR幅值加权；UR幅值加权，用来统计不同鱼的index和行为的关系

set(0,'defaultfigurecolor','w');
% frame.per_cycle = 50;
% frame.cs_start = 26;
% frame.cs_end = 38;
% frame.us_start = 36;
% 
% %6:自发
% %15：hab
% %6*4：acq
% %6:test
% %6:spon
% 
% %1:trial num;2:start trial;3:end trial
% trial.spon_bef=[6*3 1 6*3];
% trial.hab = [15 trial.spon_bef(3)+1 trial.spon_bef(3)+15];
% trial.acq_block_num=4;
% trial.acq_block_trial=6;
% trial.acq = [trial.acq_block_trial*trial.acq_block_num trial.hab(3)+1 trial.hab(3)+trial.acq_block_trial*trial.acq_block_num];
% trial.test =[6 trial.acq(3)+1 trial.acq(3)+6];
% trial.spon_aft=[6*3 trial.test(3)+1 trial.test(3)+6*3];
% trial.total =trial.spon_bef(1)+trial.hab(1)+trial.acq(1)+trial.test(1)+trial.spon_aft(1);
% fs_ca=0.4;
% fs_behavior=60;
%re_startpoint=re_startpoint_sd;startpoint=startpoint_sd;

[inputname,inputpath]=uigetfile('G:\data\.xlsx','location');
[behaviorname,behaviorpath]=uigetfile('G:\data\.mat','behavior');
[actname,actpath]=uigetfile('G:\data\.mat','activities_aft_process');
load([behaviorpath,behaviorname]);
load([actpath,actname]);
outputpath=uigetdir('G:\data','outputpath');


[~,sheet,~]=xlsfinfo([inputpath,inputname]);
statistic_CS_max=[];onset_logic_acq=[];statistic_US_max=[];
for kk=1:size(sheet,2)
    disp(sheet{kk});
    outpath=checkpath([outputpath '\' sheet{kk}]);
    num=xlsread([inputpath,inputname],sheet{kk});
    %num=xlsread(['H:\20180623_fish1\分区\分区_20180623_fish1.xlsx'],'Olfactory');
    if isempty(num)
        warning([sheet{kk} 'is empty']);continue;
    end
    acti=[];CS_max=[];
    for ii=1:size(num,1) %计算rois均值的onset
        acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
        %acti(51:100,ii)=(activities_new{num(ii,1),1}(51:100,num(ii,2))-mean(activities_new{num(ii,1),1}(75:100,num(ii,2))))/mean(activities_new{num(ii,1),1}(75:100,num(ii,2)));
    end
    aa=reshape(mean(acti,2),frame.per_cycle,trial.total)';
    base=getbaseline_cut3sd_strict_for_new_paradigm(aa,1:frame.cs_start-1,trial);
    m=base.rep_m;sd=base.rep_sd;
    baseline=m+2*sd;maxthres=m+3*sd;
    rasing_width=[1 100];
    decling_width=[4 100];%decling_width+rasing_width-1<8
    event=getevent_Trianglefitting_strict_for_new_paradigm(mean(acti,2),baseline,maxthres,rasing_width,decling_width,frame,trial,base);
    
    onset=[];
    onset(:,1)=ceil(event.ind/frame.per_cycle);
    onset(:,2)=mod(event.ind,frame.per_cycle);onset(find(onset(:,2)==0),:)=[];
    onset(find(onset(:,1)<trial.acq(2)),:)=[]; onset(find(onset(:,1)>trial.acq(3)),:)=[];%只留acq阶段CR的结果
    onset(find(onset(:,2)<frame.cs_start),:)=[];onset(find(onset(:,2)>frame.us_start-1),:)=[];
    onset_logic_acq{kk,1}=zeros(trial.acq(1),1);
    onset_logic_acq{kk,1}(onset(:,1)-trial.acq(2)+1)=true;
    onset_logic_acq{kk,2}=[[trial.acq(1):-1:1]/trial.acq(1)]';%weight
    onset_logic_acq{kk,3}=sum(onset_logic_acq{kk,1}.*onset_logic_acq{kk,2});
    onset_logic_acq{kk,4}=sheet{kk};
    m=mean(acti,2);
    % figure,plot(mean(acti,2));hold on;scatter(event.ind,m(event.ind))
    %figure,plot(acti(:,3))
    for ii=1:size(num,1)%%CR
        acti(:,ii)=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
        aa=reshape(acti(:,ii)',frame.per_cycle,trial.total)';
        CS_max(:,ii)=max(aa(:,frame.cs_start:frame.cs_end),[],2);
        CS_max(trial.acq(2):trial.acq(3),ii)=max(aa(trial.acq(2):trial.acq(3),frame.cs_start:frame.us_start-1),[],2);
        
        acti_durcs(:,ii)=reshape(aa(trial.acq(2):trial.acq(3),frame.cs_start:frame.us_start-1)',[],1);%trial.hab+1:trial.hab+trial.acq,
        if ~isempty(find(CS_max(:,ii)>1))
            CS_max(:,ii)=ones(1,trial.total)*100;
        end
        %figure,plot(CS_max(:,ii))
    end
    %CS_max_hab_test=mean(CS_max([trial.hab(2):trial.hab(3) trial.test(2):trial.test(3)],ii),2);
    %CS_max_acq=mean(CS_max(trial.acq(2):trial.acq(3),ii),2);
    %CS_max(:,3)=[];CS_max(:,16)=[];CS_max(:,16)=[];
    CS_max(:,find(CS_max(1,:)==100))=[];
    %figure,plot(CS_max)
    statistic_CR_au_acq=mean(CS_max(trial.acq(2):trial.acq(3),:),2);
    statistic_CS_max{kk,1}=mean(CS_max,2);
    statistic_CS_max{kk,2}=acti;
    statistic_CS_max{kk,3}=statistic_CR_au_acq;
    statistic_CS_max{kk,4}=[[1:trial.acq(1)]/trial.acq(1)]';%weight
    statistic_CS_max{kk,5}=sum(statistic_CS_max{kk,3}.*statistic_CS_max{kk,4});
    statistic_CS_max{kk,6}=mean(acti_durcs,2);
    statistic_CS_max{kk,7}=sheet{kk};
    statistic_CS_max{kk,8}=mean(CS_max(trial.hab(2):trial.hab(3),:),2);
    statistic_CS_max{kk,9}=mean(CS_max(trial.test(2):trial.test(3),:),2);
%     sheetfish='1009fish1';
%     filename='G:\data\summary_CR';
%     zz=2;
%     xlswrite(filename,mean(CS_max(trial.hab(1):trial.hab(2),:),2)',sheetfish,['A' num2str(zz)])
    
    
    US_max=[];
    for ii=1:size(num,1)%UR
        acti=activities_preCS_dfdf_aftcorrect{num(ii,1),1}(:,num(ii,2));
        aa=reshape(acti,frame.per_cycle,trial.total)';
        US_max(:,ii)=max(aa(trial.acq(2):trial.acq(3),frame.us_start:frame.cs_end+12),[],2);
    end
    statistic_UR_au_acq=mean(US_max,2);
    statistic_US_max{kk,1}=US_max;
    statistic_US_max{kk,2}=statistic_UR_au_acq;
    statistic_US_max{kk,3}=[[1:trial.acq(1)]/trial.acq(1)]';%weight
    statistic_US_max{kk,4}=sum(statistic_US_max{kk,2}.*statistic_US_max{kk,3});
    statistic_US_max{kk,5}=sheet{kk};
end

save([outputpath '\statistic'],'statistic_CS_max','onset_logic_acq','statistic_US_max');