clc;clear all;close all
res=[0.66,0.66,10];
%% identify regressor
%1: CS processing(high/no corr to CR performance)
%2: motion
%3: teaching
is_addiomission=1;
p='H:\1.Test US\2.Tail free¡ª¡ªData from 117\20230309\fish1\';
load(fullfile(p,'CR_ind_summary.mat'),'CS_up_ind');
load(fullfile(p,'activities_aft_process.mat'));
load(fullfile(p,'brain_region_related_statistic.mat'));
load(fullfile(p,'env.mat'));
load(fullfile(p,'act.mat'),'activities');
load(fullfile(p,'behavior\behavior_from_Results_of_alltheta.mat'));

[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara([]);
trial.acq_block_interval=10;
trial.test(2)=trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num)*trial.acq_block_interval+1;
trial.test(3)=trial.test(2)+trial.test(1)-1;
trial.spon_aft=[trial.spon_aft(1) ,trial.test(3)+1,trial.test(3)+trial.spon_aft(1)];
frame.us_dur=0.001*ones(1,trial.acq_block_num)*fs.behavior;
trial.total=trial.total+trial.acq_block_interval*(trial.acq_block_num);
if is_addiomission==1
    addomission=[[trial.hab(3)+8*trial.acq_block_trial+(8-1)*(trial.acq_block_interval)]*frame.per_cycle+1:(trial.hab(3)+8*trial.acq_block_trial+(8-1)*(trial.acq_block_interval)+1)*frame.per_cycle];
    act_addomision=activities(:,addomission)';
    act_addomision_dfdf=(act_addomision-repmat(mean(act_addomision,1),size(act_addomision,1),1))./repmat(mean(act_addomision,1),size(act_addomision,1),1);
else
    act_addomision=[];act_addomision_dfdf=[];
end
figure,imagesc(act_addomision_dfdf',[-0.00 0.01])
T_non_spon=[];
for ii=1:trial.acq_block_num+2
    switch ii
        case 1
            t=(trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3))*frame.per_cycle;
        case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
            t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
                (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle ;
        case trial.acq_block_num+2
            if is_addiomission==1
                t=[(trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle]+frame.per_cycle;
            else
            t=(trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle;
            end
    end
     T_non_spon=[T_non_spon,t];
end
activities_raw=activities(:,T_non_spon)';
load(fullfile(p,'para.mat'));


iscutmov=1;
ind=CS_up_ind{2,iscutmov};
for ii=1:size(CS_up_ind,1)
     ind=intersect(CS_up_ind{ii,iscutmov},ind);
end
act_stim=activities_preCS_dfdf_aftcorrect(:,ind);
figure, plot_rawtrace_trials(mean(act_stim,2),[],fs,frame,trial,startpoint,1);
loc=env.supervoxel(ind,1:3);
loc_region=brain_region_id(ind,2);
figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),2,'k');hold on;scatter3(loc(:,1),loc(:,2),loc(:,3),16,'r','filled');
%% find neurons
regress1=zeros(5,trial.total*frame.per_cycle);regress1(1,:)=mean(act_stim,2);
% [stregressors,~,~,~] = GetMotorRegressor(regress1,2);
% stAvrcorr_all=[];
% stim_output=[stregressors(1,1).im];
% figure,plot(regress1(1,:));hold on;
% plot(stim_output,'r');legend('raw','regressor');
[stimcorr,~] = MotorSourceCorrelation(activities_preCS_dfdf_aftcorrect',regress1(1,:),[]);
figure,hist(stimcorr);hold on;
thr=mean(stimcorr)+3*std(stimcorr);
line([thr,thr],[0 2000],'color','r','linewidth',3);
CS_related_ind=find(stimcorr>thr);
CS_related_act=activities_preCS_dfdf_aftcorrect(:,CS_related_ind);
figure,plot(regress1(1,:));hold on;
plot(mean(CS_related_act,2),'r');legend('regressor','corr > thr');
figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),2,'k');hold on;
scatter3(env.supervoxel(CS_related_ind,1),env.supervoxel(CS_related_ind,2),env.supervoxel(CS_related_ind,3),16,'r','filled');
CS_related_locregion=brain_region_id(CS_related_ind,2); unique(CS_related_locregion)
figure, plot_rawtrace_trials(mean(CS_related_act,2),[],fs,frame,trial,startpoint,1);
figure, plot_rawtrace_trials(mean(activities_raw(:,CS_related_ind),2),[],fs,frame,trial,startpoint,1);

figure,plot(mean(act_addomision_dfdf(:,CS_related_ind),2));
figure,imagesc(act_addomision_dfdf(:,CS_related_ind)',[-0.00,0.01])
%% cut US decrease
A_r=reshape(activities_preCS_dfdf_aftcorrect,frame.per_cycle,trial.total,[]);
cr=squeeze(mean(A_r(frame.cs_start:frame.us_start-1,trial.acq(2):trial.acq(3),:),1));
ur=squeeze(mean(A_r(frame.us_start:frame.us_start+1,trial.acq(2):trial.acq(3),:),1));
ur=cr-ur;
ur_df=ur(2:end,:)-ur(1:end-1,:);
ur_df_norm=ur_df;ur_df_norm(find(ur_df>0))=0;ur_df_norm(find(ur_df<0))=-1;
a=sum(ur_df_norm,1);
US_related_ind=find(sum(ur_df_norm,1)<-15);
figure, plot_rawtrace_trials(mean(activities_preCS_dfdf_aftcorrect(:,US_related_ind),2),[],fs,frame,trial,startpoint,1);
figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),2,'k');hold on;
scatter3(env.supervoxel(US_related_ind,1),env.supervoxel(US_related_ind,2),env.supervoxel(US_related_ind,3),16,'r','filled');

%% statistic
