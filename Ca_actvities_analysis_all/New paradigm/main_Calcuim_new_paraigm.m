clc;clear all;close all;
global fs
global time
global frame
global frameb
global trial
global re_startpoint
global startpoint
global y_3sd

set(0,'defaultfigurecolor','w');
[actname,actpath]=uigetfile('J:\newfile\.mat','activities new');
[loctname,locpath]=uigetfile([actpath '.xls'],'location');
[behaviorname,behaviorpath]=uigetfile([actpath '.mat'],'behavior');
[envname,envpath]=uigetfile([actpath '.mat'],'env');
outpath=uigetdir(actpath,'outputpath');
% [areaname,areapath]=uigetfile(['E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\*.xls'],'area_event');
% XLrange=[]
% while ~ischar(XLrange)
%     display('XLrange is not a char');
%     XLrange=input('input the col name to write "area"(must be char) \n ') %{'A','K','U','AE','AO','AY','BI'};
% end

[~,~,frame,~,trial,~,~,T_non_spon,T_spon]=setpara_spon();
load([actpath,actname]); load([envpath,envname]);
activities=activities_cutback2;
activities_r=[activities(:,1),activities,activities(:,end)];
activities_spon=zeros(size(activities_r,1),size(T_spon{1},2),trial.acq_block_num+2);
if size(activities_r,2)~=trial.total*frame.per_cycle
    warning(['wrong in addomisssion','---------act length:',num2str(size(activities_r,2))]);
end
if activities_r>trial.total*frame.per_cycle
    is_addiomission=1;
    frame_omission=(trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+1:...
        (trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+frame.per_cycle;
    activities_c=activities_r(:,setdiff([1:size(activities_r,2)],frame_omission));
else
    is_addiomission=0;activities_c=activities_r;
end

for ii=1:trial.acq_block_num+2
    activities_spon(:,:,ii)=activities_c(:,T_spon{ii});
end
save([actpath '/act_spon.mat'],'activities_spon','-v7.3');

activities_c=activities_c(:,T_non_spon);
[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara_cut_spon([behaviorpath,behaviorname]);
%activities_c(:,1:1500)=[];activities_c=activities_c(:,1:2100);
%figure, plot_rawtrace_trials(mean(activities_c(:,:),1)',[],fs,frame,trial,[],1);
save([actpath '\para'],'fs','time','frame','frameb','trial','re_startpoint','startpoint','y_3sd','ishuc','is_addiomission','-v7.3');

cut_move=true;
ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
ind_cut_trial=unique(ind_cut_trial)
if ishuc
    activities_new=activities_c';
    [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf_huc(activities_new);
    area=[];rawacti=[];
    %[rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
    %     outputpath=checkpath([outpath,'\subregoin_supervol']);
    %     event_plot_new_paradigm_huc([locpath,loctname],activities_preCS_dfdf_aftcorrect,activities_event_preCS,outputpath);
    outputpath=checkpath([outpath '/' 'neuron_increased intergrate dfdf']);
    area_mean=[];
    %[p,q,ind_increased]=increased_intergrate_dfdf(activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,env,outputpath);
else
    outputpath=checkpath([outpath,'\subregoin_fromROI']);
    [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf(activities_new);
    %[rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
    %!!!!!!!!!!!!!!!!!!!!
    [rawacti,area]=calculate_integtate_dfdf_main_control(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
    % [actdfname,actdfpath]=uigetfile('G:\data\.mat','activities_aft_process');
    % load([actdfpath,actdfname]);
    %area_mean=event_plot_new_paradigm([locpath,loctname],activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,cut_move,ind_cut_trial,outputpath,'off');
    area_mean=event_plot_new_paradigm_control([locpath,loctname],activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,cut_move,ind_cut_trial,outputpath,'off');
    if ~isempty([areapath,areaname])
    fields = fieldnames(area_mean)
    fishname=actpath;%25:38(32:46)
    fishname=[fishname(1:8) fishname(10:end)];
    write_A_to_xlsfile(area_mean,[areapath,areaname],XLrange,fishname);
    end
end

save([actpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf','area','area_mean','rawacti','-v7.3');


% % % %behavior onset
% % % outputpath=behaviorpath;
% % % a=load([behaviorpath,behaviorname]);
% % % % h=figure;
% % % % title(behaviorname,'Interpreter','none','fontsize',13);
% % % % plot_behavior_onset(a.delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
% % % % %plot_behavior_onset(delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
% % % % saveas(h,[outputpath '\'  'behavior onset' '.tif']);
% % % % savefig(h,[outputpath '\'  'behavior onset']);
% % % % close(h)
% % % 
date = behaviorname(1:8);
fishname = behaviorname(10:14);
[p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint,fishname,date,outputpath,1);
save([actpath '\behavior statitic'],'p_value','bef_cond','dur_cond','aft_cond','-v7.3');
% hh=figure;
% plot_behavior_statistics(fs,frame,frameb,trial,re_startpoint)
% saveas(hh,[outputpath '\'  'behavior statistics' '.tif']);
% savefig(hh,[outputpath '\'  'behavior statistics']);
%close(hh);
% % % [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint,[],[],[],1);

%% cut spontanouse
% %spon
% T_spon={};
% for ii=1:trial.acq_block_num+2
%     switch ii
%         case 1
%             t=(trial.spon_bef(2)-1)*frame.per_cycle+1:(trial.spon_bef(3))*frame.per_cycle;
%         case mat2cell([2:trial.acq_block_num],1,ones(1,trial.acq_block_num-1))
%             t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
%                 (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
%         case trial.acq_block_num+1
%             t=(trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
%                 (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-1)*(trial.acq_block_interval))*frame.per_cycle ;
%             if is_addiomission t=t+frame.per_cycle;end
%         case trial.acq_block_num+2
%             t=(trial.test(3))*frame.per_cycle+1:(trial.total)*frame.per_cycle;
%             if is_addiomission t=t+frame.per_cycle;end
%     end
%     T_spon{ii}=t;
% end
% %spon
% T_non_spon=[];
% for ii=1:trial.acq_block_num+2
%     switch ii
%         case 1
%             t=(trial.hab(2)-1)*frame.per_cycle+1:(trial.hab(3))*frame.per_cycle;
%         case mat2cell([2:trial.acq_block_num+1],1,ones(1,trial.acq_block_num))
%             t=(trial.hab(3)+(ii-2)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle+1 : ...
%                 (trial.hab(3)+(ii-1)*trial.acq_block_trial+(ii-2)*(trial.acq_block_interval))*frame.per_cycle ;
%         case trial.acq_block_num+2
%             if is_addiomission==1
%                 t=[(trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle]+frame.per_cycle;
%             else
%             t=(trial.test(2)-1)*frame.per_cycle+1:(trial.test(3))*frame.per_cycle;
%             end
%     end
%      T_non_spon=[T_non_spon,t];
% end
