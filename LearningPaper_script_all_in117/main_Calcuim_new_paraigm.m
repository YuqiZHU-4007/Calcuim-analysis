clc;clear all;close all;
set(0,'defaultfigurecolor','w');
global fs
global time
global frame
global frameb
global trial
global re_startpoint
global startpoint
global y_3sd
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
load([savepath '\Path']);

for batchi=1:4
    for fishi=1:length(Path{batchi})
        actpath=Path{batchi}{fishi};
        if ~isfile(fullfile(actpath,'activities_aft_process.mat'))
            behaviorname='behavior_from_Results_of_alltheta.mat';
            outpath=actpath
            try
                load([actpath,'act.mat'],'activities_cutback2'); load([actpath,'env.mat']);
                is_omission=1;is_opt=1;
                [~,~,frame,~,trial,Time,omission_trial,~,~,~,ishuc]=setpara_spon([],is_omission,is_opt);
%                 activities=activities_cutback2;
                activities_r=[activities(:,1),activities,activities(:,end)];
                if mod(size(activities_r,2),frame.per_cycle)~=0 %gai
                    warning(['wrong in addframe','---------act length:',num2str(size(activities_r,2))]);
                end
                if size(activities_r,2)>trial.total*frame.per_cycle
                    is_addiomission=1;
                    frame_omission=(trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+1:...
                        (trial.hab(3)+trial.acq_block_num*trial.acq_block_trial+(trial.acq_block_num-1)*(trial.acq_block_interval))*frame.per_cycle+frame.per_cycle;
                    activities_c=activities_r(:,setdiff([1:size(activities_r,2)],frame_omission));
                else
                    is_addiomission=0;activities_c=activities_r;
                end
                if size(activities_c,2)~=trial.total*frame.per_cycle
                    warning(['wrong in addomisssion','---------act length:',num2str(size(activities_c,2))]);
                end
                
                activities_spon=zeros(size(activities_r,1),size(Time.T_spon_ca{1},2),trial.acq_block_num+2);
                for ii=1:trial.acq_block_num+2
                    activities_spon(:,:,ii)=activities_c(:,Time.T_spon_ca{ii});
                end
                save([actpath '/act_spon.mat'],'activities_spon','-v7.3');
                clear('activities');clear('activities_r');clear('activities_spon');clear('activities_cutback1');clear('activities_cutback2');
                activities_c=activities_c(:,Time.T_non_spon_ca);
                [fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara_cut_spon([actpath,behaviorname]);
                save([actpath '\para'],'fs','time','frame','frameb','trial','re_startpoint','startpoint','y_3sd','ishuc','is_addiomission','-v7.3');
                
%                 cut_move=true; % ¶¼ÊÇhuc
%                 ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
%                 ind_cut_trial=unique(ind_cut_trial);
                if ishuc
                    activities_new=activities_c';
                    [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf_huc(activities_new);
                    area=[];rawacti=[];
                    outputpath=checkpath([outpath '/' 'neuron_increased intergrate dfdf']);
                    area_mean=[];
                else
                    outputpath=checkpath([outpath,'\subregoin_fromROI']);
                    [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf(activities_new);
                    [rawacti,area]=calculate_integtate_dfdf_main_control(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
                    area_mean=event_plot_new_paradigm_control([],activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,cut_move,ind_cut_trial,outputpath,'off');
                    if ~isempty([areapath,areaname])
                        fields = fieldnames(area_mean)
                        fishname=actpath;%25:38(32:46)
                        fishname=[fishname(1:8) fishname(10:end)];
                        write_A_to_xlsfile(area_mean,[areapath,areaname],XLrange,fishname);
                    end
                end
                save([actpath '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf','area','area_mean','rawacti','-v7.3');
                date = behaviorname(1:8);
                fishname = behaviorname(10:14);
                [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint,fishname,date,outputpath,1);
                save([actpath '\behavior statitic'],'p_value','bef_cond','dur_cond','aft_cond','-v7.3');
            catch
                warning(['wrong in ',actpath])
            end
        end
    end
end
