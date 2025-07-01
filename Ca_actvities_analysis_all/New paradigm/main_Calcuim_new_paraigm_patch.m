%batch
%20190328
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

%batch code -----------------------------
actpathes=batch_code('txt of <activties_new> path') ;
loc_pathes=batch_code('txt of <location> path') ;
beh_pathes=batch_code('txt of <behavior> path') ;
save_pathes=batch_code('txt of <savepath> path') ;
para_pathes=batch_code('txt of <para> path') ;
[areaname,areapath]=uigetfile(['E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328\*.xls'],'area_lowest_all');
%areapath='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all learner-shake response-cutmove-20190328\control-CS\area_lowest_all';

nbatch = length(actpathes);
if isfile([areapath '\ind_cut_trial_total.mat'])
    load([areapath '\ind_cut_trial_total.mat']);
else
    ind_cut_trial_total={};
end
for batchi=[6,9];%1:nbatch;
    actpath= actpathes{batchi};
    display(actpath);
    locpath=loc_pathes{batchi};
    behaviorpath=beh_pathes{batchi};
    outpath=save_pathes{batchi};
    load(actpath);
    %set para.
    %[fs,time,frame,frameb,trial,re_startpoint,startpoint,y_3sd,ishuc]=setpara(behaviorpath);
    parapath=para_pathes{batchi};
    load(parapath);
    outputpath=checkpath([outpath,'\subregoin_fromROI']);
    
    cut_move=false;
    ind_cut_trial=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
    ind_cut_trial_total{batchi}=unique(ind_cut_trial);
    if ishuc
        load(locpath);
        activities_new=activities';
        [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf_huc(activities_new);
        rawacti=[];area=[];
        %[rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
        %     outputpath=checkpath([outpath,'\subregoin_supervol']);
        %     event_plot_new_paradigm_huc([locpath,loctname],activities_preCS_dfdf_aftcorrect,activities_event_preCS,outputpath);
        %outputpath=checkpath([outpath '/' 'neuron_increased intergrate dfdf']);
        area_mean=[];
        save([actpath(1:end-8) '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf','area','area_mean','rawacti','-v7.3');
        activities=activities';
        activities=[activities(1,:);activities;activities(size(activities,1),:);];
        n=randi(size(activities,1),[1,8]);
        figure,plot_rawtrace_trials(activities_preCS_dfdf_aftcorrect(:,n),[],fs,frame,trial,[],1);title('aft process')
        figure,plot_rawtrace_trials(activities(:,n),[],fs,frame,trial,[],1);title('raw trace')
    else
        [activities_event_preCS,activities_preCS_dfdf_aftcorrect,activities_preCS_dfdf]=Ca_analyse_for_new_paradigm_computedfdf(activities_new);
        [rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
        %[rawacti,area]=calculate_integtate_dfdf_main_hab5(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
        
        %plot
        % [actdfname,actdfpath]=uigetfile('G:\data\.mat','activities_aft_process');
        % load([actdfpath,actdfname]);
        area_mean=event_plot_new_paradigm([locpath],activities_preCS_dfdf_aftcorrect,activities_event_preCS,area,cut_move,ind_cut_trial,outputpath,'off');
        save([actpath(1:end-18) '\activities_aft_process'],'activities_event_preCS','activities_preCS_dfdf_aftcorrect','activities_preCS_dfdf','area','area_mean','rawacti','-v7.3');
        fields = fieldnames(area_mean);
        %XLrange=char(65+(batchi-1)*length(fields));
        XLrange={'A','K','U','AE','AO','AY','BI'};
        fishname=actpath(32:46);%25:38
        fishname=[fishname(1:8) fishname(10:end)];
        %     if batchi==1
        %         delete(areapath)
        %     end
        write_A_to_xlsfile(area_mean,[areapath,areaname],XLrange{batchi},fishname);
        %write_A_to_xlsfile(area_mean,[areapath,areaname],'AO',fishname);
    end
    
    %behavior onset
     [fpath,fname,ext]=fileparts(behaviorpath);
    outputpath=fpath;
    a=load(behaviorpath);
    h=figure;
    title(fname,'Interpreter','none','fontsize',13);
    plot_behavior_onset(a.delta_r_bef1,y_3sd,fs,frame,frameb,trial,re_startpoint);
    saveas(h,[outputpath '\'  'behavior onset' '.tif']);
    savefig(h,[outputpath '\'  'behavior onset']);
    close(h)
    
    date = fname(1:8);
    fishname = fname(10:14);
    [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(a.delta_r_bef1,re_startpoint,fishname,date,outputpath);
    save([outputpath '\behavior statitic'],'p_value','bef_cond','dur_cond','aft_cond','-v7.3');
end
save([areapath '\ind_cut_trial_total'],'ind_cut_trial_total');