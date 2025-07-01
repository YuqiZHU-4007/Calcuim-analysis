clc;clear all;close all;
set(0,'defaultfigurecolor','w');
file_path='X:\calcium data 20230224\';
num{1}=[1,3,4,6,8,7,9,10];%[1,2,3,4,6,7,8,9,10];[1,3,4,6,8,9];
num{2}=[2,3,4,6,7,8,10];
num{3}=[1:9];
num{4}=[1,2,4,6];

mat_name=[];
Path{1}={[file_path '\20220709\fish1\',mat_name],...%差，no Cb after Learning
    [file_path '\20220802\fish2\',mat_name],...%差，不用
    [file_path '\20220803\fish2\',mat_name],...%好
    [file_path '\20220814\fish1\',mat_name],...%一般，结果不好
    [file_path '\20220915\fish1\',mat_name],...%calcium wrong，不用
    [file_path '\20230309\fish1\',mat_name],...%好
    [file_path '\20230807\fish1\',mat_name],...%%一般，acq后半段calcium奇怪
    [file_path '\20210402\fish1\',mat_name],...%好，Hb配准有点差
    [file_path '\20210709\fish2\',mat_name],...%好
    [file_path '\20210805\fish1\',mat_name],...%配准差
    };%learner
Path{2}={
    [file_path '\20210514\fish2\',mat_name],...%不用，un shuffle
    [file_path '\20210702\fish1\',mat_name],...
    [file_path '\20220831\fish1\',mat_name ],...
    [file_path '\20221108\fish1\',mat_name ],...
    [file_path '\20230309\fish2\',mat_name ],...%不用，行为像learner
    [file_path '20231110\fish1\',mat_name ],... 
    [file_path '20240320\fish1\',mat_name ],... 
    [file_path '\20240320\fish2\',mat_name ],... 
    [file_path '\20240408\fish1\',mat_name ],... %
    [file_path '\20240416\fish2\',mat_name ] %
    %     ['E:\A_Data_lightsheet\Data_huc\20190612\fish1\',mat_name ],...
    %     ['E:\A_Data_lightsheet\Data_huc\20190612\fish2\',mat_name ]
    %['H:\1.Test US\2.Tail free——Data from 117\20210617\fish2\',mat_name],...
    };%unpair
Path{3}={[file_path '\20220803\fish1\',mat_name ],...
    [file_path '\20220816\fish2\',mat_name ],...
    [file_path '\20220824\fish1\',mat_name ],...
    [file_path '\20220829\fish1\',mat_name ],...
    [file_path '\20221109\fish1\',mat_name],...
    [file_path '\20230314\fish2\',mat_name],...
    [file_path '\20230807\fish2\',mat_name],...
    [file_path '\20210510\fish1\',mat_name ],...
    [file_path '20210529\fish1\',mat_name]
    %     ['G:\data_huc_non_learner\20190604\fish2\',mat_name ],...
    %     ['G:\data_huc_non_learner\20190607\fish3\',mat_name ]
    };%non-learner
Path{4}={[file_path '\20220817\fish1\',mat_name ],...
    [file_path '\20210401\fish1\',mat_name],...
    [file_path '\20210717\fish1\',mat_name],...%配准差
    [file_path '\20210822\fish1\',mat_name],...
    [file_path '\20230308\fish1\',mat_name],...%配准不好
    [file_path '\20230808\fish1\',mat_name],...
    %       ['H:\1.Test US\2.Tail free——Data from 117\20210719\fish1\',mat_name],...
    };%faded-learner
savepath='X:\calcium data 20230224\';
save([savepath '\Path'],'Path','num','-v7.3');

for batchi=3:4
    for ii=1:length(Path{batchi})
        a=['X:\calcium data 20230224\',Path{batchi}{ii}(end-14:end)];
        Path{batchi}{ii}=a;
    end
end

clc;clear all;
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
% delete(gcp('nocreate'))
% addpath(genpath('X:\calcium data 20230224\script\'))
for ii=3:4
    path=Path{ii};
    %parpool(2);
    for jj=1:length(path)
        US_mapping_singletrial_onlyUS(path{jj})
    end
end

%% preprocess of Hb activation
clc;clear all;close all;
set(0,'defaultfigurecolor','w');
file_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing';
num{1}=[1:3];num{3}=[1:3];num{4}=[1:3];
mat_name=[];
Path_Hb{1}={[file_path '\20240125\fish2\',mat_name],...
    [file_path '\20240126\fish1\',mat_name],...%
    [file_path '\20240419\fish2\',mat_name],...%
    };%learner
Path_Hb{2}={[file_path '\20240110\fish1\',mat_name ],...a
    [file_path '\20240118\fish1\',mat_name],...%
    [file_path '\20240125\fish1\',mat_name],...%
    };%non-learner
Path_Hb{3}={[file_path '\20240117\fish1\',mat_name ],...
    [file_path '\20240118\fish2\',mat_name],...
    [file_path '\20240522\fish1\',mat_name],...
    };%faded-learner
savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
save([savepath '\Path_Hb'],'Path_Hb','num','-v7.3');


clc;clear all;
savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path']);
% delete(gcp('nocreate'))
% addpath(genpath('X:\calcium data 20230224\script\'))
for ii=[1,3,4]
    path=Path{ii};
    %parpool(2);
    for jj=1:length(path)
        %US_mapping_singletrial_hbactivation(path{jj})
        US_mapping_singletrial_fraction(path{jj})
    end
end
%% 去自发
% clc;clear all;
% Region_name='L H';isslidewin=1;
% savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
% load([savepath '\Path']);
% savetype='jpg';
% for batchi=1:3
%     for fishi=num{batchi}
%         path=Path{batchi}{fishi};%'H:\1.Test US\5.fear conditioning behavioral data\20220709\fish1\';
%         load(fullfile(path, '/activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect');
%         load(fullfile(path,'/para.mat'));
%         load(fullfile(path,'/behavior_from_Results_of_alltheta.mat'));
%         %load(fullfile(path,'singletrial','CR_ind_summary_singletrial.mat'),'ind_all_CSUS_RESPONSIVE');
%         [p_value,bef_cond,dur_cond,aft_cond]=tail_movement_by_trial_fin_cor(y_3sd,re_startpoint_sd,[],[],[],1);
% 
%         activities_preCS_dfdf_aftcorrect=activities_preCS_dfdf_aftcorrect(Time.T_non_spon_ca,:);
%         a=mean(activities_preCS_dfdf_aftcorrect,2,'omitnan');
%         h=figure;plot_rawtrace_trials(a,[],fs,frame,trial,[],1);
%         save(fullfile(path, '/activities_aft_process.mat'),'activities_preCS_dfdf_aftcorrect','-v7.3');
%     end
% end
% 

