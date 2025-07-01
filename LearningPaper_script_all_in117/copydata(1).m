clc;clear all;close all
savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
load([savepath '\Path'],'Path');
group={'Learner','Unpair','Non-learner','Faded-learner'};
savepath=checkpath('X:\学习数据_DULAB_ZYQ\Hb activation');%checkpath('L:\学习数据_DULAB_ZYQ');
csvpath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
filename={'act.mat','activities_dfdf_align.mat','activities_aft_process.mat','behavior_from_Results_of_alltheta.mat','brain_region_related_statistic.mat','para.mat','env.mat'};
%filename={'act_spon_dfdf.mat'};

for batchi=1:4
    path=Path{batchi};
    checkpath(fullfile(savepath,group{batchi}));
    for ii=1:length(path)
        pp=path{ii};
        nn=[pp(end-14:end-7),pp(end-5:end-1)];
        if ~isfile(fullfile(savepath,group{batchi},[nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv']))
            %delete(fullfile(savepath,group{batchi},[nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv']))
            disp(fullfile(csvpath,[nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv']))
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(csvpath,[nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv']),fullfile(savepath,group{batchi},[nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv']));
            if SUCCESS~=1
                warning(pp)
            end
        end
        for filei=1:length(filename)
            if ~isfile(fullfile(savepath,group{batchi},nn,filename{filei}))
                disp(fullfile(pp,filename{filei}))
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(pp,filename{filei}),fullfile(checkpath(fullfile(savepath,group{batchi},nn)),filename{filei}));
                if SUCCESS~=1
                    warning(pp)
                end
            end
        end
    end
end


%%
savepath='X:\calcium data 20230224\';
load([savepath '\Path']);
for batchi=1:4
    path=Path{batchi};
    for fishi=1:length(path)
        p=path{fishi} ;
        nn1=[p(end-14:end-7)];
        nn2=[p(end-5:end-1)];
        [SUCCESS,MESSAGE,MESSAGEID] = copyfile(fullfile(p,'activities_dfdf_align.mat'),...
            fullfile(checkpath(fullfile('X:\copy',nn1,nn2)),'activities_dfdf_align.mat'));
        if SUCCESS~=1
            warning(p)
        end
        
    end
end