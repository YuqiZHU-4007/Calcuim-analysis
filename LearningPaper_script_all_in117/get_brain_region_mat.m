clc;
clear all;
%% 正常traning

load('X:\calcium data 20230224\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');
temp_supervoxel=temp_env.env.supervoxel(:,1:3);
%temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
warped_SyN_csv_path='X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\';
% 包含下面这些
% 20190514fish3vol_env_spatialloc_warped_SyN.csv
% 20190514fish3vol_env_spatialloc_warped_SyN_add_brainregion.csv

list=dir('X:\calcium data 20230224\DS_MV_TO_DS_TEMP_adjust_location\*.csv');
savepath='X:\calcium data 20230224\';
for fishi=1:length(list)
    p=list(fishi).folder;
    nn=list(fishi).name;
    %load(fullfile(p,nn));
    nn=[nn(1:13)];
    if ~isempty(strfind(nn,'2019'))
        res=[0.66,0.66,8];
    else
        res=[0.66,0.66,10.0465];
    end
    if exist([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv'])~=0
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        supervolxeli(:,1)= supervolxeli(:,1)/res(1); supervolxeli(:,2)= supervolxeli(:,2)/res(2);supervolxeli(:,3)=supervolxeli(:,3)/res(3)+1;
        ind_supervoxel=1:size(supervolxeli,1);
%         figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),5,[0.5 0.5 0.5],'filled');hold on;
%         scatter3(env.supervoxel( ind_supervoxel,1),env.supervoxel( ind_supervoxel,2),env.supervoxel( ind_supervoxel,3),10,[1 0 0],'filled');hold on;
        if ~isempty(ind_supervoxel)
            gIX=ones(1,length(ind_supervoxel))';supervoxel_in_fish=supervolxeli(ind_supervoxel,1:3);clrmap = GetColormap('hsv_new',max(gIX));
            for fraction_type=1
                switch fraction_type
                    case 1
                        temp_vol=supervolxeli(:,1:3);
                    case 2
                        temp_vol=temp_supervoxel;
                end
                [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,brain_region_id,loc_in_region_temp,id_in_region_temp,num_in_region_in_clust,Label]=get_region_fraction_temp_preprocess(gIX,ind_supervoxel,supervoxel_in_fish,reg_mask,reg_name,temp_vol,clrmap);
                %[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,ind_supervoxel,supervoxel_in_fish,reg_mask,reg_name,temp_vol,clrmap);
                pp=checkpath(fullfile(savepath,nn(1:8),nn(9:13)));
                if exist([pp,'\brain_region_related_statistic.mat'])~=0
                    delete([pp,'\brain_region_related_statistic.mat'])
                end
                save([pp,'\brain_region_related_statistic.mat'],'loc_in_region_in_clust','index_in_region_in_clust','fraction_in_region_in_clust','brain_region_id','loc_in_region_temp','id_in_region_temp','num_in_region_in_clust','Label');
                columns = {'x', 'y', 'z', 't','index','brain_region','label'};
                M=table(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),supervolxeli(:,4),supervolxeli(:,5),brain_region_id(:,1),brain_region_id(:,2),'VariableNames', columns);
                outpath=[p,'\' nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv'];
                writetable(M,outpath);
            end
        end
    else
        warning([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
    end
end

%% Hb activation
clc;clear all
load('X:\calcium data 20230224\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('X:\calcium data 20230224\脑区分割\env.mat');
temp_supervoxel=temp_env.env.supervoxel(:,1:3);
warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\';
list=dir('X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\DS_MV_TO_DS_TEMP_adjust_location\*.csv');
% warped_SyN_csv_path='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\test\';
% list=dir('X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\test\*.csv');

savepath='X:\data\fear conditioning_ZYQ\Hb_activation\Hb_after_processing\';
for fishi=1:length(list)
    p=list(fishi).folder;
    nn=list(fishi).name;
    %load(fullfile(p,nn));
    nn=[nn(1:13)];
    if ~isempty(strfind(nn,'2019'))
        res=[0.66,0.66,8];
    else
        res=[0.66,0.66,10.0465];
    end
    if exist([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv'])~=0
        supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
        supervolxeli(:,1)= supervolxeli(:,1)/res(1); supervolxeli(:,2)= supervolxeli(:,2)/res(2);supervolxeli(:,3)=supervolxeli(:,3)/res(3)+1;
        ind_supervoxel=1:size(supervolxeli,1);
%         figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),5,[0.5 0.5 0.5],'filled');hold on;
%         scatter3(env.supervoxel( ind_supervoxel,1),env.supervoxel( ind_supervoxel,2),env.supervoxel( ind_supervoxel,3),10,[1 0 0],'filled');hold on;
        if ~isempty(ind_supervoxel)
            gIX=ones(1,length(ind_supervoxel))';supervoxel_in_fish=supervolxeli(ind_supervoxel,1:3);clrmap = GetColormap('hsv_new',max(gIX));
            for fraction_type=1
                switch fraction_type
                    case 1
                        temp_vol=supervolxeli(:,1:3);
                    case 2
                        temp_vol=temp_supervoxel;
                end
                [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,brain_region_id,loc_in_region_temp,id_in_region_temp,num_in_region_in_clust,Label]=get_region_fraction_temp_preprocess(gIX,ind_supervoxel,supervoxel_in_fish,reg_mask,reg_name,temp_vol,clrmap);
                %[loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,ind_supervoxel,supervoxel_in_fish,reg_mask,reg_name,temp_vol,clrmap);
                pp=checkpath(fullfile(savepath,nn(1:8),nn(9:13)));
                if exist([pp,'\brain_region_related_statistic.mat'])~=0
                    delete([pp,'\brain_region_related_statistic.mat'])
                end
                save([pp,'\brain_region_related_statistic.mat'],'loc_in_region_in_clust','index_in_region_in_clust','fraction_in_region_in_clust','brain_region_id','loc_in_region_temp','id_in_region_temp','num_in_region_in_clust','Label');
                columns = {'x', 'y', 'z', 't','index','brain_region','label'};
                M=table(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),supervolxeli(:,4),supervolxeli(:,5),brain_region_id(:,1),brain_region_id(:,2),'VariableNames', columns);
                outpath=[p,'\' nn,'vol_env_spatialloc_warped_SyN_add_brainregion.csv'];
                writetable(M,outpath);
            end
        end
    else
        warning([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
    end
end


%% check
clr=jet(31);a=str2double(brain_region_id(:,1));b=isnan(a);a(b)=31;
figure,scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),5,[0.5 0.5 0.5],'filled');hold on;
for ii=1:length(index_in_region_in_clust);
    ind=index_in_region_in_clust{ii};
    scatter3(supervolxeli(ind,1),supervolxeli(ind,2),supervolxeli(ind,3),10,clr(a(ind),:),'filled');hold on;
end
