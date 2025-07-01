function get_brain_region_related_statistic(p)
res=[0.66,0.66,10];
load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';



%p='H:\1.Test US\2.Tail free——Data from 117\20220709\fish1\';
load(fullfile(p,'para.mat'));
load(fullfile(p,'env.mat'));
nn=[p(end-14:end-7),p(end-5:end-1)];
supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
ind=1:size(supervolxeli,1);
if ~isempty(ind)
    gIX=ones(1,length(ind))';supervoxel_in_fish=supervolxeli(ind,1:3);clrmap = GetColormap('hsv_new',max(gIX));
    for fraction_type=2
        switch fraction_type
            case 1
                temp_vol=supervolxeli(:,1:3);
            case 2
                temp_vol=temp_supervoxel;
        end
        [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervoxel_in_fish,reg_mask,reg_name,reg_loc,temp_vol,temp_env,clrmap,nn);
        save(fullfile(p,'brain_region_related_statistic.mat'),'Label','brain_region_id','fraction_in_region_in_clust','loc_in_region_cell','loc_in_region_in_clust','num_in_region_in_clust');
    end
end

%unique(strcmp(brain_region_id2(:,1),brain_region_id(:,1)))