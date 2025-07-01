load('G:\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_nonlearner_cutmove_20190921\path\path.mat');
name='activities_aft_process';
for tt=1:size(actbatch,1)
    [fpath,fname,ext]=fileparts(actbatch{tt});
    if ~strcmp(fname,name)
        fname=name;
        actbatch{tt}=fullfile(fpath,[fname,ext]);
    end
end

K=[2,3,4,6,9];id=[];
for kk=K
    load(actbatch{kk});
    load(parabatch{kk}); load(envbatch{kk});
    
    act=activities_preCS_dfdf_aftcorrect;
    act=reshape(act,frame.per_cycle,trial.total,[]);
    load([savebatch{kk} '\brain arearegion_mask.mat']);
    [loc_in_region_in_clust,~,~,id_in_region_cell]=get_region_fraction(reg_mask,reg_name,reg_loc,ones(length(env.supervoxel),1),[1:length(env.supervoxel)],[env.supervoxel(:,2),env.supervoxel(:,1),env.supervoxel(:,3)],[1 0 0],false);
    figure,scatter(env.supervoxel(:,2),env.supervoxel(:,1));hold on ; scatter(env.supervoxel(id_in_region_cell{16},2),env.supervoxel(id_in_region_cell{16},1))
    acti_per_region=struct;
    for rr=1:size(reg_name,2)
        a=mean(act(:,:,id_in_region_cell{rr}),3);a=a(:,1:42);
        acti_per_region=setfield(acti_per_region,strrep(reg_name{rr},' ','_'),a);
    end
    save([savebatch{kk} '\acti_per_region' strrep(savebatch{kk}(end-13:end),'\','_') '.mat'],'acti_per_region');
    figure,plot(acti_per_region.R_Pallium_middle)
    loc_per_region=struct;
    name=fieldnames(reg_loc);figure,
    for rr=1:size(name,1)
        a=getfield(reg_loc,name{rr});a=mean(a,2)';
        loc_per_region=setfield(loc_per_region,name{rr},a);
        scatter(a(:,1),a(:,2),a(:,3));hold on;
    end
    axis equal;
    save([savebatch{kk} '\loc_per_region' strrep(savebatch{kk}(end-13:end),'\','_') '.mat'],'loc_per_region');
    %         acti_per_region_spon=struct;
    %     for rr=1:size(reg_name,2)
    %         a=mean(act(:,:,id_in_region_cell{rr}),3); a=a(6:frame.cs_start-1,trial.acq(2):trial.acq(3));a=reshape(a,size(a,1)*trial.acq_block_trial,[]);
    %         figure,plot(a)
    %         acti_per_region_spon=setfield(acti_per_region_spon,strrep(reg_name{rr},' ','_'),a);
    %     end
    %     save([savebatch{kk} '\acti_per_region_spon' savebatch{kk}(end-13:end)],'acti_per_region_spon');
end