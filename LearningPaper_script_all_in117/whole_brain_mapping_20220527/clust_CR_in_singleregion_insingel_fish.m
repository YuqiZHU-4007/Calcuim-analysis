%% load temp
clc;clear all;close all;
res=[0.66,0.66,10];radium=floor(30);
load('H:\3.Juvenile reference brain\registration to templete\脑区分割\segmentation_file_0525_DSregion_mask.mat');
load('H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\Path.mat');
temp_env=load('H:\3.Juvenile reference brain\registration to templete\脑区分割\env.mat');
warped_SyN_csv_path='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\regis_results\after_regist\DS_MV_TO_DS_TEMP_adjust_location\';
savepath='H:\3.Juvenile reference brain\registration to templete\data\All_cal_summary\Pre_Post_CR\';
temp_supervoxel(:,1)=temp_env.env.supervoxel(:,1)*res(1);temp_supervoxel(:,2)=temp_env.env.supervoxel(:,2)*res(2);temp_supervoxel(:,3)=(temp_env.env.supervoxel(:,3)-1)*res(3);

%% clust regioni
Lable={'L OB','R OB','L P','R P',...
    'L H','R H','L TeO','R TeO','L Np','R Np','L Otr','R Otr','TL','PO',...
    'Pr','PT','Th','rH','iH','cH','mPT','T',...
    'L TS','R TS','Va','CCe','L gT','R gT','R1','R2'};
for batchi=2;length(Path);
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi};nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
        if exist(fullfile(path_fishi,'brain_region_related_statistic.mat'),'file')
            try
                load(fullfile(path_fishi,'activities_aft_process.mat'));
                load(fullfile(path_fishi,'para.mat'));
                load(fullfile(path_fishi,'env.mat'));
                load(fullfile(path_fishi,'brain_region_related_statistic.mat'));
                supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                activities_preCS_dfdf_aftcorrect_zscore=normalize(activities_preCS_dfdf_aftcorrect,1,'zscore');
                region_id=[9,10;11,12;23,24;14,0;15,0;16,0;17,0;27,28;22,0;21,0;];
                for regioni=1:size(region_id,1)
                    a=region_id(regioni,:); a(find(a==0))=[];
                    if length(a)>1
                        Region=strcat(Label{region_id(regioni,1)},'-',Label{region_id(regioni,2)});
                        neuron_id_in_regioni=find(strcmp(brain_region_id(:,2),Label{region_id(regioni,1)}) | strcmp(brain_region_id(:,2),Label{region_id(regioni,2)}));
                    else
                        Region=Label{region_id(regioni,1)};neuron_id_in_regioni=find(strcmp(brain_region_id(:,2),Region));
                    end
                    disp(strcat('Runnning...',path_fishi,'...',Region));
                    %!!!!!! savepath
                    savepath=checkpath(fullfile(path_fishi,'clust_durCond',Region));
                    supervolxel_regioni=supervolxeli(neuron_id_in_regioni,:);
                    act_regioni=activities_preCS_dfdf_aftcorrect_zscore(:,neuron_id_in_regioni);
                    figure;
                    scatter3(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),10,[0.5 0.5 0.5],'filled');hold on;
                    scatter3(supervolxel_regioni(:,1),supervolxel_regioni(:,2),supervolxel_regioni(:,3),10,'r','filled');hold on;legend(Region);axis equal;
                    %figure,    scatter3(supervolxel_regioni(:,1),supervolxel_regioni(:,2),supervolxel_regioni(:,3),10,'filled');hold on;
                    %% para
                    frame_ind_acq=frame.cs_start:frame.us_start-1;%%%%%%%%%%%%%%%%%%%!!!!!!!
                    trial_ind_acq=trial.acq(2):trial.acq(3);
                    A=act_regioni(1:frame.per_cycle*trial.total,:);A_r=reshape(A,frame.per_cycle,trial.total,[]);
                    act_all_acq_block=reshape(A_r(:,trial_ind_acq,:),frame.per_cycle*length(trial_ind_acq),length(neuron_id_in_regioni));
                    M_0=A_r(frame_ind_acq,trial_ind_acq,:);                    act_all_acq_block_CR_AUC=squeeze(sum(M_0,1));
                    M_cluster=reshape(M_0,length(frame_ind_acq)*length(trial_ind_acq),length(neuron_id_in_regioni))';
                    trial_ind_hab=trial.hab(2):trial.hab(3);frame_ind_hab=frame.cs_start:frame.cs_end;
                    act_all_hab=reshape(A_r(:,trial_ind_hab,:),frame.per_cycle*length(trial_ind_hab),length(neuron_id_in_regioni));
                    M_cluster_hab=reshape(A_r(frame_ind_hab,trial_ind_hab,:),length(frame_ind_hab)*length(trial_ind_hab),length(neuron_id_in_regioni))';
                    
                    trial_ind_test=trial.test(2):trial.test(3);frame_ind_test=frame.cs_start:frame.cs_end;
                    act_all_test=reshape(A_r(:,trial_ind_test,:),frame.per_cycle*length(trial_ind_test),length(neuron_id_in_regioni));
                    M_cluster_test=reshape(A_r(frame_ind_test,trial_ind_test,:),length(frame_ind_test)*length(trial_ind_test),length(neuron_id_in_regioni))';
                    
                    figure,subplot(1,2,1),plot(mean(M_cluster,2))
                    Corr_map=corr(M_cluster');subplot(1,2,2),imagesc( Corr_map,[-0.2 0.2]);colorbar;
                    cIX=1:length(neuron_id_in_regioni);
                    save([savepath,'/clust_results.mat'],'Region','neuron_id_in_regioni','supervolxel_regioni','act_regioni','act_all_acq_block_CR_AUC',...
                        'M_cluster','Corr_map','act_all_acq_block',...
                        'act_all_hab','act_all_test','M_cluster_hab' ,'M_cluster_test',...
                        'trial_ind_hab','trial_ind_acq','trial_ind_test','frame_ind_hab','frame_ind_acq','frame_ind_test','-v7.3');
                    dr_type_name={'PCA','t-SNE'};
                    clus_type_name={'hireatch','k-means'};M_cluster_dr=[];
                    %1 降维+层次聚类
                    %                     for dr_type=1:2
                    %                         for clut_type=1:2
                    %                             savepathh=checkpath(fullfile(savepath,[dr_type_name{1},'-',clus_type_name{2}]));
                    %                             if exist(fullfile(savepathh,'logfile.txt'),'file')
                    %                                 delete(fullfile(savepathh,'logfile.txt'));
                    %                             end
                    %                             diary(fullfile(savepathh,'logfile.txt'));
                    %                             [number_cluster,savedclusternum,~,M_cluster_dr,Y_cluster,explained_all,n_cluster]=myclusering(M_cluster',dr_type,clut_type);
                    %                             [gIX,explained]=myclusering_final(M_cluster',M_cluster_dr,clut_type,savedclusternum,5);
                    %                             save([savepathh,'/clust_results.mat'],'M_cluster_dr','cIX','gIX','number_cluster','savedclusternum','explained',...
                    %                                 'n_cluster','-v7.3');
                    %                         end
                    %                     end
                    %2 k-means
                    savepathh=checkpath(fullfile(savepath,'k-means'));
                    if ~exist([savepathh,'/clust_results.mat'],'file')
                        n_minist=8;n_max=32;%savedclusternum;
                        [gIXx,~,number_cluster,~]=find_best_k_in_range(M_cluster,n_minist:n_max);
                        [gIX,explained]=myclusering_final([],M_cluster,2,number_cluster,5);
                        save([savepathh,'/clust_results.mat'],'M_cluster','cIX','gIX','number_cluster','explained','-v7.3');
                    end
                    % 3 k-means of correlation map
                    %                     savepathh=checkpath(fullfile(savepath,'k-means of correlation map'));
                    %                     [gIXx,~,number_cluster,~]=find_best_k_in_range(Corr_map,n_minist:n_max);
                    %                     [gIX,explained]=myclusering_final([],Corr_map,2,number_cluster,5);
                    %                     save([savepathh,'/clust_results.mat'],'M_cluster','cIX','gIX','number_cluster','explained','-v7.3');
                    %
                    %auto cluster
                    try
                        addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
                        savepathh=checkpath(fullfile(savepath,'auto cluster'));
                        if ~exist([savepathh,'/clust_results.mat'],'file')
                            [gIXx,~,number_cluster,~]=find_best_k_in_range(M_cluster, n_minist:n_max);
                            cIX=1:length(neuron_id_in_regioni);
                            cIXx=cIX;cIX_reg =cIX;%ind;
                            if isempty(number_cluster)
                                number_cluster=20;
                            end
                            masterthres=0.5;para=struct;para=setfield(para,'k1',number_cluster);para=setfield(para,'merge',masterthres);para=setfield(para,'cap',masterthres);para=setfield(para,'reg1',masterthres);para=setfield(para,'reg2',masterthres);
                            para=setfield(para,'minSize',3);
                            for ii=1:4
                                [cIX,gIX] = AutoClustering(cIXx,gIXx,M_cluster,cIX_reg,1,para,1,masterthres);
                                %                         ind=[1:5:length(cIX)];%randi(length(cIX),[floor(length(cIX)/5),1]);
                                cIXx=cIX;
                                gIXx=gIX;
                            end
                            number_cluster=length(unique(gIX));
                            save([savepathh,'/clust_results.mat'],'M_cluster','cIX','gIX','number_cluster','-v7.3');
                        end
                    catch
                        warning(['Error in',savepathh]);continue;
                    end
                end
            catch
                warning(['Error in ',path_fishi]);continue;
            end
        end
        close all;
    end
end
%% get brain_region_id singel_fish
for batchi=1:length(Path)
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi}
        %         if ~exist(fullfile(p,'brain_region_related_statistic.mat'),'file')
        try
            %batchi=1;fishi=1;
            %load(fullfile(p,'activities_aft_process.mat'));
            load(fullfile(path_fishi,'para.mat'));
            load(fullfile(path_fishi,'env.mat'));
            nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
            supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
            % supervolxeli(:,3)=[430-supervolxeli(:,3)];
            %% get region_id
            surpervolxel=supervolxeli;gIX=[ones(size(surpervolxel,1),1)];
            clrmap=jet(10);
            [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,surpervolxel,reg_mask,reg_name,reg_loc,temp_supervoxel,temp_env,clrmap,nn);
            save(fullfile(path_fishi,'brain_region_related_statistic.mat'),'loc_in_region_in_clust','fraction_in_region_in_clust','loc_in_region_cell','num_in_region_in_clust','brain_region_id','Label','-v7.3');
        catch
            warning([path_fishi,' not extract brain region']);
        end
        %         end
    end
end

%% write to csv
for batchi=1:length(Path)
    for fishi=1:length(Path{batchi})
        path_fishi=Path{batchi}{fishi}
%        if ~exist(fullfile(path_fishi,'vol_env_spatialloc_warped_SyN_add_brainregion.csv'),'file')
            try
                nn=[path_fishi(end-14:end-7),path_fishi(end-5:end-1)];
                supervolxeli=readmatrix([warped_SyN_csv_path ,nn,'vol_env_spatialloc_warped_SyN.csv']);
                load(fullfile(path_fishi,'brain_region_related_statistic.mat'))
                columns = {'x', 'y', 'z', 't','index','brain_region_id','label'};
                if length(brain_region_id(:,1))~=length(supervolxeli(:,1))
                    brain_region_id(end+1:length(supervolxeli(:,1)),1)='';
                end
                M=table(supervolxeli(:,1),supervolxeli(:,2),supervolxeli(:,3),supervolxeli(:,4),supervolxeli(:,5),brain_region_id(:,1),brain_region_id(:,2),'VariableNames', columns);
                outpath=fullfile(path_fishi,'vol_env_spatialloc_warped_SyN_add_brainregion.csv');
                writetable(M,outpath,'Encoding','UTF-8');
            catch
                warning([path_fishi,' not write csv']);
            end
%         else
%             warning([fullfile(path_fishi,'vol_env_spatialloc_warped_SyN_add_brainregion.csv'),'------already exist'])
%         end
    end
end
           
