%% 按学习过程中CR聚类
savepathh=checkpath('H:\3.Juvenile reference brain\registration to templete\data\2021_ind\cluter_durCond\');
%单条鱼
addpath(genpath('F:\DUlab\FC analyse\FishExplorer'));
load('H:\1.Test US\2.Tail free——Data from 117\20210805\fish1\para.mat')
ref_win=ceil(0/fs.ca+1):frame.cs_start-1;%取时间段的baseline作为参考判断event
time=([1:frame.per_cycle]-frame.cs_start)*fs.ca;
area_win_hab=frame.cs_start:frame.cs_end-1;
stimCS=zeros(1,frame.per_cycle*trial.total(1));
for ii=1:size(stimCS,2)/frame.per_cycle
    stimCS((ii-1)*frame.per_cycle+frame.cs_start:(ii-1)*frame.per_cycle+frame.cs_end)=0;%23
end
stimUS=ones(1,frame.per_cycle*trial.acq(1))*3;
for ii=trial.acq(2):trial.acq(3)
    stimUS((ii-1)*frame.per_cycle+frame.us_start:(ii-1)*frame.per_cycle+frame.us_start+2)=0;%16
end
%% 
inputpath='H:\3.Juvenile reference brain\registration to templete\data\2021_ind\';
for zz=1:length(seg1)
    for  hh=[1 2 4]
        if strcmp(seg1{zz},'UNRESPONSIVE') && strcmp(seg2{hh},'shift to other')
            continue;
        end
        M_cluster=[];act_all_acq_block=[];supervoxel={};supervoxel_i=[];M_cluster_fishi=[];
        for batchi=1:3
            path=Path{batchi};aa=string(path)';aa=strrep(aa,'\','');
            csv_file=[checkpath(fullfile(inputpath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);
            for name=string(unique(M.label))'
                if find(contains(name,'2019'))
                    continue;
                elseif find(contains(name,'20210430fish1'))
                    continue;
                end
                fish_i= find(contains(aa,name));p=path{fish_i};
                ind_in_fishi=M.ind(find(contains(M.label,name)));
                load(fullfile(p,'para.mat'));load(fullfile(p,'env.mat'));
                trial_ind_acq=trial.acq(2):trial.acq(3);frame_ind=frame.cs_start:frame.us_start-1;
                if length(trial_ind_acq)~=8*3
                    continue;
                end
                load(fullfile(p,'activities_aft_process.mat'));
                Z =normalize(activities_preCS_dfdf_aftcorrect,1);
                A=Z(1:frame.per_cycle*trial.total,ind_in_fishi);A_r=reshape(A,frame.per_cycle,trial.total,[]);
                M_cluster=cat(2,M_cluster,reshape(A_r(frame_ind,trial_ind_acq,:),[],length(ind_in_fishi)));
                act_all_acq_block=cat(2,act_all_acq_block,reshape(A_r(:,trial_ind_acq,:),[],length(ind_in_fishi)));
                supervoxel_i=cat(1,supervoxel_i,[M.x(find(contains(M.label,name))) M.y(find(contains(M.label,name))) M.z(find(contains(M.label,name)))]);
                M_cluster_fishi=cat(2,M_cluster_fishi,[batchi*ones(1,length(ind_in_fishi));repmat(name,1,length(ind_in_fishi)),ind_in_fishi]);
            end
        end
        save([savepath,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'M_cluster','act_all_acq_block','supervoxel_i','M_cluster_fishi','-v7.3');  
    end
end
    
    
%% re-load
inputpath='H:\3.Juvenile reference brain\registration to templete\data\2021_ind\';
savepath='H:\3.Juvenile reference brain\registration to templete\data\2021_ind\cluter_durCond\';
for zz=1:length(seg1)
    for  hh=[1 2 4]
        if strcmp(seg1{zz},'UNRESPONSIVE') && strcmp(seg2{hh},'shift to other')
            continue;
        end
        load([savepath,'cluster_type_',seg1{zz},seg2{hh},'.mat']);
        M_cluster_fishi=[];
        for batchi=1:3
            path=Path{batchi};aa=string(path)';aa=strrep(aa,'\','');
            csv_file=[checkpath(fullfile(inputpath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\',seg2{hh},'.csv']
            M=readtable(csv_file);
            for name=string(unique(M.label))'
                if find(contains(name,'2019'))
                    continue;
                elseif find(contains(name,'20210430fish1'))
                    continue;
                end
                ind_in_fishi=M.ind(find(contains(M.label,name)));
                %A{kk}=[batchi*ones(1,length(ind_in_fishi));repmat(name,1,length(ind_in_fishi));ind_in_fishi'];kk=kk+1;
                M_cluster_fishi=cat(2,M_cluster_fishi,[batchi*ones(1,length(ind_in_fishi));repmat(name,1,length(ind_in_fishi));ind_in_fishi']);
            end
        end
        save([savepath,'cluster_type_',seg1{zz},seg2{hh},'.mat'],'M_cluster','act_all_acq_block','supervoxel_i','M_cluster_fishi','-v7.3');
    end
end

%% PCA
