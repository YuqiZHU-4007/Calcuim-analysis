function [act_all,act_all_hab,act_all_acq,act_all_acq_block_mean,act_all_tst,env_all,index_all]=Get_all_act(K,colnum,cut_move,is_eatra_act_all,actbatch,envbatch,behavbatch,parabatch,savepath_all)
%extra all acti
load(parabatch{K(1)});
num_trial_hab=6;num_trial_tst=6;
% n=200000;
% act_all=nan(n*size(actbatch,1),colnum);index_all=nan(n*size(actbatch,1),1);
% act_all_hab=nan(frame.per_cycle,num_trial_hab,n*size(actbatch,1));
% act_all_acq=nan(frame.per_cycle,trial.acq_block_trial,trial.acq_block_num,n*size(actbatch,1));
% act_all_acq_block_mean=nan(frame.per_cycle,trial.acq_block_num,n*size(actbatch,1));
% act_all_tst=nan(frame.per_cycle,num_trial_tst,n*size(actbatch,1));
%env_all.supervoxel=nan(n*size(actbatch,1),3);

act_all=[];index_all=[];
act_all_hab=[];
act_all_acq=[];
act_all_acq_block_mean=[];
act_all_tst=[];
env_all.supervoxel=[];
%act_all_tst_mean=[];act_all_hab_mean=[];
%area_all=struct;
for kk=K;%1:size(actbatch,1);
    load(actbatch{kk});
    load(behavbatch{kk});
    load(parabatch{kk});
    act=activities_preCS_dfdf_aftcorrect';
    index_all=cat(1,index_all, [ones(1,size(act',2))*kk]');
    load(envbatch{kk});
    env_all.supervoxel=cat(1,env_all.supervoxel,env.supervoxel);
    if is_eatra_act_all
        if size(act,2)>colnum %size(act,2)~=size(act_all,1) && ~isempty(act_all)
            act(:,colnum+1:end)=[];
            warning([actbatch{kk} ' acti was cut into ' num2str(size(act,2))]);
        end
        act_all=cat(1,act_all,act);
    end
    i_ind=3;
    %     ind_cut_trial_preCS=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start-5*fs.behavior) & re_startpoint(:,2)<(frameb.cs_start)),1);
    %     ind_cut_trial_preCS=unique(ind_cut_trial_preCS);
    %     ind_cut_trial_CS=re_startpoint(find(re_startpoint(:,2)>=(frameb.cs_start) & re_startpoint(:,2)<(frameb.us_start-1)),1);%%%%%%%%%%%%%%%%
    %     ind_cut_trial_CS=unique(ind_cut_trial_CS);
    %     ind_cut_trial=union(ind_cut_trial_preCS,ind_cut_trial_CS);
    %    [rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial);
    %     if ~exist('rawacti') && i_ind>1 && ~exist('area')
    %         [rawacti,area]=calculate_integtate_dfdf_main(activities_preCS_dfdf_aftcorrect,cut_move,trial,frame,fs,ind_cut_trial_preCS);
    %     end
    for ii=1:i_ind
        switch ii
            case 1
                trial_ind=trial.hab(2):trial.hab(2)+num_trial_hab-1; %%trial.test(2):trial.test(2)+2
                act_CS= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,[],size(act,1));
                %act_all_hab(:,:,end+1:end+size(act_CS_hab,3))=act_CS_hab;
                act_all_hab=cat(3,act_all_hab,act_CS);
                %act_CS_hab_mean=reshape(mean(act_CS,2),frame.per_cycle,[]);
                %act_all_hab_mean=[act_all_hab_mean,act_CS_hab_mean];
            case 2
                trial_ind=trial.acq(2):trial.acq(3); %%trial.test(2):trial.test(2)+2
                act_CS= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,trial.acq_block_trial,trial.acq_block_num,size(act,1));
                act_all_acq=cat(4,act_all_acq,act_CS);
                act_all_acq_block_mean=cat(3,act_all_acq_block_mean,squeeze(mean(act_CS,2)));
%                 for jj=1:size(act',2)
%                     act_all_acq_block_mean(:,:,size(act_all_acq_block_mean,3)+1)=rawacti.acq_mean_blocks{jj};
%                 end
                %                 name=fieldnames(area);
                %                 for jj=1:length(name)
                %                     if ~isfield(area_all,name{jj})
                %                         area_all=setfield(area_all,name{jj},[]);
                %                     end
                %                     a=getfield(area,name{jj});
                %                     aa=getfield(area_all,name{jj});
                %                     area_all=setfield(area_all,name{jj},[aa a]);
                %                 end
            case 3
                trial_ind=trial.test(2):trial.test(2)+num_trial_tst-1; %%trial.test(2):trial.test(2)+2
                act_CS= reshape(act(:,(trial_ind(1)-1)*frame.per_cycle+1:trial_ind(end)*frame.per_cycle)',frame.per_cycle,[],size(act,1));
                act_all_tst=cat(3,act_all_tst,act_CS);
                %act_CS_hab_mean=reshape(mean(act_CS,2),frame.per_cycle,[]);
                %act_all_tst_mean=[act_all_tst_mean,act_CS_hab_mean];
        end
        %act_CS_hab_mean=zeros(frame.per_cycle,size(act,1));
    end
end
%act_all(isnan(act_all(:,1)))=[];
%act_all_acq_block_mean(:,:,1)=[];act_all_tst(:,:,1)=[];
save([savepath_all '\act'],'act_all','index_all','env_all','act_all','act_all_hab','act_all_acq','act_all_acq_block_mean','act_all_tst','-v7.3');

