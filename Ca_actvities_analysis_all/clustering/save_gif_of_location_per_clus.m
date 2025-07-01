savepathi='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\figures_masterthres0.6_01';
savepathi='K:\2.Fear conditioning_huc_dimming\autoclustering_results\figures_masterthres0.6_01';

%%  per fish per cluster
clrmap_iii=clrmap;
num_clust=[3,14];cIX_iii=cIX;gIX_iii=gIX;a=[];
for ii=num_clust
    for kk=unique(index_all);%fish_num
        [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
        load(envbatch{kk});
        ind_cls_in_this_fish=i_cIX;ind_add=[];
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
        ind_cls_in_this_fish=i_cIX(ind_add);
        colorind=unique(gIX_iii(ind_cls_in_this_fish));
        showspv=mapback_show_zyq_20190729(repmat(clrmap_iii(colorind,:),length(ind_cls_in_this_fish),1), env_all.supervoxel,cIX_iii(ind_cls_in_this_fish),env.vol(:,:,1:24),'filledcircle');
        %repmat(clrmap_iii(colorind,:),length(ind_cls_in_this_fish),1)
        %h=show_spv_GUI(showspv);%seqwrite(showspv,[savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)) '\'],'gif');
        filename=[savepathi,'\fish' num2str(kk,'%02d') '_' num2str(ii,'%02d'),'.gif']; % Specify the output file name
        for idx = 1:size(showspv,4)
            [A,map] = rgb2ind(showspv(:,:,:,idx),256);
            if idx == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
            end
        end
    end
end

%% per fish all cluster
clrmap_iii=clrmap;
num_clust=[3,14];cIX_iii=cIX;gIX_iii=gIX;a=[];
for kk=unique(index_all)
    [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
    load(envbatch{kk});
    ind_cls_in_this_fish=i_cIX;ind_add=[];colorind=[];
    for ii=num_clust
        ind_add=[ind_add;find(gIX_iii(i_cIX)==ii)];
        colorind=[colorind;ii*ones(size(find(gIX_iii(i_cIX)==ii)))];
        %ind_cls_in_this_fish( find(gIX_iii(i_cIX)~=ii))=[];
    end
    ind_cls_in_this_fish=i_cIX(ind_add);
    showspv=mapback_show_zyq_20190729(clrmap_iii(colorind,:), env_all.supervoxel,cIX_iii(ind_cls_in_this_fish),env.vol(:,:,1:24),'filledcircle');
    %figure;h=show_spv_GUI(showspv);%seqwrite(showspv,[savepath 'stim avr corr mapping' num2str(trial_ind(1)) '-' num2str(trial_ind(end)) '\'],'gif');
    filename=[savepathi,'\fish' num2str(kk,'%02d') '_' num2str(num_clust,'%02d'),'.gif']; % Specify the output file name
    for idx = 1:size(showspv,4)
        [A,map] = rgb2ind(showspv(:,:,:,idx),256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end
