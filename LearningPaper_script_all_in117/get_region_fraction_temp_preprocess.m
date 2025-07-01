function [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,brain_region_id,loc_in_region_temp,id_in_region_temp,num_in_region_in_clust,Label]=get_region_fraction_temp_preprocess(gIX,ind_supervoxel,supervoxel,reg_mask,reg_name,temp_supervoxel,clrmap)


supervoxel_raw=supervoxel;index_in_region_in_clust=nan(length(reg_name),1);
loc_in_region_in_clust=nan(length(reg_name),1);brain_region_id={}; Label=[];
fraction_in_region_in_clust=nan(length(reg_name),1);
num_in_region_in_clust{1}=nan(length(reg_name),1);
num_in_region_in_clust{2}=nan(length(reg_name),1);
num_in_region_in_clust{3}=nan(length(reg_name),1);
loc_in_region_temp={};id_in_region_temp={};
if ~isempty(supervoxel)
    clr= clrmap(unique(gIX),:);
%     [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust,Label]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
%         gIX,supervoxel,temp_supervoxel,clr,false,nn);
    [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,brain_region_id,loc_in_region_temp,id_in_region_temp,num_in_region_in_clust,Label]=get_region_fraction_temp_brief(gIX,ind_supervoxel,supervoxel,temp_supervoxel,reg_mask,reg_name,clr,false);
    %% brain_id_map
%     figure,
%     clr=jet(31);
%     a=str2double(brain_region_id(:,1));b=isnan(a);a(b)=31;
%     scatter3(supervoxel_raw(a,1),supervoxel_raw(a,2),supervoxel_raw(a,3),10,clr(a(:),:),'filled');hold on;
%     %  load('X:\calcium data 20230224\20210709\fish2\env.mat')
%     % figure,scatter3(env.supervoxel(:,1),env.supervoxel(:,2),env.supervoxel(:,3),5,[0.5 0.5 0.5],'filled');hold on;
%     for ii=1:length(index_in_region_in_clust);
%         ind=index_in_region_in_clust{ii};
%         %          scatter3(supervoxel_raw(:,1),supervoxel_raw(:,2),supervoxel_raw(:,3),5,[0.5 0.5 0.5],'filled');hold on;
%         %          scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,clr(a(ind),:),'filled');hold on;
%         %                   figure,scatter3(supervoxel_raw(:,1),supervoxel_raw(:,2),supervoxel_raw(:,3)  ,5,[0.5 0.5 0.5],'filled');hold on;
%         %          scatter3(loc_in_region_in_clust{ii}(:,1), loc_in_region_in_clust{ii}(:,2), loc_in_region_in_clust{ii}(:,3),10,clr(a(ind),:),'filled');hold on;
%         scatter3(env.supervoxel(ind,1),env.supervoxel(ind,2),env.supervoxel(ind,3),10,clr(a(ind),:),'filled');hold on;
%     end
    % for ii=1:length(index_in_region_in_clust);
    %     a=index_in_region_in_clust{ii};ind=[];
    %     for zz=1:length(a)
    %         ind(zz,1)=find(supervoxel_raw(:,1)==supervoxel(a(zz),1) & supervoxel_raw(:,2)==supervoxel(a(zz),2) & supervoxel_raw(:,3)==supervoxel(a(zz),3));
    %     end
    %     brain_region_id(ind,1)=ii;
    %     brain_region_id(ind,2)=string(Label{ii});
    %     scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,repmat(clr(ii,:),length(ind),1),'filled');hold on;
    % end
    % title(nn);legend(Label);axis equal;
    % close all;
end
end