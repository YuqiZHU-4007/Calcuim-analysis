function [loc_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,num_in_region_in_clust,brain_region_id,Label]=get_region_fraction_temp_preprocess(gIX,supervoxel,reg_mask,reg_name,reg_loc,temp_supervoxel,temp_env,clrmap,nn)
if ~isempty(strfind(nn,'2019'))
    res2=[0.66,0.66,8];
elseif ~isempty(strfind(nn,'2021')) | ~isempty(strfind(nn,'2022'))
    res2=[0.66,0.66,10];
end
supervoxel_raw=supervoxel;
loc_in_region_in_clust=nan(length(reg_name),1);
fraction_in_region_in_clust=nan(length(reg_name),1);loc_in_region_cell=[];num_in_region_in_clust{1}=nan(length(reg_name),1);
num_in_region_in_clust{2}=nan(length(reg_name),1);num_in_region_in_clust{3}=nan(length(reg_name),1);
if ~isempty(supervoxel)
    gIX(find(supervoxel(:,1)<=0 | supervoxel(:,1)>size(temp_env.env.vol,1)*res2(1)))=[];
    supervoxel(find(supervoxel(:,1)<=0 | supervoxel(:,1)>size(temp_env.env.vol,1)*res2(1)),:)=[];
end
if ~isempty(supervoxel)
    gIX(find(supervoxel(:,2)<=0 | supervoxel(:,2)>size(temp_env.env.vol,2)*res2(2)))=[];
    supervoxel(find(supervoxel(:,2)<=0 | supervoxel(:,2)>size(temp_env.env.vol,2)*res2(2)),:)=[];
end
if ~isempty(supervoxel)
    gIX(find(supervoxel(:,3)<=0 | supervoxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)))=[];
    supervoxel(find(supervoxel(:,3)<=0 | supervoxel(:,3)>(size(temp_env.env.vol,3)-1)*res2(3)),:)=[];
end
if ~isempty(supervoxel)
    clr= clrmap(unique(gIX),:);
    [loc_in_region_in_clust,index_in_region_in_clust,fraction_in_region_in_clust,loc_in_region_cell,~,num_in_region_in_clust,Label]=get_region_fraction_temp(reg_mask,reg_name,reg_loc,...
        gIX,supervoxel,temp_supervoxel,clr,false,nn);
end
%% brain_id_map
brain_region_id=string;    
figure,
clr=jet(30);
scatter3(supervoxel_raw(:,1),supervoxel_raw(:,2),supervoxel_raw(:,3),10,[0.5 0.5 0.5],'filled');hold on;
for ii=1:length(index_in_region_in_clust);
    a=index_in_region_in_clust{ii};ind=[];
    for zz=1:length(a)
        ind(zz,1)=find(supervoxel_raw(:,1)==supervoxel(a(zz),1) & supervoxel_raw(:,2)==supervoxel(a(zz),2) & supervoxel_raw(:,3)==supervoxel(a(zz),3));
    end
    brain_region_id(ind,1)=ii;
    brain_region_id(ind,2)=string(Label{ii});
    scatter3(supervoxel_raw(ind,1),supervoxel_raw(ind,2),supervoxel_raw(ind,3),10,repmat(clr(ii,:),length(ind),1),'filled');hold on;
end
title(nn);legend(Label);axis equal;
close all;

end