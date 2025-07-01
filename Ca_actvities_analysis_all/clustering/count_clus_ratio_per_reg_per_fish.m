cIX_iii=cIX;gIX_iii=gIX;clrmap_iii=clrmap;a=[]; isloadmask=true;
ratio_fish=zeros(length(unique(gIX_iii)),10,size(actbatch,1));
num_cell_region=[];
savepathi='E:\A_Data_lightsheet\Data_vmat\Summary\Calcium summary\integrate dfdf\all_learner_huc_cutmove_20190619\huc_learner_cutmove_20190617\path\mask\';
num_reg=5;
for kk=unique(index_all);%fish_num
    load(envbatch{kk});
    if isloadmask
        load([savepathi '\region_mask' num2str(kk,'%02d')]);
        %mask=ones(2048,2048);
    else
        temp=max(env.vol(:,:,1:24),[],3);
        %h=figure;imshow(temp,[min(temp(:)) max(temp(:))]);
        %num_reg=input('input number of region\n');
        %close(h);
        mask=getmask_imfreehand(env.vol(:,:,1:24),num_reg);
        save([savepathi '\region_mask' num2str(kk,'%02d')],'mask');
        figure,imshow(mask',[min(mask(:)) max(mask(:))]);
    end
    maskk=reshape(mask,1,[]);
    num_reg=length(unique(mask))-1;
    [ind_event_in_this_fish,i_index_all,i_cIX]=intersect(find(index_all==kk) ,cIX_iii);
    ind_env_in_this_fish=find(index_all==kk);
    ind_cls_in_this_fish=i_cIX;
    for ii=unique(gIX(ind_cls_in_this_fish))';%clus_num
        for jj=1:num_reg
            ind=sub2ind(size(mask),env.supervoxel(:,2),env.supervoxel(:,1));
            ind_region_in_this_fish=find(mask(ind)==jj);
%                         figure,imagesc(mask==jj);hold on;
%                         scatter(env.supervoxel(ind_region_in_this_fish,1),env.supervoxel(ind_region_in_this_fish,2));hold on
%                         scatter(env_all.supervoxel(ind_env_in_this_fish(ind_region_in_this_fish),1),env_all.supervoxel(ind_env_in_this_fish(ind_region_in_this_fish),2));hold on
            ratio_fish(ii,jj,kk)=length(intersect(cIX(ind_cls_in_this_fish(find(gIX(ind_cls_in_this_fish)==ii))),ind_env_in_this_fish(ind_region_in_this_fish)));
%                         figure,imagesc(mask==jj);hold on;
%                         scatter(env_all.supervoxel(cIX(ind_cls_in_this_fish(find(gIX(ind_cls_in_this_fish)==ii))),1),env_all.supervoxel(cIX(ind_cls_in_this_fish(find(gIX(ind_cls_in_this_fish)==ii))),2));hold on
%                         scatter(env_all.supervoxel(ind_env_in_this_fish(ind_region_in_this_fish),1),env_all.supervoxel(ind_env_in_this_fish(ind_region_in_this_fish),2));hold on
        end
    end
        for jj=1:num_reg
            ind=sub2ind(size(mask),env.supervoxel(:,2),env.supervoxel(:,1));
            num_cell_region(jj,1,kk)=length(find(mask(ind)==jj));
        end
end

ratio_fish(:,6:10,:)=[];ratio_fish(:,:,1)=[];ratio_fish(:,:,8)=[];
a=sum(ratio_fish,3)/length(cIX);
a=sum(ratio_fish,3);b=sum(num_cell_region,3);
a=a./repmat(b',size(a,1),1);

%%
a=sum(ratio_fish(:,:,[3,4,5,7,8,9]),3);a(:,6:10)=[];
b=sum(ratio_fish(:,:,[12,13]),3);b(:,6:10)=[];
a=a./repmat([sum(num_cell_region(:,:,[3,4,5,7,8,9]),3)]',length(unique(gIX)),1);
b=b./repmat([sum(num_cell_region(:,:,[12,13]),3)]',length(unique(gIX)),1);
% a=a./repmat([sum(num_cell_region(:,:,:),3)]',length(unique(gIX)),1);
% b=b./repmat([sum(num_cell_region(:,:,:),3)]',length(unique(gIX)),1);
% a=a./sum([sum(num_cell_region(:,:,[3,4,5,7,8,9]),3)]');
% b=b./sum([sum(num_cell_region(:,:,[12,13]),3)]');
% a=a./repmat([sum(num_cell_region(:,:,[3,9]),3)]',length(unique(gIX)),1);
% b=b./repmat([sum(num_cell_region(:,:,[13]),3)]',length(unique(gIX)),1);
max([a(:);b(:)])
ylim=[0 0.25];
figure,subplot(1,2,1);imagesc(a,ylim);colorbar;
subplot(1,2,2);imagesc(b,ylim);colorbar;

