function temp_plot(map_supervoxel,temp_supervoxel,temp,seg,ploti,savepath,nn,colormap)
%map_supervoxel:[x y z]

if ~isempty(strfind(nn,'2019'))
    res=[0.66,0.66,8];
elseif ~isempty(strfind(nn,'2021'))
    res=[0.66,0.66,10];
end
switch ploti
    case 1
        clr=[1 0 0];
        showspv=[];showspv(:,:,1,:)=temp;showspv(:,:,2,:)=temp;showspv(:,:,3,:)=temp;showspv=showspv/255;
        for ii=1:size(temp,3)
            id=find(floor(map_supervoxel(:,3)/res(3))==ii);A=[map_supervoxel(id,1)./res(1),map_supervoxel(id,2)./res(2),8*ones(length(id),1)];
            slice = insertShape(showspv(:,:,:,ii), 'filledcircle', A,'color',clr);
            showspv(:,:,:,ii)=slice;
        end
        %show_spv_GUI(showspv);title(seg);
        save_img(showspv,seg,savepath);%seqwrite(showspv, fullfile(checkpath(fullfile(outpath,seg_batchi{batchi},seg1{zz},seg2{hh})),'\mapback.tif'),'tif')
    case 2
        colormap=colormap+1;
        clr=hot(ceil(ceil(max(colormap)*1.8)));clr_2=clr(colormap,:);
        showspv=[];showspv(:,:,1,:)=temp;showspv(:,:,2,:)=temp;showspv(:,:,3,:)=temp;showspv=showspv/255;
        for ii=1:size(temp,3)
            id=find(floor(map_supervoxel(:,3)/res(3))==ii);A=[map_supervoxel(id,1)./res(1),map_supervoxel(id,2)./res(2),8*ones(length(id),1)];
            slice = insertShape(showspv(:,:,:,ii), 'filledcircle', A,'color',clr_2(id,:));
            showspv(:,:,:,ii)=slice;
        end
        %show_spv_GUI(showspv);title(seg);
        save_img(showspv,seg,savepath);
    case 3
        showspv_zproject=squeeze(sum(showspv,4));size(showspv_zproject)
        figure,imshow(showspv_zproject,[min(showspv_zproject(:)),max(showspv_zproject(:))]);title([seg1{zz},seg2{hh}])
    case 4
        v=[90,90];
        radium=floor(15/0.66);
        count_spatial_location_density(temp_supervoxel,temp_supervoxel(:,:), radium);view(v(1),v(2));title(seg);
end