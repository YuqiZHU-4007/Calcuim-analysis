function showspv=mapback_to_vol_zyq_20180520(env,opt,index,drawtype)

[spvzind, spvregion] = supervoxel_index(env, opt);
spv = env.supervoxel;

%% segmentation result
% svol = zeros([env.height env.width env.depth]);
% if strcmp(drawtype, 'dot')
%     showspv = svol;
% elseif strcmp(drawtype, 'circle')
%     showspv = zeros(env.height, env.width, 3, env.depth);
% end
showspv=[];
for zi = 1:env.depth
    svol = rescalegd(im2double(env.vol(:,:,zi)), [1/10000 1/10000]);
    slice = svol;
    pti = spvzind{zi};
    pti=intersect(index,pti);
    if strcmp(drawtype, 'dot')
        %slice(ind) = 1;  %1;
        %ind = sub2ind(size(slice), spv(pti, 2), spv(pti, 1));
        for ii = 1:numel( pti)
            spi = pti(ii);
            slice(spvregion{spi}, spv(spi,3)) = 1;
        end
        showspv(:,:,zi) = slice;
    elseif strcmp(drawtype, 'circle')
        slice = insertShape(slice, 'circle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1],'linewidth',3);
        showspv(:,:,:,zi) = slice;
    end
end
end