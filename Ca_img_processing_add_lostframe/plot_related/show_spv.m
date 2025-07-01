function showspv = show_spv(env, opt, drawtype)
% showspv = show_spv(env, opt, drawtype = 'dot')
%
% drawtype can be 'dot' or 'circle'
% 

[spvzind, spvregion] = supervoxel_index(env, opt);
spv = env.supervoxel;

%% segmentation result
svol = rescalegd(im2double(env.vol), [1/10000 1/10000]);
% svol = zeros([env.height env.width env.depth]);
if strcmp(drawtype, 'dot')
    showspv = svol;
elseif strcmp(drawtype, 'circle')
    showspv = zeros(env.height, env.width, 3, env.depth);
end

for zi = 1:env.depth
    slice = svol(:,:,zi);
    pti = spvzind{zi};
    if strcmp(drawtype, 'dot')
        ind = sub2ind(size(slice), spv(pti, 2), spv(pti, 1));
        slice(ind) = 1;  %1;
        showspv(:,:,zi) = slice;
    elseif strcmp(drawtype, 'circle')
        slice = insertShape(slice, 'circle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1]);
        showspv(:,:,:,zi) = slice;
    end
end

