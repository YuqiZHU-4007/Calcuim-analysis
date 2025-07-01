function spot_spv(pti, env)
% spot_spv(pti, env)

zi = env.supervoxel(pti, 3);
slice = normim(env.vol(:,:,zi), [1/10000, 1/10000]);

slice = insertShape(slice, 'circle', [env.supervoxel(pti, 1), env.supervoxel(pti, 2), floor(env.supervoxel(pti, 5)) - 1], 'color', 'r');

figure; imshow(slice);

