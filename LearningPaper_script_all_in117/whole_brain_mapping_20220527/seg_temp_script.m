env.vol = uint16(nrrdread('H:\3.Juvenile reference brain\registration to templete\template4_20210709_Cut_Ds_0.66_0.66_10.nrrd'));
opt.imgpath='H:\3.Juvenile reference brain\registration to templete\template4_20210709_Cut_Ds_0.66_0.66_10.nrrd';
env.T=size(env.vol,3);env.depth=size(env.vol,3);
opt.savepath='H:\3.Juvenile reference brain\registration to templete\';
%% f-b segmentation         001.0

[env.volmask, env.mask] = background_segmentation_nucleus(env, opt);
tiffwrite(im2uint8(env.volmask), [opt.savepath '/volmask.tiff']);
tiffwrite(im2uint8(cat(2, rescalegd(env.vol, [0 0]), env.volmask)), [opt.savepath '/inspect_volmask.tiff']);

%% supervoxel segmentation

[env.supervoxel, scoremaps, radmaps] = supervoxel_segmentation_n(env, opt);

%% index supervoxel
[spvzind, spvregion] = supervoxel_index(env, opt);
save([opt.savepath '/supervoxel_index.mat'], 'spvzind', 'spvregion');
%% save segmentation result
showspv = show_spv(env, opt, 'dot');
svpath = checkpath([opt.savepath '/spvseg_dot/']);
seqwrite(showspv, svpath);
showspv = show_spv(env, opt, 'circle');
svpath = checkpath([opt.savepath '/spvseg_circle/']);
seqwrite(showspv, svpath);
%     save([opt.savepath '/env.mat'], 'env', 'opt', '-v7.3');
save([opt.savepath '/env.mat'], 'env', 'opt');
    
