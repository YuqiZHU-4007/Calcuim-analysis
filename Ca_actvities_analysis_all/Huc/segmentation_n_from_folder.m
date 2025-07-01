%单独分割
inputpath='C:\Users\pc\Desktop\20190521_331FVH_seedfish04_DXF_G.tif';
opt.maxrad=7;
opt.minrad=2;
opt.thres=[];

% listdir = dir([inputpath '/*.tif']);
% for iii = 1:length(listdir)
%     filename = listdir(iii).name;
%     %fileid = str2double(filename(24:end-4));
%     env.vol(:,:,iii) = tiffread([inputpath  '/' filename]);
% end
env.vol=tiffread(inputpath);
env.height=size(env.vol(:,:,:),1);
env.width=size(env.vol(:,:,:),2);
env.depth=size(env.vol(:,:,:),3);
show_spv_GUI(env.vol);
%backgroung_segmentarion_n
ZZ=1:size(env.vol,3);
multithreshcount = 2;    
slice = max(im2double(env.vol(:,:,ZZ)), [], 3);
slice = rescalegd(slice, [0, 1/100]);  
pslice = imfilter(slice, fspecial('gaussian', opt.maxrad*2, opt.maxrad/2));
level = multithresh(pslice, 2);
level = level(1);
mask = im2bw(pslice, level);
figure; imshow(mask);
%env.volmask(:,:,ZZ)=mask;
smvol = zeros(env.height, env.width, env.depth);
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    pslice = imfilter(slice, fspecial('gaussian', opt.minrad*2, opt.minrad/2));
    smvol(:,:,zi) = pslice;
end
try
    level = multithresh(smvol(:), 3);
catch
    level = multithresh(smvol(:), 2);
end
level = level(1);
volmask = false(size(env.vol));
for zi = 1:env.depth
    slice = smvol(:,:,zi);
    bw = slice > level;
    bw = bw & mask;
    volmask(:,:,zi) = bw;
end
env.volmask=volmask;env.mask=mask;
show_spv_GUI(env.volmask);
%supervoxel_segmentation_n
minrad = opt.minrad;
maxrad = opt.maxrad;
thres = opt.thres;
supervoxel = [];
scoremaps = zeros(env.height, env.width, env.depth);
radmaps   = zeros(env.height, env.width, 2, env.depth);
show = env.vol;
%normvol = normim(env.vol, [4/2048/2048 4/2048/2048]);  % normalize together
normvol = normim(env.vol, [0 0]);  % normalize together
for zi = 1:env.depth
    zi
    slice = normvol(:,:,zi);    
    % smooth  smooth width should be smaller than cell diameter
    gs = fspecial('gaussian', maxrad, 1);
    smslice = imfilter(slice, gs);
    %% detect
    radii = minrad:maxrad;
    [points, scoremap, radmap] = detect_by_ring_2d_n(smslice, radii, thres, env.volmask(:,:,zi));
    scoremaps(:,:,zi) = scoremap;
    radmaps(:,:,:,zi) = radmap;    
    %% rois
    [yy, xx] = find(points);
    center = [xx, yy];
    rad = mean(radmap, 3);
    rad = rad(points);
    rin = radmap(:,:,1);
    rin = rin(points);
    rex = radmap(:,:,2);
    rex = rex(points);
    scores = scoremap(points);   
    %% remove overlap
    size(center,1) %%分割出多少个
    if ~isempty(center)
    idx = remove_overlap(center, rin, scores);
    length(idx)
    center = center(idx,:);
    rad = rad(idx);
    rex = rex(idx);
    rin = rin(idx);   
    end
    %%
    supervoxel = [supervoxel; center zi*ones(length(rin),1) rin rex];
    slice = insertShape(zeros(size(env.vol(:,:,zi))), 'circle', [supervoxel(:, 1),supervoxel(:, 2), supervoxel(:, 4)],'linewidth',3);
    slice(~env.volmask(:,:,zi))=0;
     show_spv_GUI(slice);
    imwrite(slice(:,:,1),fullfile(['C:\Users\pc\Desktop\ ','slice', num2str(zi,'%03d') ,'.tif']))
end
env.supervoxel=supervoxel;env.scoremaps=scoremaps; env.radmaps=radmaps;
showspv = show_spv(env, opt, 'circle');
seqwrite(showspv, checkpath([inputpath '/spvseg_circle/']));
save([inputpath '/supervoxel_index.mat'], 'env', 'opt','showspv',  '-v7.3');