function [supervoxel, scoremaps, radmaps] = supervoxel_segmentation_cytosol(env, opt)
% supervoxel_segmentation(env, opt)
%   
% outputs:
%   supervoxel: a n_supervoxel * 5 matrix, each row correspond to one supervoxel
%       a supervoxel is represented in a 5 element vector [x y z inner_radius outer_radius]
%   scoremap:   an array with the size of the volume
%       each position represents the score of it being a cell center
%   radmaps:    an array with the size [height, width, 2, depth]
%       each position represents the optimal radius (inner and outer), if this
%       position being a cell center
% inputs:
%   env: environment variable (usually named env) of the pipeline script
%       which should contain:
%           env.height, env.width, env.depth, env.vol and env.volmask
%   opt: option variable (usually named opt) of the pipeline script
%        which should contain:
%           minrad: min radius of cell
%           maxrad: max radius of cell
%           thres:  threshold of a segmentation, 
%               the lower the thres, the more supervoxel will be found
%               default is [] for auto estimation
%   

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
    %slice = rescalegd(im2double(env.vol(:,:,zi)), [4/2048/2048 4/2048/2048]);  % caution: inline parameters
    slice = normvol(:,:,zi);
    
    % smooth  smooth width should be smaller than cell diameter
    gs = fspecial('gaussian', maxrad, 1);
    smslice = imfilter(slice, gs);
    
    % contrast enhancement
    %numtiles = ceil(size(slice)/10);
    %slice = adapthisteq(slice, 'numTiles', numtiles);
    %slice = adapthisteq(slice);
    
    %% detect
    radii = minrad:maxrad;
    [points, scoremap, radmap] = detect_by_ring_2d_cytosol(smslice, radii, thres, env.volmask(:,:,zi));
    scoremaps(:,:,zi) = scoremap;
    radmaps(:,:,:,zi) = radmap;
% % % %     figure; plot(sort(scoremap(:)));
%     
%     rr = slice; gg = slice; bb = slice;
%     rr(points) = 1;
%     gg(points) = 0;
%     bb(points) = 0;
%     %slice = cat(3, rr, gg, bb);
%     
%     s =  10;
%     numtiles = ceil(size(slice)/s)
%     slice_ah = adapthisteq(slice, 'numTiles', numtiles);
%     
%     %% ws
%     conn = 4;
%     cc = bwconncomp(points, conn);
%     L = watershed_meyer(slice, conn, cc);
% 
%     rr(L==0) = 1;
%     gg(L==0) = 1;
%     bb(L==0) = 1;
%     showz = cat(3, rr, gg, bb);
%     figure; imshow(showz);
    
    %showz = slice;
    %showz(L == 0) = 255;

    %% mask
    %points = points & env.volmask(:,:,zi);
    
    if sum(points(:)) == 0
        continue;
    end
    
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
    size(center,1)
    idx = remove_overlap(center, rin, scores);
    length(idx)
    center = center(idx,:);
    rad = rad(idx);
    rex = rex(idx);
    rin = rin(idx);
    
    %%
    supervoxel = [supervoxel; center zi*ones(length(rin),1) rin rex];
    
    %% for display
% % % %     cind = sub2ind(size(slice), center(:,2), center(:,1));
% % % %     tshow = show(:,:,zi);
% % % %     tshow(cind) = 4096;
% % % % %     show(:,:,zi) = tshow;
% % % %     figure; imshow(normim(tshow));
end    



