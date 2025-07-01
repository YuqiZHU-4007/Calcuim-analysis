function [volmask, mask] = background_segmentation(env, opt)
% [volmask, mask] = background_segmentation(env, opt)

% caution: inline parameters

multithreshcount = 2;

%% union mask
% slice = mean(im2double(env.vol, 3);
% level = multithresh(slice, 20);
slice = max(im2double(env.vol), [], 3);
slice = rescalegd(slice, [0, 1/100]);  

pslice = imfilter(slice, fspecial('gaussian', opt.maxrad*2, opt.maxrad/2));
pslice = imdilate(pslice, strel('disk', opt.minrad));

level = multithresh(pslice, 5); 
level = level(2);
mask = im2bw(pslice, level);
mask = imclose(mask, strel('disk', 10));
mask = imfill(mask, 'holes');
mask = imfillholes(mask, [256, 256], 8);
% figure; imshow(pslice); figure; imshow(mask);

smvol = zeros(env.height, env.width, env.depth);
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    slice = normim(slice);
    slice = imdilate(slice, strel('disk', opt.minrad));
    smvol(:,:,zi) = slice;
end

level = multithresh(smvol(:), multithreshcount);

volmask = false(size(env.vol));
for zi = 1:env.depth
    slice = smvol(:,:,zi);
    bw = im2bw(slice, level(1));
    bw = imfillholes(bw, [256, 256], 8);
    bw = bw & mask;
    volmask(:,:,zi) = bw;
end

%% oversmooth scheme
%{
%% smoothed vol
smvol = zeros(env.height, env.width, env.depth);
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    slice = rescalegd(slice, [0 0]);  % caution: inline parameters
    
    winsize = 2 * round(sqrt(sum(mask(:))) / 4) - 1;  % caution: inline parameters
    pslice = imfilter(slice, fspecial('gaussian', winsize, (winsize - 1)/4));
    pslice = rescalegd(pslice, [0 0]);

    smvol(:,:,zi) = pslice;
end

%% threshold again
% prefer retrival
volmask = false(size(env.vol));
for zi = 1:env.depth
    slice = smvol(:,:,zi);
    level = multithresh(slice, 3);
    bw = im2bw(slice, level(1));
    bw = imerode(bw, strel('disk', opt.maxrad));  % this is for oversmooth
    bw = imfillholes(bw, [100, 100], 8);
    bw = bw & mask;
    volmask(:,:,zi) = bw;
end
% level = multithresh(smvol(:), 2);
% bw = im2bw(smvol(:), level(1));
% bw = reshape(bw, size(smvol));
% bw = bw & repmat(mask, [1, 1, env.depth]);
% volmask = bw;
%}

%{
% threshold oversmooth
% tend to under segment but even
volmask = false(size(env.vol));
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    slice = rescalegd(slice, [0 0]);  % caution: inline parameters
    
    winsize = 2 * round(sqrt(sum(mask(:))) / 4) - 1;  % caution: inline parameters
    pslice = imfilter(slice, fspecial('gaussian', winsize, (winsize - 1)/4));
    pslice = rescalegd(pslice, [0 0]);
    
    % erode
    %ft = strel('disk', opt.maxrad*2+1);  %TODO: should be estimated by gaussian
    %pslice = imerode(pslice, ft);
    
    %t = sort(pslice(~mask));
    %level = t(round(length(t)*0.9));
    level = multithresh(pslice, 2);
    bw = im2bw(pslice, level(1));
    bw = imfill(bw, 'holes');
    bw = imfillholes(bw, [100, 100], 8);  % caution: inline parameters
    volmask(:,:,zi) = bw;
end
%}

%{
% threshold oversmooth
% tend to under segment but even
volmask = false(size(vol));
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    slice = rescalegd(slice, [0 0]);  % caution: inline parameters
    
    %numtiles = ceil(size(slice)/5);
    %pslice = adapthisteq(slice, 'numTiles', numtiles);
    
    %ft = strel('disk', opt.maxrad).getnhood;
    %pslice = ordfilt2(pslice, sum(ft(:)), ft);
    
    pslice = imfilter(pslice, fspecial('gaussian', 512, 128));
    pslice = rescalegd(pslice, [0 0]);
    
    %level = max(pslice(~mask));
    level = graythresh(pslice);
    bw = im2bw(pslice, level);
    bw = imfill(bw, 'holes');
    bw = imfillholes(bw, [100, 100], 8);  % caution: inline parameters
    volmask(:,:,zi) = bw;
end



    for ii = 10:5:50
        numtiles = ceil(size(slice)/ii);
        t = adapthisteq(slice, 'numTiles', numtiles);
        figure; imshow(t);
        title(ii);
    end
%}
%{
% simple
% tend to under segment
volmask = false(size(vol));
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));
    pslice = histeq(slice);
    level = graythresh(pslice);
    bw = im2bw(pslice, level);
    volmask(:,:,zi) = bw;
end
%}
    
%{
% seg by active contour

% over segment at the dark edge
% and slow ...
% and you never know what the result would be like

volmask = false(size(vol));
for zi = 1:env.depth
    zi
    slice = im2double(vol(:,:,zi));
    closedslice = imclose(slice, strel('disk', opt.maxrad*2-1));
    bw = activecontour(closedslice, mask, 1000, 'Chan-Vese', 1);
    bw = imfillholes(bw, [100, 100], 8);
    volmask(:,:,zi) = bw;
end
%}

%{
%% active contour 2
for zi = 1:T-1
    slice = rescalegd(im2double(vol(:,:,zi)), [1/10000 1/10000]);
    % threshold
    %{
    levels = zeros(20, 1);
    for ii = 1:20
        t = multithresh(slice, ii);
        levels(ii) = t(1);
    end
    level = mean(levels);
    bw = im2bw(slice, level);
    bw = imfill(bw, 'holes');
    %}
    % ac
    level = multithresh(slice, 20);
    level = level(1);
    bw = im2bw(slice, level);
    bw = imfill(bw, 'holes');
    %bw = activecontour(slice, mask, 1000);
    bw = imfillholes(bw, [10000, 10000], 8);
    volmask(:,:,zi) = bw;
end
%}

%{
% lower than the estimated background intensities
% tend to oversegment
for zi = 1:env.depth
    slice = rescalegd(im2double(env.vol(:,:,zi)), [1/10000 1/10000]);  % caution: inline parameters
    level = mean(slice(mask));
    bw = im2bw(slice, level);
    bw = imfill(bw, 'holes');
    bw = imfillholes(bw, [100, 100], 8);
    volmask(:,:,zi) = bw;
end
%}







