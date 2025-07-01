function [volmask, mask] = background_segmentation_nucleus(env, opt)
% [volmask, mask] = background_segmentation(env, opt)

% caution: inline parameters

multithreshcount = 2;

% check all params of multithresh

%% union mask
% slice = mean(im2double(env.vol, 3);
% level = multithresh(slice, 20);
slice = max(im2double(env.vol), [], 3);
slice = rescalegd(slice, [0, 1/100]);  %变成灰度图，[0, 1/100]取值范围，x%以下为黑，y%以上为白

pslice = imfilter(slice, fspecial('gaussian', opt.maxrad*2, opt.maxrad/2));%fspecial：Create predefined 2-D filter
 h=figure('visible','on');subplot(2,1,1);imshow(pslice);
%高斯低通滤波器是利用高斯核的一个2维的卷积算子，用于图像模糊化（去除细节和噪声）
% pslice = imdilate(pslice, strel('disk', opt.minrad));

level = multithresh(pslice, 6);%Calculate a single threshold value for the image.imquantize to convert image A into an image with N + 1 discrete levels
level = level(1);%(level(1)+level(2))/2;
mask = im2bw(pslice, level);%im2bw is not recommended. Use imbinarize instead.
% 二值化图像greater than level with the value 1 (white) and replaces all other pixels with the value 0 (black). 
% mask = imclose(mask, strel('disk', 10));
% mask = imfill(mask, 'holes');
% mask = imfillholes(mask, [256, 256], 8);
% figure; imshow(pslice); 
subplot(2,1,2);imshow(mask);%saveas(h,[opt.savepath '/mask1.tif']);

smvol = zeros(env.height, env.width, env.depth);
for zi = 1:env.depth
    slice = im2double(env.vol(:,:,zi));slice = rescalegd(slice, [0, 1/100]); 
    pslice = imfilter(slice, fspecial('gaussian', opt.minrad*2, opt.minrad/2));
    smvol(:,:,zi) = pslice;
end
% figure; hist(smvol(:))
% figure; plot(sort(smvol(:)))
try
    level = multithresh(smvol(:),6) ;%用于错误检查，if当前行出错，直接跳转到catch后执行
catch
    level = multithresh(smvol(:), 2);
end
%format long
%levell =level(2);%

volmask = false(size(env.vol));
for zi = 1:env.depth
    slice = smvol(:,:,zi);
   %bw = im2bw(pslice, level(1));
   %levell = multithresh(smvol(:,:,zi), 4);levell =(levell(2)+levell(3))/2;
   if zi>=36
       levell=level(4);
   else
       levell=(level(2)+level(3))/2;
   end
    bw = slice > levell;
    %unique(bw)
    %figure,imshow(bw,[0 1]);
    bw = bw & mask;
    volmask(:,:,zi) = bw;
end
zi=34;
figure,subplot(2,1,1),imshow(smvol(:,:,zi))
subplot(2,1,2);imshow(volmask(:,:,zi));

%{
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
%}

