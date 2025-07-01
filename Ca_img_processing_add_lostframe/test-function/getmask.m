function [volmask,mask]=getmask(slice,loc_nucleus,opt,env)
if loc_nucleus==1
    slice1 = rescalegd(slice, [0, 1/100]);  %��ɻҶ�ͼ��[0, 1/100]ȡֵ��Χ��x%����Ϊ�ڣ�y%����Ϊ��
    pslice = imfilter(slice1, fspecial('gaussian', opt.maxrad*2, opt.maxrad/2));%fspecial��Create predefined 2-D filter
    %��˹��ͨ�˲��������ø�˹�˵�һ��2ά�ľ�����ӣ�����ͼ��ģ������ȥ��ϸ�ں�������
    level = multithresh(pslice, 2);%Calculate a single threshold value for the image.imquantize to convert image A into an image with N + 1 discrete levels
    level = level(1);
    mask = im2bw(pslice, level);
else
    slice1 = rescalegd(slice, [0, 1/100]);
    pslice = imfilter(slice1, fspecial('gaussian', opt.maxrad*2, opt.maxrad/2));
    pslice = imdilate(pslice, strel('disk', opt.minrad));
    level = multithresh(pslice, 2);
    level = level(1);
    mask = im2bw(pslice, level);
    mask = imclose(mask, strel('disk', 10));
    mask = imfill(mask, 'holes');
    mask = imfillholes(mask, [256, 256], 8);
    mask = imdilate(mask, strel('disk', env.width/16+1));%%����mask
end

if loc_nucleus==1
    pslice = imfilter(slice, fspecial('gaussian', opt.minrad*2, opt.minrad/2));
    smvol= pslice;
    try
        level = multithresh(smvol(:), 3);%���ڴ����飬if��ǰ�г���ֱ����ת��catch��ִ��
    catch
        level = multithresh(smvol(:), 2);
    end
    level = level(1);
    volmask = false(size(slice));
    bw = smvol > level;
    bw = bw & mask;
    volmask = bw;
else
    multithreshcount = 2;
    smvol = zeros(size(slice));
    slice1 = im2double(slice);
    slice1 = normim(slice1);
    slice1 = imdilate(slice1, strel('disk', opt.minrad));
    smvol = slice1;
    level = multithresh(smvol(:), multithreshcount);
    
    volmask = false(size(slice));
    bw = im2bw(smvol, level(1));
    bw = imfillholes(bw, [256, 256], 8);
    bw = bw & mask;
    volmask= bw;
end