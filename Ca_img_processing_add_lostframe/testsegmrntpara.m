function testsegmrntpara(env,opt)
opttest=opt;
opttest.minrad=floor(5/2)+1;
opttest.maxrad=floor(13/2)+1;
opttest.thres=[];
str='5_13';
minrad = opttest.minrad;
maxrad = opttest.maxrad;
thres = opttest.thres;
ti=1;
zi=17;
n=env.depth;
slice = max(im2double(env.vol), [], 3);
slice = rescalegd(slice, [0, 1/100]);
pslice = imfilter(slice, fspecial('gaussian', opttest.maxrad*2, opttest.maxrad/2));
level = multithresh(pslice, 2);
level = level(1);
mask = im2bw(pslice, level);

smvol = zeros(env.height, env.width, env.depth);
for zi = 1:n%envtest.depth
    slice = im2double(env.vol(:,:,zi));
    pslice = imfilter(slice, fspecial('gaussian', opttest.minrad*2, opttest.minrad/2));
    smvol(:,:,zi) = pslice;
end

try
    level = multithresh(smvol(:), 3);
catch
    level = multithresh(smvol(:), 2);
end
level = level(1);%%level和层数有关！！！vlomaask的效果和层数有关.min/manrad并不影响level
% level = 0.002;

volmask = false(size(env.vol));
for zi = 1:n%envtest.depth
    slice = smvol(:,:,zi);
    %     bw = im2bw(pslice, level(1));
    bw = slice > level;
    bw = bw & mask;
    volmask(:,:,zi) = bw;
%     figure,subplot(1,3,1),imshow(env.vol(:,:,zi),'DisplayRange',[0 500]);
%     subplot(1,3,2),imshow(volmask(:,:,zi));
%     subplot(1,3,3),imshow(env.volmask(:,:,zi));%%
    %saveastiff(im2uint8(volmask(:,:,zi)), ['Z:\data\fear conditioning\test' '/volmasktest' num2str(zi),'.tiff']);
end

figure,subplot(1,3,1),imshow(env.vol(:,:,1),'DisplayRange',[0 500]);
subplot(1,3,2),imshow(volmask(:,:,1));
subplot(1,3,3),imshow(env.volmask(:,:,1));

supervoxel = [];
scoremaps = zeros(env.height, env.width, env.depth);
radmaps   = zeros(env.height, env.width, 2, env.depth);

normvol = normim(env.vol, [0 0]);  % normalize together
for zi = 1:n%env.depth
    zi
    %slice = rescalegd(im2double(env.vol(:,:,zi)), [4/2048/2048 4/2048/2048]);  % caution: inline parameters
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
    size(center,1)
    idx = remove_overlap(center, rin, scores);
    length(idx)
    center = center(idx,:);
    rad = rad(idx);
    rex = rex(idx);
    rin = rin(idx);
    
    %%
    supervoxel = [supervoxel; center zi*ones(length(rin),1) rin rex];
end

figure;imshow(env.vol(:,:,zi),'DisplayRange',[100 500])
hold on
calcium_13 =supervoxel(find(supervoxel(:,3)==zi),:);
for i=1:size(calcium_13,1)
    rectangle('position',[calcium_13(i,1)-calcium_13(i,5) calcium_13(i,2)-calcium_13(i,5) calcium_13(i,5)*2 calcium_13(i,5)*2],'curvature',[1 1],'EdgeColor','b'); hold on
    rectangle('position',[calcium_13(i,1)-calcium_13(i,4) calcium_13(i,2)-calcium_13(i,4) calcium_13(i,4)*2 calcium_13(i,4)*2],'curvature',[1 1],'EdgeColor','r'); hold on
end
%saveas(gcf,'Z:\data\fear conditioning\test/volmasktest1.tiff');
h=getimage(gcf);
imwrite(h,['Z:\data\fear conditioning\test/volmasktest' str '.tiff'])

a(:,:,1)=env.vol(:,:,zi)';a(:,:,2)=env.vol(:,:,zi)';a(:,:,3)=env.vol(:,:,zi)';
for i=1:size(calcium_13,1)
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)+calcium_13(i,4),1)=1;
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4),1)=1;
    a(calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),1)=1;
    a(calcium_13(i,1)-calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),1)=1;
    
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)+calcium_13(i,4),2)=0;
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4),2)=0;
    a(calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),2)=0;
    a(calcium_13(i,1)-calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),2)=0;
    
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)+calcium_13(i,4),3)=0;
    a(calcium_13(i,1)-calcium_13(i,4):calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4),3)=0;
    a(calcium_13(i,1)+calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),3)=0;
    a(calcium_13(i,1)-calcium_13(i,4),calcium_13(i,2)-calcium_13(i,4):calcium_13(i,2)+calcium_13(i,4),3)=0;
end
figure;imshow(a);imwrite(a',['Z:\data\fear conditioning\test/volmasktest' num2str(1), str '.tiff'])
%saveastiff(im2uint8(a), ['Z:\data\fear conditioning\test' '/volmasktest' num2str(1),'.tiff']);

%scatter(calcium_13(:,1),calcium_13(:,2),calcium_13(:,5).^2,'r');hold on
%scatter(calcium_13(:,1),calcium_13(:,2),calcium_13(:,4).^2,'b');hold on
% for i=1:size(calcium_13,5)
% scatter(calcium_13(:,1),calcium_13(:,2),calcium_13(:,5),'r');hold on
% end[2 4 2 2]
