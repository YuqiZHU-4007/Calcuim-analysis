function [im]=heartbeat_test1(imagedi,listdir,ii,heartreg)
%diff;logic and;swell;腐蚀

imgt=Tiff([imagedir '\' listdir(ii).name],'r');
imgraw=read(imgt);
img=imgraw(heartreg{1},heartreg{2});
%img = imfilter(img, fspecial('gaussian', 2, 2));%fspecial：Create predefined 2-D filter
figure,subplot(3,3,1),imshow(img,[0 500]);

imgt2=Tiff([imagedir '\' listdir(ii+20).name],'r');
imgraw2=read(imgt2);
img2=imgraw2(heartreg{1},heartreg{2});
subplot(3,3,2),imshow(img2,[0 500]);

imgt3=Tiff([imagedir '\' listdir(ii+400).name],'r');
imgraw3=read(imgt3);
img3=imgraw3(heartreg{1},heartreg{2});
subplot(3,3,3),imshow(img3,[0 500]);

im1=img2-img;im1=im2double(im1);im1 = rescale(im1);  %变成灰度图，[0, 1/100]取值范围，x%以下为黑，y%以上为白
thresh = graythresh(im1);     %自动确定二值化阈值
im1= im2bw(im1,thresh);       %对图像二值化
subplot(3,3,4),imshow(im1);

im2=img3-img;im2=im2double(im2);im2 = rescale(im2);  %变成灰度图，[0, 1/100]取值范围，x%以下为黑，y%以上为白
thresh = graythresh(im2);     %自动确定二值化阈值
im2= im2bw(im2,thresh);       %对图像二值化
subplot(3,3,5),imshow(im2);

im=im1 & im2;
subplot(3,3,6),imshow(im)

im = imfilter(im, fspecial('gaussian',2,2));subplot(3,3,7);imshow(im)
im = imdilate(im, strel('disk', 6));subplot(3,3,8);imshow(im)
im = imerode(im,strel('disk', ceil(size(im,1)/40)));subplot(3,3,9);imshow(im)