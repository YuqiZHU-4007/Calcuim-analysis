function [im]=heartbeat_test1(imagedi,listdir,ii,heartreg)
%diff;logic and;swell;��ʴ

imgt=Tiff([imagedir '\' listdir(ii).name],'r');
imgraw=read(imgt);
img=imgraw(heartreg{1},heartreg{2});
%img = imfilter(img, fspecial('gaussian', 2, 2));%fspecial��Create predefined 2-D filter
figure,subplot(3,3,1),imshow(img,[0 500]);

imgt2=Tiff([imagedir '\' listdir(ii+20).name],'r');
imgraw2=read(imgt2);
img2=imgraw2(heartreg{1},heartreg{2});
subplot(3,3,2),imshow(img2,[0 500]);

imgt3=Tiff([imagedir '\' listdir(ii+400).name],'r');
imgraw3=read(imgt3);
img3=imgraw3(heartreg{1},heartreg{2});
subplot(3,3,3),imshow(img3,[0 500]);

im1=img2-img;im1=im2double(im1);im1 = rescale(im1);  %��ɻҶ�ͼ��[0, 1/100]ȡֵ��Χ��x%����Ϊ�ڣ�y%����Ϊ��
thresh = graythresh(im1);     %�Զ�ȷ����ֵ����ֵ
im1= im2bw(im1,thresh);       %��ͼ���ֵ��
subplot(3,3,4),imshow(im1);

im2=img3-img;im2=im2double(im2);im2 = rescale(im2);  %��ɻҶ�ͼ��[0, 1/100]ȡֵ��Χ��x%����Ϊ�ڣ�y%����Ϊ��
thresh = graythresh(im2);     %�Զ�ȷ����ֵ����ֵ
im2= im2bw(im2,thresh);       %��ͼ���ֵ��
subplot(3,3,5),imshow(im2);

im=im1 & im2;
subplot(3,3,6),imshow(im)

im = imfilter(im, fspecial('gaussian',2,2));subplot(3,3,7);imshow(im)
im = imdilate(im, strel('disk', 6));subplot(3,3,8);imshow(im)
im = imerode(im,strel('disk', ceil(size(im,1)/40)));subplot(3,3,9);imshow(im)