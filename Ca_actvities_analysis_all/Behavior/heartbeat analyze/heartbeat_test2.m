%I=rgb2gray(test);%im2uint8(vol);

imgt=Tiff([imagedir '\' listdir(ii).name],'r');
imgraw=read(imgt);
img=imgraw(heartreg{1},heartreg{2});
I=img;
%h=fspecial('gaussian',8);fI=imfilter(I,h,'replicate');  %备选参数3,5,8
fI=medfilt2(I,[3 3]);%中值滤波去噪声
%fI=wiener2(I,[3 3]);
%h=fspecial('gaussian',8);fI=imfilter(fI,h,'replicate');%fI=imdilate(fI, strel('disk', 1));
figure,
subplot(1,2,1);imshow(I,[0 500]);title('original image');
subplot(1,2,2);imshow(fI,[0 500]);title('medfilt');

%采用’laplacian’算子:
laplacianSample = uint8(filter2(fspecial('laplacian'),fI));%'gaussian', 'sobel', 'prewitt', 'laplacian', 'log', 'average', 'unsharp', 'disk', 'motion'
figure,imshow(fI+laplacianSample);
BW1=edge(fI,'Canny',0.2);mask = imfill(BW1, 'holes');
mask = imfillholes(mask, 8);figure,imshow(mask);
BW = im2bw(fI+laplacianSample,0.15);BW1 = bwperim(BW,8);mask = imfill(BW1, 'holes');
figure,subplot(1,2,1);imshow(mask);
subplot(1,2,2);imshow(I);
tiffwrite(im2uint8(BW1),'C:\Users\pc\Desktop\zhf\after_pro.tif');%保存图片


%   骨骼图像增强
im=I;im=im2double(im);
%---原始图像
figure,subplot(2,4,1);
imshow(im,[0 500]);
title('1：原始图像');
%---为图2，使用模板为[-1,-1,-1;-1,8,-1;-1,-1,-1]的滤波器对原图像进行拉普拉斯操作，为了便于显示，对图像进行了标定，这一步先对图像进行初步的锐化滤波。
subplot(2,4,2);
h=fspecial('laplacian',0.1);%h =[-1,-1,-1;-1,8,-1;-1,-1,-1];
im1 =imfilter(im,h);
imshow(im1,[0 max(max(im1))]);
title('2：拉普拉斯操作后图像');
%---图3，由于使用的模板如上，让常数c=1，简单的将原图和图2相加就可以得到一幅经过锐化过的图像。
subplot(2,4,3);
im2=im-im1;
imshow(im2,[0 max(max(im2))])
title('3:1图和2图相加后图像');
%---图4，对原图像试用Sobel梯度操作，分量gx为[-1,-2,-1;0,0,0;1,2,1],而分量gy为[-1,0,1;-2,0,2;-1,0,1]的模板。
subplot(2,4,4);
h=fspecial('prewitt');%hx=[-1,-2,-1;0,0,0;1,2,1]; %生产sobel垂直梯度模板hy=[-1,0,1;-2,0,2;-1,0,1]; %生产sobel水平梯度模板
gradx=filter2(h,im1-im2,'same');
gradx=abs(gradx); %计算图像的sobel垂直梯度
grady=filter2(h',im1-im2,'same');
grady=abs(grady); %计算图像的sobel水平梯度
im3=gradx+grady;  %得到图像的sobel梯度
imshow(im3,[0 max(max(im3))]);
title('4:1图sobel梯度处理后图像');
% ---图5，使用大小为5*5的一个均值滤波器得到平滑后的Sobel梯度图像。
subplot(2,4,5);
h1 = fspecial('average',3);im4 =imfilter(im3,h1);
%im4 =medfilt2(im3,[3 3]);%
imshow(im4,[0 max(max(im4))]);
title('5:使用5*5均值滤波器平滑后的sobel图像');
% --图6，将拉普拉斯图像（即图3）与平滑后的梯度图像（即图5）进行点乘。
subplot(2,4,6);
% im5=immultiply(im2,im4);
im5=im2.*im4;
imshow(im5,[0 max(max(im5))]);
title('6:3图和5图相乘相乘的掩蔽图像');
% --图7，将乘积图像（即图6）与原图像相加就产生一幅需要的锐化图像。
subplot(2,4,7);
im6=im+im5;
imshow(im6,[0 max(max(im6))]);
title('7:1图和6图求和得到的锐化图像');
% --图8，我们希望扩展灰度范围，对图7进行幂率变换处理，r=0.5，c=1，然后即可对图像进行幂率变换
subplot(2,4,8);
gamma=0.8;
c=1;
im7=abs(c.*im6.^gamma);
imshow(im7,[0 max(max(im7))]);
title('8:图7进行幂率变换后的最终图像');

fI=medfilt2(im7,[4 4]);figure,imshow(fI);%中值滤波去噪声
BW1=edge(fI,'Canny',0.4);mask = imfill(BW1, 'holes');
mask = imfillholes(mask, 8);figure,imshow(mask);
BW = im2bw(fI,0.4);BW1 = bwperim(BW,8);mask = imfill(BW1, 'holes');figure,imshow(mask);

