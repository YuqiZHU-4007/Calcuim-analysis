%I=rgb2gray(test);%im2uint8(vol);

imgt=Tiff([imagedir '\' listdir(ii).name],'r');
imgraw=read(imgt);
img=imgraw(heartreg{1},heartreg{2});
I=img;
%h=fspecial('gaussian',8);fI=imfilter(I,h,'replicate');  %��ѡ����3,5,8
fI=medfilt2(I,[3 3]);%��ֵ�˲�ȥ����
%fI=wiener2(I,[3 3]);
%h=fspecial('gaussian',8);fI=imfilter(fI,h,'replicate');%fI=imdilate(fI, strel('disk', 1));
figure,
subplot(1,2,1);imshow(I,[0 500]);title('original image');
subplot(1,2,2);imshow(fI,[0 500]);title('medfilt');

%���á�laplacian������:
laplacianSample = uint8(filter2(fspecial('laplacian'),fI));%'gaussian', 'sobel', 'prewitt', 'laplacian', 'log', 'average', 'unsharp', 'disk', 'motion'
figure,imshow(fI+laplacianSample);
BW1=edge(fI,'Canny',0.2);mask = imfill(BW1, 'holes');
mask = imfillholes(mask, 8);figure,imshow(mask);
BW = im2bw(fI+laplacianSample,0.15);BW1 = bwperim(BW,8);mask = imfill(BW1, 'holes');
figure,subplot(1,2,1);imshow(mask);
subplot(1,2,2);imshow(I);
tiffwrite(im2uint8(BW1),'C:\Users\pc\Desktop\zhf\after_pro.tif');%����ͼƬ


%   ����ͼ����ǿ
im=I;im=im2double(im);
%---ԭʼͼ��
figure,subplot(2,4,1);
imshow(im,[0 500]);
title('1��ԭʼͼ��');
%---Ϊͼ2��ʹ��ģ��Ϊ[-1,-1,-1;-1,8,-1;-1,-1,-1]���˲�����ԭͼ�����������˹������Ϊ�˱�����ʾ����ͼ������˱궨����һ���ȶ�ͼ����г��������˲���
subplot(2,4,2);
h=fspecial('laplacian',0.1);%h =[-1,-1,-1;-1,8,-1;-1,-1,-1];
im1 =imfilter(im,h);
imshow(im1,[0 max(max(im1))]);
title('2��������˹������ͼ��');
%---ͼ3������ʹ�õ�ģ�����ϣ��ó���c=1���򵥵Ľ�ԭͼ��ͼ2��ӾͿ��Եõ�һ�������񻯹���ͼ��
subplot(2,4,3);
im2=im-im1;
imshow(im2,[0 max(max(im2))])
title('3:1ͼ��2ͼ��Ӻ�ͼ��');
%---ͼ4����ԭͼ������Sobel�ݶȲ���������gxΪ[-1,-2,-1;0,0,0;1,2,1],������gyΪ[-1,0,1;-2,0,2;-1,0,1]��ģ�塣
subplot(2,4,4);
h=fspecial('prewitt');%hx=[-1,-2,-1;0,0,0;1,2,1]; %����sobel��ֱ�ݶ�ģ��hy=[-1,0,1;-2,0,2;-1,0,1]; %����sobelˮƽ�ݶ�ģ��
gradx=filter2(h,im1-im2,'same');
gradx=abs(gradx); %����ͼ���sobel��ֱ�ݶ�
grady=filter2(h',im1-im2,'same');
grady=abs(grady); %����ͼ���sobelˮƽ�ݶ�
im3=gradx+grady;  %�õ�ͼ���sobel�ݶ�
imshow(im3,[0 max(max(im3))]);
title('4:1ͼsobel�ݶȴ����ͼ��');
% ---ͼ5��ʹ�ô�СΪ5*5��һ����ֵ�˲����õ�ƽ�����Sobel�ݶ�ͼ��
subplot(2,4,5);
h1 = fspecial('average',3);im4 =imfilter(im3,h1);
%im4 =medfilt2(im3,[3 3]);%
imshow(im4,[0 max(max(im4))]);
title('5:ʹ��5*5��ֵ�˲���ƽ�����sobelͼ��');
% --ͼ6����������˹ͼ�񣨼�ͼ3����ƽ������ݶ�ͼ�񣨼�ͼ5�����е�ˡ�
subplot(2,4,6);
% im5=immultiply(im2,im4);
im5=im2.*im4;
imshow(im5,[0 max(max(im5))]);
title('6:3ͼ��5ͼ�����˵��ڱ�ͼ��');
% --ͼ7�����˻�ͼ�񣨼�ͼ6����ԭͼ����ӾͲ���һ����Ҫ����ͼ��
subplot(2,4,7);
im6=im+im5;
imshow(im6,[0 max(max(im6))]);
title('7:1ͼ��6ͼ��͵õ�����ͼ��');
% --ͼ8������ϣ����չ�Ҷȷ�Χ����ͼ7�������ʱ任����r=0.5��c=1��Ȼ�󼴿ɶ�ͼ��������ʱ任
subplot(2,4,8);
gamma=0.8;
c=1;
im7=abs(c.*im6.^gamma);
imshow(im7,[0 max(max(im7))]);
title('8:ͼ7�������ʱ任�������ͼ��');

fI=medfilt2(im7,[4 4]);figure,imshow(fI);%��ֵ�˲�ȥ����
BW1=edge(fI,'Canny',0.4);mask = imfill(BW1, 'holes');
mask = imfillholes(mask, 8);figure,imshow(mask);
BW = im2bw(fI,0.4);BW1 = bwperim(BW,8);mask = imfill(BW1, 'holes');figure,imshow(mask);

