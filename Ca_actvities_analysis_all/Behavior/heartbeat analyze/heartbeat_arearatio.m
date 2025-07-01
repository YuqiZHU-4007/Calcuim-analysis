%heartbeat analyze
[imagedir]=uigetdir('G:\','behavior images floder path');
listdir = dir([imagedir '/*.tif']);
heartreg={[194:340],[1064:1261]};%[x_left,x_right],[y_up,y_bottom]
 if isempty(listdir)
     error('no images in this floder')
 end
     for ii = 1:length(listdir)
           filename = listdir(ii).name;
           %img=tiffread([imagedir '\' filename]);
           imgt=Tiff([imagedir '\' listdir(ii).name],'r');
           imgraw=read(imgt);
           img=imgraw(heartreg{1},heartreg{2});
           %img = imfilter(img, fspecial('gaussian', 2, 2));%fspecial：Create predefined 2-D filter
           figure,subplot(3,3,1),imshow(img,[0 500]);
           
           %imwrite(im,'C:\Users\pc\Desktop\im1.tiff');
           
           %图像增强
           A=fspecial('average',n); %生成系统预定义的3X3滤波器
           n=2;
           modelx=[-1 0;0 1];
           modely=[0 -1;1 0];
           Iend=img;
           Idouble=double(I1);
           for i=1:size(Idouble,1)-5+n
               for j=1:size(Idouble,2)-5+n
                   area=Idouble(i:i+n-1,j:j+n-1);
                   x=area.*modelx;
                   y=area.*modely;
                   grad=max(abs(sum(x(:))),abs((sum(y(:)))));
                   Iend(i,j)=grad;
               end
           end
           figure,imshow(Iend);
           
     end
