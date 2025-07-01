clc;clear all;
[imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
savepath=uigetdir(imgpath,'savedir');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\';
T=25;
fmt='tif';
Compression='lzw';%deflate

tic;
zpath={};
for jj=1:T
    zpath{jj,1} = checkpath([savepath '/raw tiff/z' num2str(jj, '%02d')]);
end   
toc;

imio = dcimgreader([imgpath imgname], T, 1);  %, opt.t0); % auto estimate t
t=imio.nt;
tic;
for ii=1:t
    frame=imio.readvolume(ii);
    for  jj=1:T
         imwrite(frame(:,:,jj), [zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt],'Compression',Compression);
%          frame_lzw=imread([zpath{jj,1} '/' num2str(ii, '%04d') '.' fmt]);
%          com(ii,jj)=isequal(frame(:,:,jj), frame_lzw);
    end
end
toc;
imio.close();
toc;

% %imwrite(frame(:,:,jj), [zpath{jj,1} '/' 'packbits' '.tif'],'Compression','packbits');      
% imwrite(frame(:,:,jj), [zpath{jj,1} '/' 'deflate' '.tif'],'Compression','deflate');    
% imwrite(frame(:,:,jj), [zpath{jj,1} '/' 'lzw' '.tif'],'Compression','lzw');  