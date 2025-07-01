clear all;
close all
%% set para
Dir =uigetdir('/mnt/X/data/fear_conditioning/huc','tifpath');% insert path to tiff stack here;'/mnt/W/data/fear_conditioning/huc'
disp(Dir);
h5path=[Dir,'_h5path','.h5'];
listdir=dir([Dir '/*.tif']);
y= imread([listdir(1).folder,'/',listdir(1).name]);datasize=size(y);Datatype=class(y);
T=length(listdir);
datasetname='/default';
compress=4;
%figure;imshow(y,[min(y(:)) max(y(:))]);
%% select ROI
selectFOV=false;
if selectFOV
halfFOV=false;%datasize=filp(datasize);
pos=ceil([datasize(2)/2 1 datasize(2)/2 datasize(1)*3/5]);%pos=ceil([datasize(2)/2 1 datasize(2)/2 datasize(1)/2]);
if halfFOV
    fov = pos;
else
    %figure;imshow(y,[min(y(:)) max(y(:))]);
    roi = imrect(gca,pos);wait(roi);
    fov=getPosition(roi);fov=ceil(fov);
    fov=[fov(1) max(fov(2),1) min(datasize(2)-fov(1),fov(3)) min(datasize(1)-fov(2),fov(4))];
    %delect(roi);
end
else
    fov=[1 1 datasize(2)-1 datasize(1)-1];
end
a=y(fov(2):fov(2)+fov(4),fov(1):fov(1)+fov(3));
%figure;imshow(a,[min(a(:)) max(a(:))]);
%% trans
chunksize=[fov(4) fov(3) 1];
tic;
if exist(h5path)
    delete(h5path);
end
h5create(h5path,datasetname,[fov(4) fov(3) Inf],'Datatype',Datatype,...
    'chunksize',chunksize,'Deflate',compress);
for t=1:T
    yy=imread([Dir '/',num2str(t,'%05d'), '.tif'],'PixelRegion',{[fov(2),fov(2)+fov(4)-1],[fov(1),fov(1)+fov(3)-1]});
    if t==1
        %figure;imshow(yy,[min(yy(:)) max(yy(:))]);
    end
    start=[1 1 t];
    h5write(h5path,datasetname,yy,start,[fov(4),fov(3),1]); 
end
% save first raw image & FOV position
h5create(h5path,'/fst rawimage',[datasize(1) datasize(2) 1],'Datatype',Datatype,...
    'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
h5write(h5path,'/fst rawimage',y); 
h5create(h5path,'/FOV pos',[size(fov,1) size(fov,2) 1],'Datatype',class(fov),...
    'chunksize',[size(fov,1) size(fov,2) 1],'Deflate',compress);
h5write(h5path,'/FOV pos',fov); 
h5disp(h5path)
toc;

%% compare
% ti=370;
% data_h5 = h5read(h5path,datasetname,[1 1 ti],[fov(4),fov(3),1]);
% %show_spv_GUI(data_h5);
% %figure,imshow(data_h5 )
% data_tif=imread([Dir '/',num2str(ti,'%05d'), '.tif'],'PixelRegion',{[fov(2),fov(2)+fov(4)-1],[fov(1),fov(1)+fov(3)-1]});
% %figure,imshow(data_tif)
% % show_spv_GUI(data_tif);   
% length(find((data_h5-data_tif)~=0))