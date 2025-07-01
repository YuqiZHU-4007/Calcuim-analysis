clc;clear all;
selectFOV=true;%datasize=filp(datasize);

[imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.seq','seq');
h5path=[[imgpath imgname],'_h5path','.h5'];
compress=4;
if exist(h5path)
    delete(h5path);
end
[y, nframes] = readtailseq([imgpath imgname],1);
%% selec FOV
datasize=size(y);
pos=ceil([datasize(2)/2 1 datasize(2)/2 datasize(1)*3/5]);
if ~selectFOV
    fov = [1 1 size(y,2) size(y,1)];
else
    figure;imshow(y,[min(y(:)) max(y(:))]);
    roi = imrect(gca,pos);wait(roi);
    fov=getPosition(roi);fov=ceil(fov);
    fov=[fov(1) max(fov(2),1) min(datasize(2)-fov(1),fov(3)) min(datasize(1)-fov(2),fov(4))];
    %delect(roi);
end
y=y(fov(2):fov(2)+fov(4),fov(1):fov(1)+fov(3));
figure;imshow(y,[min(y(:)) max(y(:))]);
Datatype=class(y);datasize=size(y);
%% save
h5create(h5path,['/' num2str(ii, '%04d')],[datasize(1) datasize(2) nframes],'Datatype',Datatype,...
    'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
tic;
for ii=1:nframes
    [frame, ~] = readtailseq([imgpath imgname], ii);
    %figure;imshow(frame,[min(frame(:)) max(frame(:))]);
    frame=frame(fov(2):fov(2)+fov(4),fov(1):fov(1)+fov(3));
    start=[1 1 ii];
    h5write(h5path,['/' num2str(ii, '%04d')],frame,start,[datasize(1),datasize(2),1]);
end
% save first raw image & FOV position
[y, nframes] = readtailseq([imgpath imgname],1);
h5create(h5path,'/fst rawimage',[datasize(1) datasize(2) 1],'Datatype',Datatype,...
    'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
h5write(h5path,'/fst rawimage',y); 
h5create(h5path,'/FOV pos',[size(fov,1) size(fov,2) 1],'Datatype',class(fov),...
    'chunksize',[size(fov,1) size(fov,2) 1],'Deflate',compress);
h5write(h5path,'/FOV pos',fov); 
h5disp(h5path)
toc;
