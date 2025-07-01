clc;clear all;
T=25; % delete the last frame of each volume
[imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
h5path=[[imgpath imgname],'_h5path','.h5'];
compress=4;
if exist(h5path)
    delete(h5path);
end
imio = dcimgreader([imgpath imgname], T, 1);  %, opt.t0); % auto estimate t
y=imio.readframe(1,1);Datatype=class(y);datasize=size(y);
t=imio.nt;
tic;
for ii=1:t
    frame=imio.readvolume(ii);frame(:,:,T)=[];
    h5create(h5path,['/' num2str(ii, '%04d')],[datasize(1) datasize(2) T],'Datatype',Datatype,...
        'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
    h5write(h5path,['/' num2str(ii, '%04d')],frame);
end
toc;
imio.close();
toc;
