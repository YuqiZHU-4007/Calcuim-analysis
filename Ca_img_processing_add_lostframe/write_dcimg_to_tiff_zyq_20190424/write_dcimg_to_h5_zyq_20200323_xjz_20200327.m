clc;clear all;
file_str='dcimg';

fliedir=uigetdir('path of dcimg folder');
[txtname,txtpath]=uigetfile('C:\standard brain\*.txt','txt of elapsedTime');
imgpathes=scanDir_types(fliedir,file_str);

fid=fopen([txtpath txtname],'wt');
nbatch = length(imgpathes);
T=[25,30];%8093;9937

process.T_lessframe=[];
process.z_lessframe=[];
process.T_moreframe=[];
process.z_moreframe=[];
for batchi=1:nbatch
    imgpath=imgpathes{batchi}
    T1=T(batchi); % delete the last frame of each volume
%     [imgname,imgpath]=uigetfile('C:\standard brain\20180901_31dpf\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
    h5path=[imgpath(1:end-6),'_h5path'];
    mkdir(h5path);
    compress=4;

    imio = dcimgreader(imgpath, T1, 1,process);  %, opt.t0); % auto estimate t
    y=imio.readframe(1,1);Datatype=class(y);datasize=size(y);
    t=imio.nt;
    tic;
    for ii=1:t
        frame=imio.readvolume(ii);frame(:,:,T1)=[];
        if exist(fullfile(h5path,['/' num2str(ii, '%06d'),'.h5']))
            delete(fullfile(h5path,['/' num2str(ii, '%06d'),'.h5']));
        end
        h5create(fullfile(h5path,['/' num2str(ii, '%06d'),'.h5']),'/default',[datasize(1) datasize(2) size(frame,3)],'Datatype',Datatype,...
            'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
        
        h5write(fullfile(h5path,['/' num2str(ii, '%06d'),'.h5']),'/default',frame);
        %         h5create(h5path,['/' num2str(ii, '%04d')],[datasize(1) datasize(2) (T1-1)],'Datatype',Datatype,...
        %             'chunksize',[datasize(1) datasize(2) 1],'Deflate',compress);
        %         h5write(h5path,['/' num2str(ii, '%04d')],frame);
    end
    toc;
    imio.close();
    toc;
end
