% if ~libisloaded('tmcamcon')
%     loadlibrary('tmcamcon.dll', 'tmcamcon.h');
% end
clc;clear all;
[imgname,imgpath]=uigetfile('E:\A_Data_lightsheet\standard brain\*.dcimg','dcimg');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
savepath=uigetdir(imgpath,'savedir');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\';
T=25;
poolsize=4;


poolobj = parpool(min(poolsize,T));%gcp('nocreate');
poolsize = poolobj.NumWorkers
blocksize = ceil(T / poolsize); 
tic;
spmd(poolsize)
    imio = dcimgreader([imgpath imgname], T, 1);  %, opt.t0); % auto estimate t
    %imio.set_subsamplet(1:imio.nframes);
    zi_in_this_thread = (labindex-1)*blocksize+1:min([labindex*blocksize T]);
    for ii = 1:length(zi_in_this_thread)
        %frame=imio.readframe( ti, zi);
        frame=imio.readslice(zi_in_this_thread(ii));
        zpath = checkpath([savepath '/raw tiff/z' num2str(zi_in_this_thread(ii), '%02d')]);
        %imwrite(uint16(frame), [zpath '/' num2str(i, '%02d') '.tif']);
        seqwrite(frame, zpath, 'tif')
    end
    imio.close();
end
toc;
delete(poolobj)

%tif to tif
% end
clc;clear all;
[imgpath]=uigetdir('H:\品系管理\Huc_h2b_g7f\20200805\fish1_light-sheet\10x_z vol\rec1119557024\');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\rec186952424.dcimg';
savepath=uigetdir(imgpath,'savedir');%'E:\A_Data_lightsheet\standard brain\20180901_31dpf\';
listdir=dir(imgpath);listdir=listdir(3:end);
T=25;
tic;
for ii=1:length(listdir)
    frame=imread(fullfile(listdir(ii).folder,listdir(ii).name));
    frame_id=str2num(listdir(ii).name(15:end-4));
    z=mod(frame_id,T);
    zpath = checkpath([savepath '/raw tiff/z' num2str(z, '%02d')]);
    imwrite(frame, [zpath '/' num2str(frame_id, '%04d') '.' 'tif'],'Compression','lzw');
end
toc;


