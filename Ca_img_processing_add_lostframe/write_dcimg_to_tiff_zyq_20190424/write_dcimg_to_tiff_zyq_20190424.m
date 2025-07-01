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




