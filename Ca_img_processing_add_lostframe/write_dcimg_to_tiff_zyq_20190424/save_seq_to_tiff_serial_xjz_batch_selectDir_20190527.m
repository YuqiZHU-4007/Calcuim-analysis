%% save dcimg to lzw
clear all;

fmt='tif';
Compression='lzw';%deflate
file_str='seq';

filedir=uigetdir('path of seq folder');
% [imgname, dirpath] = uigetfile('*.seq');
imgpathes=scanDir_types(filedir,file_str);
% imgpath = [dirpath imgname];
% save_path = checkpath([dirpath '/' imgname]);
nbatch = length(imgpathes);

for batchi=1:nbatch
    imgpath=imgpathes{batchi};
    savepath=imgpath(1:end-4);
    savepath = checkpath(savepath);
    [frame1, nframes] = readtailseq(imgpath,1);
    for ii = 1:nframes
        frameid1 = ii;
        [frame1, nframes] = readtailseq(imgpath, frameid1);
        imwrite(frame1,[savepath '/' num2str(ii, '%05d') '.' fmt],'Compression',Compression)
    end
end

%% save tif to lzw
clear all;

fmt='tif';
Compression='lzw';%deflate
file_str='tif';

filedir=uigetdir('path of tif folder');
% [imgname, dirpath] = uigetfile('*.seq');
imgpathes = dir(filedir);
imgpathes = imgpathes(3:end);
% imgpathes=scanDir_types(filedir,file_str);
% imgpath = [dirpath imgname];
% save_path = checkpath([dirpath '/' imgname]);
nbatch = length(imgpathes);

for batchi=1:(nbatch-1)
    imgname = imgpathes(batchi).name;
    imgpath = imgpathes(batchi).folder;
    savepath = [imgpath, '\', imgname];
%     savepath = checkpath(savepath);
%     [frame1, nframes] = readtailseq(imgpath,1);
    frameinfo = dir(savepath);
    frameinfo = frameinfo(3:end);
    nframes = length(frameinfo);
    for ii = 1:nframes
        frameid1 = ii;
        frame1 = imread([frameinfo(ii).folder, '\',frameinfo(ii).name]);
        imwrite(frame1,[savepath '/' num2str(ii, '%04d') '.' fmt],'Compression',Compression)
    end
end