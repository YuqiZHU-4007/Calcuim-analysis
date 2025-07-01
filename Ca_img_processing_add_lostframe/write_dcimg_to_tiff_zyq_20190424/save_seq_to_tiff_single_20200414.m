clear all;

fmt='tif';
Compression='lzw';%deflate

[imgname dirpath]=uigetfile('*.seq','path of seq folder');
% [imgname, dirpath] = uigetfile('*.seq');
imgpath = [dirpath imgname];
[~,name,~]=fileparts(imgpath );
save_path = checkpath([dirpath '/' name]);
tic;
[frame1, nframes] = readtailseq(imgpath,1);
nframes=36000;
for ii = 1:nframes
    frameid1 = ii;
    [frame1, ~] = readtailseq(imgpath, frameid1);
    imwrite(frame1,fullfile(save_path ,[num2str(ii, '%07d') '.' fmt]),'Compression',Compression)
end
toc