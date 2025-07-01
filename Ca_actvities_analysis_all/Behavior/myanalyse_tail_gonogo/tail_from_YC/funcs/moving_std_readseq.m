function [stdwin, stdwintail] = moving_std(imgpath, winsize, tailind)
% [stdwin, stdwinmask] = moving_std(fpaths, winsize, masks)
%

hwin    = ceil((winsize-1)/2);

% frame1 = imread(fpaths{1});
frame1=imread(fullfile(imgpath,num2str(1,'%07d')),'tif');nframes= length(dir(imgpath))-2;
[sizey, sizex] = size(frame1);
npix = sizey * sizex;
% nframes = length(imgpath);

tstack = zeros(sizey, sizex, winsize);
for fi = 1:winsize
    frame=imread(fullfile(imgpath,num2str(fi,'%07d')),'tif');nframes= length(dir(imgpath))-2;
    %frame = medfilt2(frame, [3 3]);
    frame = double(frame);
    % normalize ?
    frame = frame ./ sum(frame(:));  %%removing illumination difference across time caused by lightsheet scanning
    tstack(:,:,fi) = frame;
end

s_tstack = sum(tstack, 3);
s2_tstack = sum(tstack.^2, 3);
leftptr   = 1;

stdwin = zeros(nframes, 1);
stdwintail = zeros(nframes, 1);
%stdwintail2 = zeros(nframes, 1);

for fi = hwin+1:nframes-hwin
    if mod(fi, 10000) == 1, disp(fi); end  % print
    %     sd = std(tstack, [], 3);  % todo: addative comput
    
    sd = sqrt(winsize/(winsize-1)*(s2_tstack/winsize - ((s_tstack/winsize).^2)));
    stdwin(fi) = mean(sd(:));
    stdwintail(fi) = mean(sd(tailind));
    %stdwintail2(fi) = mean(sd(tailind2));
    
    %     frame = imread(imgpath{min([nframes, fi+hwin+1])});
    frame=imread(fullfile(imgpath,num2str(min([nframes, fi+hwin+1]),'%07d')),'tif');nframes= length(dir(imgpath))-2;
    %frame = medfilt2(frame, [3 3]);
    frame = double(frame);
    frame = frame ./ sum(frame(:));
    s_tstack = s_tstack - tstack(:,:,leftptr) + frame;
    s2_tstack = s2_tstack - tstack(:,:,leftptr).^2 + frame.^2;
    tstack(:,:,leftptr) = frame;
    
    leftptr = mod(leftptr, winsize) + 1;
end

stdwin = abs(stdwin);
stdwintail = abs(stdwintail);



