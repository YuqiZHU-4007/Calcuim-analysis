
% todo: make all the settings above
% copy out moved images

%% load images
ext = '.tif';

nframes = 0;
outdir = 'Z:\light-sheet\jturn\20151205_Beha_fish_2';
for ii = 1:10
    dirpath = [outdir '/20151205_Beha_fish_2_trial_' num2str(ii)];
    nframes = nframes + length(dir(dirpath)) - 2;
end
fpaths = cell(nframes, 1);
k = 1;
for ii = 1:10
    dirpath = [outdir '/20151205_Beha_fish_2_trial_' num2str(ii)];
    for jj = 1:3000
        fpaths{k} = [dirpath '/RecordedImage_GO-5000M-USB__' num2str(jj-1, '%.3d') ext];
        k = k + 1;
    end
end

%% options
savepath = checkpath('Z:/light-sheet/jturn/test_detect/');
npoints = 9;   % 16 segments
% % % % params for find_tail
% % % controlpoint = 8;  % 
% % % 
% % % params.prolong_factor = 1.5;
% % % params.turn_scope = 40;
% % % params.turn_step = 2;
% % % params.max_width = 11;
% % % params.min_width = 3;

winsize = 301;
boxwidth = 32;
stdthres = 2.4e-7;


frame1 = 255 - imread(fpaths{2});
[sizey, sizex] = size(frame1);
npix = sizey * sizex;
% figure; imshow(frame1);

%% inspect
% % % ind = 11870:12450;
% % % % suspect region
% % % ind = 45320:45850;    % moved
% % % ind = 126200:126900;  % swinged
% % % ind = 191900:192500;  % moved
% % % ind = 194800:195800;  % swinged
% % % ind = 82702:86300;
% % % figure; 
% % % hdim = imshow(frame1);
% % % for fi = ind(:)'
% % %     frame = imread(fpaths{fi});
% % %     set(hdim, 'cdata', frame);
% % %     drawnow;
% % %     xlabel(fi);
% % % end

%% interactive
% params
%enhframe1 = mat2gray(frame1, vlim);
%enhframe1 = adapthisteq(frame1);
selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);


%% movement detection
hwin    = ceil((winsize-1)/2);
tailind = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth);
tailind2 = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth / 2);
% caution: inline params

tstack = zeros(sizey, sizex, winsize);
for fi = 1:winsize
    frame = 255 - imread(fpaths{fi});
    frame = medfilt2(frame, [3 3]);
    frame = double(frame);
    % normalize ?
    frame = frame ./ sum(frame(:));
    tstack(:,:,fi) = frame;
end

s_tstack = sum(tstack, 3);
s2_tstack = sum(tstack.^2, 3);
leftptr   = 1;

tic;
stdwin = zeros(nframes, 1);
stdwintail = zeros(nframes, 1);
stdwintail2 = zeros(nframes, 1);

for fi = hwin+1:nframes-hwin
    if mod(fi, 10000) == 4000, disp(fi); toc; end  % print 
%     sd = std(tstack, [], 3);  % todo: addative comput
    
    sd = sqrt(winsize/(winsize-1)*(s2_tstack/winsize - ((s_tstack/winsize).^2)));
    stdwin(fi) = mean(sd(:));
    stdwintail(fi) = mean(sd(tailind));
    stdwintail2(fi) = mean(sd(tailind2));
    
    frame = 255 - imread(fpaths{min([nframes, fi+hwin+1])});
    frame = medfilt2(frame, [3 3]);
    frame = double(frame);
    frame = frame ./ sum(frame(:));
    s_tstack = s_tstack - tstack(:,:,leftptr) + frame;
    s2_tstack = s2_tstack - tstack(:,:,leftptr).^2 + frame.^2;
    tstack(:,:,leftptr) = frame;
    
    leftptr = mod(leftptr, winsize) + 1;
end
duration_computestd = toc;

figure; plot(abs(stdwintail))

% stdthres = mean(abs(stdwintail))

moving = (stdwintail > stdthres);
movlabel = bwlabel(imclose(moving(:)', ones(1, winsize)));
movind = find(moving);

save([savepath '/movement_detection.mat'], 'stdwintail', 'stdwin', 'stdwintail2', 'tailind', 'tailind2', 'selx', 'sely', 'frame1', 'movind', 'movlabel');


