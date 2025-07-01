
% todo: make all the settings above
% copy out moved images

%% load images
[fpaths, nframes] = ask_for_images;  % caution: image formats/names dependent

%% options
dirpath = fileparts(fpaths{1});
savepath = [dirpath '/../test_tail/']; %% /.. stands for parental folder
% savepath = checkpath('Z:/data/20150625/rec5338806_sig/tail_new/');
npoints = 17;   % 16 segments
% % % % params for find_tail
% % % controlpoint = 8;  % 
% % % 
% % % params.prolong_factor = 1.5;
% % % params.turn_scope = 40;
% % % params.turn_step = 2;
% % % params.max_width = 11;
% % % params.min_width = 3;

winsize = 301; %%temporal window for movement detection
boxwidth = 50; %%the spatial zone for tail searching
stdthres = 0.6e-7;

%

frame1 = imread(fpaths{2});
[sizey, sizex] = size(frame1);
npix = sizey * sizex;


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

%% interactive%%getting tail axis poisition by drawing
%enhframe1 = mat2gray(frame1, vlim);
%enhframe1 = adapthisteq(frame1);
selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);

%% movement detection  %% for indices of all points in the searching zone
hwin    = ceil((winsize-1)/2);
tailind = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth);
tailind2 = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth / 4);
% caution: inline params

tstack = zeros(sizey, sizex, winsize);
for fi = 1:winsize
    frame = imread(fpaths{fi});
    frame = medfilt2(frame, [3 3]);
    frame = double(frame);
    % normalize ?
    frame = frame ./ sum(frame(:)); %%removing illumination difference across time caused by lightsheet scanning
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
    
    frame = imread(fpaths{min([nframes, fi+hwin+1])});
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

moving = (stdwintail > stdthres);
movlabel = bwlabel(imclose(moving(:)', ones(1, winsize)));
movind = find(movlabel > 0);

save([savepath '/movement_detection.mat'], 'stdwintail', 'stdwin', 'stdwintail2', 'tailind', 'tailind2', 'selx', 'sely', 'frame1', 'movind', 'movlabel');


%{
%% movedframes
nframes = length(stdwintail);
nvol = size(dfdf, 2);  % needs calc info

% caution: inline params
tp_calc = 0.4 * 1000;
fs_calc = 1000/tp_calc;
tp_tail = 3.333;
fs_tail = 1000/tp_tail;

fs_tail = 299.68;
tp_tail = 1000/fs_tail;

fs_tail = 240;
tp_tail = 1000/fs_tail;

corfi = interp1(linspace(0, tp_tail*nframes, nframes), 1:nframes, linspace(0, tp_calc*nvol, nvol), 'nearest');

corwin  = floor(fs_tail / fs_calc / 2) * 2 + 1;
hcorwin = (corwin - 1) / 2;
movthres = 0.25;

moving = (stdwintail > 0.5e-7);
% movlabel = bwlabel(moving(:));

movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution
movind = find(moving);

imat = veccolon(corfi - hcorwin, 1, corwin);
imat(imat<1) = 1;
imat(imat>length(moving)) = length(moving);

movets = sum(movlabel(imat)>0, 2)/corwin > movthres;
movets = movets';

figure; imagesc(movets)

%%% 
% % % sig = dfdf(88810,:);
% % % movets = movets(1:nvol)';
% % % figure; plot(sig);
% % % hold on;
% % % plot(movets/50, 'r');
% % % 

%% inspect
figure; 
hdim = imshow(frame1);
for fi = find(movlabel==3)
    frame = imread(fpaths{fi});
    set(hdim, 'cdata', frame);
    drawnow;
    xlabel(fi);
    ylabel(movlabel(fi));
end

save movement_detection_062513.mat movets movind moving movlabel stdwintail stdwin stdwintail2 tailind tailind2 selx sely frame1 ...
    npoints theta0
% caution: stationary save path

%% copy out moved files
cppath = checkpath('D:\vc\data\20150625\fish 02');
% cppath = checkpath([imdir '/../movedd']);
cppath = checkpath('Z:\lightsheet\fish 02\fish 02');
for ii = 20000:length(movind)
    ii
    fi = movind(ii);
    [~, fname, ext] = fileparts(fpaths{fi});
    fname = [fname ext];
    if exist(fullfile(cppath, fname)) == 0
%         system(['copy "' fullfile(fpaths{fi}) '" "' fullfile(cppath) '"']);
        copyfile(fpaths{fi}, cppath);
%         movefile(fpaths{fi}, cppath);
%         frame = imread(fpaths{fi});
%         imwrite(frame, fullfile(cppath, fname));  % can't believe it's faster.
    end
        %    clear functions
end

fid = fopen('./movind.txt', 'w');
fprintf(fid, '%d\n', movind);
fclose(fid);
%}