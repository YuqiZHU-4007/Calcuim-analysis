close all;clear all;clc
% copy out moved images
code_dir='Z:\Behavior\myanalyse_tail_gonogo\tail_from_YC'
addpath(genpath(code_dir));
%TODO: save params

%% load images
[imgname, dirpath] = uigetfile('*.seq','','W:\data\fear conditioning_ZYQ\Fear conditioning\20220309\fish2\*.seq');
imgpath = [dirpath imgname]

% [~, name, ext] = fileparts(imgname);
savepath = checkpath([dirpath '/' imgname '_find_tail/'])

%% options for movement detections
% dirpath = fileparts(fpaths{1});
% [~, dirname] = fileparts(dirpath);
% savepath = checkpath([dirpath '/../' imgname '_find_tail/']); %% /.. stands for parental folder
% savepath = checkpath('Z:/data/20150625/rec5338806_sig/tail_new/');

npoints = 9;   % npoints-1 segments
frameid1 = 1000;   % frameid used for draw the stationary tail
controlpoint = 1;  % 

% conserved options --------------------------------------------

winsize = 241;  % temporal window for movement detection, better be odd
boxwidth =20;  % spatial zone for movement detection
stdthres = 0;   % 0 for auto estimation


%% options for find tail
isplot = false;
% npoints = 13;   % npoints-1 segments
% params for find_tail

params.prolong_factor = 1.5;
params.turn_scope = 45;
params.turn_step = 3;
params.max_width = 12;   % radius width
params.min_width = 2;
vlim = [0 255];  % better inspect 

% constants
%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% frame1 = imread(fpaths{2});
[frame1, nframes] = readtailseq(imgpath, frameid1);
[sizey, sizex] = size(frame1);
npix = sizey * sizex;
%figure; imshow(frame1)
% figure; imshow(normim(frame1))

%% interactive  %%getting tail axis poisition by drawing
selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);
%figure; imshow(frame1);hold on;plot(selx, sely, 'r-*'); 
save([savepath '/settings_and_constants.mat'], 'npoints', ...
    'winsize', ...
    'boxwidth', ...
    'isplot', ...
    'controlpoint', ...
    'params', ...
    'vlim', ...
    'frame1', ...
    'selx', 'sely');

%% movement detection  %% for indices of all points in the searching zone
disp 'detecting movements'

%hwin    = ceil((winsize-1)/2);
boxwidth =20;
tailind = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth);
a=frame1;a(setdiff(1:size(a,1)*size(a,2),tailind))=0;figure,imshow(a);
%tailind2 = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], boxwidth / 4);
% caution: inline params

tic;
[stdwin, stdwintail] = moving_std_readseq_20171019(imgpath, winsize, tailind);
%figure,plot(stdwin);hold on;plot(stdwintail)
duration_moving_std = toc

% % % moving = false(1, nframes);
% % % moving(movind) = true;
% % % movlabel = bwlabel(moving);
[moving] = estimate_threshold_for_movement_detection(stdwintail, winsize);
%moving=ones(1,244800);
%movlabel = bwlabel(imclose(moving(:)', ones(1, winsize)));
%moving(40:90)=true;moving(193610:193617)=true;
movlabel = bwlabel(moving);
movind = find(movlabel > 0);

save([savepath '/movement_detection.mat'], 'stdwintail', 'stdwin', 'tailind', 'selx', 'sely', 'frame1', 'movind', 'movlabel');
%'tailind2', 'selx', 'sely', 'frame1', 'movind', 'movlabel');

% % % % save out moved frames for inspection
% % % temppath = checkpath([savepath '/moved_frames/']);
% % % for fi = movind(:)'
% % %     if mod(fi, 1000) == 1, disp(fi); end;
% % %     frame = readtailseq(imgpath, fi);
% % %     imwrite(frame, [temppath '/' num2str(fi, '%.6d') '.tif']);
% % % end
%% background estimation
disp 'estimating background'
% pixdistrib moved
tic;
pixdistrib = zeros(sizey*sizex, 256);%count value of per pixels across frames.e.g.pixdistrib(1,1)=4;number of pixel(1,1)=1 is 4
for fi = movind(:)'
    if mod(fi, 1000) == 1, disp(fi); end;
    %frame = imread(fpaths{fi});
    frame = readtailseq(imgpath, fi);
    frame = double(frame);
    frame = frame ./ sum(frame(:));
    frame = normim(frame, [0, 1/100]);
    frame = round(frame * 255);
    
    ind = sub2ind(size(pixdistrib), 1:npix, frame(:)'+1); %%%
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
duration_pixdistrib = toc

% pixdistrib_back = pixdistrib;
% pixdistrib = pixdistrib_back;

figure; bar(sum(pixdistrib))
% vlim = [0, 99];

% estimate background
[estback, quantv] = estimate_background_from_distrib(pixdistrib, vlim, 0.05);
backim = reshape(quantv, size(frame1));
%backim = reshape(estback, size(frame1));

figure; imshow(uint8(round(backim)));
% figure; imshow(normim(backim, [0 0]));

save([savepath '/backim.mat'], 'backim', 'pixdistrib', 'vlim', 'movind');
% save([savepath '/backim2_movind_4301-end.mat'], 'backim', 'pixdistrib', 'movind');

%% background for each case
% todo: imerode movind for a stable movlabel
%%% tind = find(movlabel==38);
%%% im3d = zeros(sizey, sizex, length(movind));
%%% for ii = 1:length(tind)
%%%     fi = tind(ii);
%%%     frame = imread(fpaths{fi});
%%%     im3d(:,:,ii) = frame;
%%% end
%%% 
%%% pixv = reshape(im3d, [sizex*sizey, length(tind)]);
%%% pixd = hist(pixv', 0:255)';
%%% 
%%% figure; bar(sum(pixd))
%%% vlim = [0, 99];
%%% 
%%% estback = estimate_background_from_distrib(pixd, vlim);
%%% 
%%% backim = reshape(estback, size(frame1));

%% find tail: todo: must have guess across time
% % % moving = false(1, nframes);
% % % moving(movind) = true;
% % % moving = imerode(moving, strel('disk', 108));
% % % movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution
%movind=1:nframes;

alltailpos = zeros(npoints, 2, length(movind));

% show
figure;
hdim = imshow(zeros(size(frame1)));
hold on;
hdline = plot(selx, sely, 'r-*');

% initial values
theta0 = atan2d(sely(2)-sely(1), selx(2)-selx(1)) - atan2d(0, 1);  % [-180, 180]
theta0 = - theta0;  % find_tail use right-handed coord
thetap = theta0;
pos0   = [selx(1) sely(1)];
cumlen = cumsum(sqrt(diff(selx).^2 + diff(sely).^2));
taillen = cumlen(end);
seglen = cumlen(end)/(npoints-1);
params.basetheta = theta0;

% main loop
svpath = checkpath([savepath '/find_tail/']);
for ii = 1:length(movind)  
    if mod(ii, 1000) == 1, disp(ii); end
    fi = movind(ii);
    
%     frame = imread(fpaths{fi});
    frame = readtailseq(imgpath, fi);
    frame = double(frame);
    frame = frame ./ sum(frame(:));
    frame = normim(frame, [0, 1/100]);
    frame = round(frame .* 255);  
    
    %engim = abs(frame - backim);
    %[posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen, params);  % fixed init point
    engim = frame;
    [posarray, thetaarray] = find_tail_2im(engim, backim, pos0, thetap, npoints, seglen, params);  % fixed init point
    
    x = posarray(:,1); y = posarray(:,2);

    % smooth tail
    windowWidth = 5;  % caution: inline params
    polynomialOrder = 3;
    smx = sgolayfilt(x, polynomialOrder, windowWidth);
    smy = sgolayfilt(y, polynomialOrder, windowWidth);
%     figure; imshow(engim);
%     hold on; plot(x, y, 'r.')
%     hold on; plot(smx, smy)

%     smx = x; smy = y;
    
    % re-sample points depend on length
    cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
    refx = interp1([0 cumlen'], smx, linspace(0, cumlen(end), npoints), 'pchip');
    refy = interp1([0 cumlen'], smy, linspace(0, cumlen(end), npoints), 'pchip');

    % plot tail movement instantly
    if isplot
        delete(hdline);
        % show results
        %enhframe = mat2gray(frame, vlim);
        %enhframe = adapthisteq(enhframe);
        %set(hdim, 'cdata', enhframe); drawnow;
        set(hdim, 'cdata', normim(engim, [0 0])); drawnow;
        hold on; hdline = plot(refx, refy, 'r-*');
        xlabel(fi);
        pause(0.01);
    end
    
    % save to movie
    %enhframe = mat2gray(frame, vlim);
    %enhframe = adapthisteq(enhframe);
    enhframe = uint8(frame);
    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe = insertShape(enhframe, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    % put info in the movie
%     enhframe = insertText(enhframe, [sizex sizey] - [64, 36], num2str(fi));
%     enhframe = insertText(enhframe, [32, 18], ['case '  num2str(movlabel(fi))]);

    % save movie
%     svpath = checkpath('Z:/tailmovie2/');
    imwrite(enhframe, [svpath '/' num2str(fi, '%.6d') '.tif']);
    
    % results
    posarray = [refx(:) refy(:)];
    thetap = thetaarray(2);
%     tte(ii,:) = posarray(end,:);
%     ttm(ii,:) = posarray(10,:);
    alltailpos(:,:,ii) = posarray;
end

% % % % depricated: more smooth of tails
% % % alltailpos_sm = alltailpos;
% % % for jj = 11:npoints
% % %     windowWidth = 7;
% % %     polynomialOrder = 3;
% % %     smx = sgolayfilt(squeeze(alltailpos(jj,1,:)), polynomialOrder, windowWidth);
% % %     smy = sgolayfilt(squeeze(alltailpos(jj,2,:)), polynomialOrder, windowWidth);
% % %     alltailpos_sm(jj,1,:) = smx;
% % %     alltailpos_sm(jj,2,:) = smy;
% % % end

save([savepath '/alltailpos.mat'], 'alltailpos');

%{
%% make movie with lines
svpath = checkpath([savepath '/tail_with_line/']);
figure;
hdim = imshow(zeros(size([frame1])));
for ii = 1:size(alltailpos,3)
    fi = movind(ii);
    fi
    frame = imread(fpaths{fi});
    frame = double(frame);
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    
    refx = squeeze(alltailpos(:,1,ii));
    refy = squeeze(alltailpos(:,2,ii));

    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe1 = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe1 = insertShape(enhframe1, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    enhframe1 = insertText(enhframe1, [sizex sizey] - [64, 36], num2str(fi));
    enhframe1 = insertText(enhframe1, [32, 18], ['case '  num2str(movlabel(fi))]);
    
    enhframe1 = insertShape(enhframe1, 'line', [pos(controlpoint,1) pos(controlpoint,2) selx(controlpoint) sely(controlpoint)], 'color', 'g');
%     refx = squeeze(alltailpos_sm(:,1,ii));
%     refy = squeeze(alltailpos_sm(:,2,ii));
% 
%     pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
%     enhframe2 = insertShape(enhframe1, 'line', posline, 'color', 'g');
%     enhframe2 = insertShape(enhframe2, 'filledcircle', [pos repmat(2, size(pos,1), 1)], 'color', 'g');
    %enhframe2 = insertText(enhframe2, [sizex sizey] - [64, 36], num2str(fi));
    %enhframe2 = insertText(enhframe2, [32, 18], ['case '  num2str(movlabel(fi))]);
    
    %enhframe = cat(2, enhframe1, enhframe2);
    enhframe = enhframe1;
%     enhframe = enhframe2;

    imwrite(enhframe, [svpath '/' num2str(ii) '.tif']);
    
% % %     % show instantly
% % %     set(hdim, 'cdata', enhframe); drawnow;
end
%}


%% compute magnitude  % check once more
basevec = [selx(controlpoint) sely(controlpoint)] - [selx(1) sely(1)];
%basevec = [1 0] - [0 0];
% b = regress(sely(1:end-1), [ones(npoints-1,1) selx(1:end-1)]);  % basevec from regression
% temp = [0 1] .* b(2) + b(1);
% basevec = [1 temp(2)] - [0 temp(1)];

movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(1); sely(1)], 1, size(alltailpos,3));

% distance of control points
dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));

% theta, anti-clockwise as positive
theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
theta(theta > 180) = theta(theta > 180) - 360;
theta = -theta;
figure; plot(theta);


%% fill moved frames with whole time
% theta00 = atan2d(sely(end) - sely(1), selx(end) - selx(1)) - atan2d(basevec(2), basevec(1));
% theta00 = -theta00;
theta00 = 0;

% % % ddd = dd; % ddd = 1 - theta;
% % % ddd(theta<0) = -ddd(theta<0);
% % % figure; plot(ddd);

alltheta = theta00 * ones(nframes, 1);
alltheta(movind) = theta;
figure; plot(alltheta)

save([savepath '/alltheta.mat'], 'alltheta', 'theta');


%{
%% correct theta (jturn)
moving = false(nframes, 1);
moving(movind) = true;
lastmove = movind(1);
correctind = 1:nframes;
for ii = 1:nframes
    if moving(ii)
       lastmove = ii;
    else
        correctind(ii) = lastmove;
    end
end
correctedtheta = alltheta(correctind);
figure; plot(correctedtheta)
%}

%{
%% moving magnitude
% frmind = 
% degree = 
nframes = length(stdwintail);
nvol = size(dfdf, 2);  % needs calc info

% caution: inline params
tp_calc = 0.6 * 1000;
fs_calc = 1000/tp_calc;
tp_tail = 3.333;
fs_tail = 1000/tp_tail;

fs_tail = 299.68;
tp_tail = 1000/fs_tail;

corfi = interp1(linspace(0, tp_tail*nframes, nframes), 1:nframes, linspace(0, tp_calc*nvol, nvol), 'nearest');

corwin  = floor(fs_tail / fs_calc / 2) * 2 + 1;
hcorwin = (corwin - 1) / 2;
movthres = 0.5;

moving = (stdwintail > 0.5e-7);
% movlabel = bwlabel(moving(:));

movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution
movind = find(moving);

imat = veccolon(corfi - hcorwin, 1, corwin);
imat(imat<1) = 1;
imat(imat>length(moving)) = length(moving);

tailmovement = zeros(1, nframes);
tailmovement(movind) = ddd;
figure; plot(tailmovement);
movmag = mean(abs(tailmovement(imat)), 2);
figure; plot(movmag);

save movemag_rec5.mat movmag tailmovement imat;

%%% 
% % % sig = dfdf(88810,:);
% % % movets = movets(1:nvol)';
% % % figure; plot(sig);
% % % hold on;
% % % plot(movets/50, 'r');
%}

%close all

%% check non mov img 
tic
[~, nframes] = readtailseq(imgpath, 1);
ind=setdiff(1:nframes,movind);
savepath = checkpath([dirpath '/' imgname '_find_non_tail/']);
for ii = ind
    frameid1 = ii;
    [frame1, ~] = readtailseq(imgpath, frameid1);
    imwrite(frame1,fullfile(savepath ,[num2str(ii, '%07d') '.' '.tif']),'Compression','lzw');
    %I=imread(fullfile(savepath ,[num2str(ii, '%07d') '.' '.tif']));
    %figure,imshow(I);
end
toc
close all;
% %% plot
% frame=[];
% frame.per_cycle=20*60;
% frame.us_start=round(13.2*60+1);
% frame.us_dur=[0.001 0.001 0.01 0.01 0.1]*60;
% trial.hab=15;
% trial.acq_block_trial=6;
% trial.acq_block_number=5;
% trial.acq_block_interval=15;
% trial.test=15;
% 
% trial.total=trial.hab+trial.test+trial.acq_block_trial*trial.acq_block_number+trial.acq_block_interval*(trial.acq_block_number-1);
% stimUS=zeros(1,frame.per_cycle*trial.total);ind_US=[];kk=1;
% for ss=1:trial.acq_block_number
%     for tt=1:trial.acq_block_trial
%         ind=(trial.hab+(tt-1)+(ss-1)*(trial.acq_block_interval+trial.acq_block_trial))*frame.per_cycle+frame.us_start;
%         stimUS(ind:ind+max(round(frame.us_dur(ss)),1)-1)=1;%23
%         ind_US(kk)=ind;kk=kk+1;
%     end
% end
% find( stimUS==1);
% [i,~]=readtailseq(imgpath, 1);
% US_I_mos(size(i,1),size(i,2),length(ind_US));
% H=figure;ju 
% n=30;
% for ii=1:length(ind_US)
%     I=zeros(size(i,1),size(i,2),n);
%     for jj=1:n
%     [I(:,:,jj), ~] = readtailseq(imgpath, ind_US(ii)+jj-1);
%     end
%     US_I_m(:,:,ii)=max(I,[],3);
%     subplot(trial.acq_block_trial,trial.acq_block_number,ii),
%     imshow(US_I_m(:,:,ii),[min(min(US_I_m(:,:,ii))) max(max(US_I_m(:,:,ii)))]);
% end
% savepath = checkpath([dirpath '/' imgname '_find_tail/']);
% savefig(H, [savepath '/US RESPONSE','.fig']);
% save([savepath '/US RESPONSE.mat'],'US_I_m');






