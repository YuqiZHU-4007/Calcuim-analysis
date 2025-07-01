
% movement detection -> background image -> tail movement

%% load images
[fpaths, nframes] = ask_for_images;  % caution: image formats/names dependent
% [fpaths, movind] = ask_for_images_moved(nframes);

%% options if not given
savepath = checkpath('./find_tail_temp/');
% savepath = checkpath('z:/data/20151110/tail_results3/');
% movind = 1:100;
% movlabel = ones(1, nframes);

isplot = false;
npoints = 13;   % 16 segments
% params for find_tail
controlpoint = 5;  % 

params.prolong_factor = 1.5;
params.turn_scope = 45;
params.turn_step = 3;
params.max_width = 12;   % radius width
params.min_width = 2;
vlim = [0 99];  % better inspect 

frame1 = imread(fpaths{2});  % generate by movement_detection
[sizey, sizex] = size(frame1);
npix = sizey * sizex;

%% inspect
% % % figure; 
% % % hdim = imshow(frame1);
% % % for fi = movind(:)'
% % %     frame = imread(fpaths{fi});
% % %     set(hdim, 'cdata', frame);
% % %     drawnow;
% % %     xlabel(fi);
% % %     ylabel(movlabel(fi));
% % % end

%% background estimation
% pixdistrib moved
tic;
pixdistrib = zeros(sizey*sizex, 256);
for fi = movind(:)'
    if mod(fi, 1000) == 0, disp(fi); end;
    frame = imread(fpaths{fi});
    ind = sub2ind(size(pixdistrib), 1:npix, double(frame(:))'+1); %%%
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
duration_pixdistrib = toc

% pixdistrib_back = pixdistrib;
% pixdistrib = pixdistrib_back;

figure; bar(sum(pixdistrib))
% vlim = [0, 99];

% estimate background
estback = estimate_background_from_distrib(pixdistrib, vlim, 0.01);
backim = reshape(estback, size(frame1));

figure; imshow(uint8(round(backim)));

save([savepath '/backim.mat'], 'backim', 'pixdistrib', 'vlim');

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

%% interactive  % todo: make this a function
% interactive
%enhframe1 = mat2gray(frame1, vlim);
%enhframe1 = adapthisteq(frame1);
selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);

%% find tail: todo: must have guess across time
% % % moving = false(1, nframes);
% % % moving(movind) = true;
% % % moving = imerode(moving, strel('disk', 108));
% % % movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution

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

% main loop
alltailpos = zeros(npoints, 2, length(movind));
svpath = checkpath([savepath '/find_tail/']);
for ii = 1:length(movind)  %53101:53520
    if mod(ii, 1000) == 0, disp(ii); end
    fi = movind(ii);
    frame = imread(fpaths{fi});
    %frame = double(frame);
    %engim = abs(frame - backim);
    %engim = normim(engim, [1/100, 1/100]);  % caution: inline params
    engim = double(frame);
    
%     figure; imshow(engim)
%     pos0
    %[posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen, params);  % fixed init point
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
        set(hdim, 'cdata', normim(engim)); drawnow;
        hold on; hdline = plot(refx, refy, 'r-*');
        xlabel(fi);
        pause(0.01);
    end
    
    % save to movie
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe = insertShape(enhframe, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    % put info in the movie
    enhframe = insertText(enhframe, [sizex sizey] - [64, 36], num2str(fi));
    enhframe = insertText(enhframe, [32, 18], ['case '  num2str(movlabel(fi))]);

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









