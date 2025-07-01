% 20151206 from YangChen
% movement detection -> background image -> tail movement 

%% load images
% % % ext = '.tif';
% % % 
% % % nframes = 0;
% % % outdir = 'H:\YinChen\Behavior_raw_data\20151201_Beha_data\20151201_Beha_fish_1';
% % % for ii = 1:10
% % %     dirpath = [outdir '/20151201_Beha_fish_1_trial_' num2str(ii)];
% % %     nframes = nframes + length(dir(dirpath)) - 2;
% % % end
% % % fpaths = cell(nframes, 1);
% % % k = 1;
% % % for ii = 1:10
% % %     dirpath = [outdir '/20151201_Beha_fish_1_trial_' num2str(ii)];
% % %     for jj = 1:3000
% % %         fpaths{k} = [dirpath '/RecordedImage_GO-5000M-USB__' num2str(jj-1, '%.3d') ext];
% % %         k = k + 1;
% % %     end
% % % end

%% options if not given
% savepath = checkpath('Z:/light-sheet/jturn/test_detect/');
% movind = 1:nframes;
% movlabel = ones(1, nframes);

isplot = false;
npoints = 9;   % 16 segments
% params for find_tail
controlpoint = 9;  % 

params.prolong_factor = 1.5;
params.turn_scope = 40;
params.turn_step = 2;
params.max_width = 11;
params.min_width = 3;

    
% % % fpaths = cell(nframes, 1);
% % % %for fi = movind(:)'    % caution: fmt dependent
% % % %    fpaths{fi} = [imdir '/' num2str(fi-1) '.bmp'];
% % % %end
% % % % for fi = movind(:)'    % caution: fmt dependent
% % % %     fpaths{fi} = [imdir '/' listdir(fi).name];
% % % % end
% % % for fi = movind(:)'    % caution: fmt dependent
% % %     fpaths{fi} = [imdir '/RecordedImage_GO-5000M-USB__' num2str(fi-1, '%.3d') ext];
% % % end

frame1 = 255 - imread(fpaths{movind(1)});  % generate by movement_detection
[sizey, sizex] = size(frame1);
npix = sizey * sizex;
% figure; imshow(frame1);

%% inspect moved
% % % figure; 
% % % hdim = imshow(frame1);
% % % for fi = movind(:)'
% % %     frame = imread(fpaths{fi});
% % %     set(hdim, 'cdata', frame);
% % %     drawnow;
% % %     xlabel(fi);
% % %     ylabel(movlabel(fi));
% % % end

%% pixdistrib moved
tic;
pixdistrib = zeros(sizey*sizex, 256);
for fi = movind(:)'
    if mod(fi, 1000) == 0, disp(fi); end;
    frame = 255 - imread(fpaths{fi});
    ind = sub2ind(size(pixdistrib), 1:npix, double(frame(:))'+1);
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
duration_pixdistrib = toc

% pixdistrib_back = pixdistrib;
% pixdistrib = pixdistrib_back;

figure; bar(sum(pixdistrib))
vlim = [0, 255];   % this should depend on the pixdistrib

estback = estimate_background_from_distrib(pixdistrib, vlim);

backim = reshape(estback, size(frame1));

figure; imshow(uint8(round(backim)));

save([savepath '/backim.mat'], 'backim', 'pixdistrib', 'vlim');

%% background for each case
% todo: imerode movind for a stable movlabel

% % % tind = find(movlabel==38);
% % % im3d = zeros(sizey, sizex, length(tind));
% % % for ii = 1:length(tind)
% % %     fi = tind(ii);
% % %     frame = imread(fpaths{fi});
% % %     im3d(:,:,ii) = frame;
% % % end
% % % 
% % % pixv = reshape(im3d, [sizex*sizey, length(tind)]);
% % % pixd = hist(pixv', 0:255)';
% % % 
% % % figure; bar(sum(pixd))
% % % vlim = [0, 99];
% % % 
% % % estback = estimate_background_from_distrib(pixd, vlim);
% % % 
% % % backim = reshape(estback, size(frame1));


%% interactive  % todo: make this a function
% enhframe1 = mat2gray(frame1, vlim);
enhframe1 = adapthisteq(frame1);
figure; imshow(enhframe1);
h = impoly();
pos = wait(h);
x = pos(:,1); y = pos(:,2);

if size(pos, 1) >= 3
    windowWidth = 3;
    polynomialOrder = 2;
    smx = sgolayfilt(x, polynomialOrder, windowWidth);
    smy = sgolayfilt(y, polynomialOrder, windowWidth);
else
    smx = x; smy = y;
end
% smx = x; smy = y;
imshow(enhframe1);
hold on;
plot(smx, smy, 'b-*');

% smooth tail
cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
taillen = cumlen(end);
selx = interp1([0 cumlen'], smx, linspace(0, cumlen(end), npoints));
sely = interp1([0 cumlen'], smy, linspace(0, cumlen(end), npoints));
seglen = cumlen(end)/(npoints-1);
pos0 = [selx(1) sely(1)];

figure; imshow(enhframe1);
hold on;
plot(selx, sely, 'r-*');

%% find tail: todo: must have guess across time
% % % % too many frames withou movement
% % % moving = false(1, nframes);
% % % moving(movind) = true;
% % % moving = imerode(moving, strel('disk', 108));
% % % movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution

% show
figure;
hdim = imshow(zeros(size(frame1)));
hold on;
hdline = plot(selx, sely, 'r-*');

theta0 = atan2d(sely(2)-sely(1), selx(2)-selx(1)) - atan2d(0, 1);  % [-180, 180]
theta0 = - theta0;  % find_tail use right-handed coord
thetap = theta0;
% % % pos0   = [selx(1) sely(1)];
% % % cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
% % % seglen = cumlen(end)/(npoints-1);

% % % tind = find(movlabel > 0);
tind = movind;
% tte = zeros(length(tind), 2);
% ttm = zeros(length(tind), 2);

alltailpos = zeros(npoints, 2, length(tind));
svpath = checkpath([savepath '/find_tail/']);
for ii = 1:length(tind)  %53101:53520 %movind(2000:end)'
    fi = tind(ii);
    fi
    frame = 255 - imread(fpaths{fi});
    frame = double(frame);
    engim = abs(frame - backim);
    engim = normim(engim, [1/100, 1/100]);  % caution: inline params
         
%     figure; imshow(engim)
%     pos0
    [posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen, params);  % fixed init point
    
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
    refx = interp1([0 cumlen'], smx, linspace(0, taillen, npoints), 'pchip');
    refy = interp1([0 cumlen'], smy, linspace(0, taillen, npoints), 'pchip');
    
    % plot tail movement instantly
    if isplot
        delete(hdline);
        % show results
        %enhframe = mat2gray(frame, vlim);
        %enhframe = adapthisteq(enhframe);
        %set(hdim, 'cdata', enhframe); drawnow;
        set(hdim, 'cdata', engim); drawnow;
        hold on; hdline = plot(refx, refy, 'r-*');
        xlabel(fi);
        pause(0.01);
    end
    
    % save to movie
    enhframe = 1 - mat2gray(frame, vlim);
    %enhframe = adapthisteq(enhframe);
    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe = insertShape(enhframe, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    % put info in the movie
    enhframe = insertText(enhframe, [sizex sizey] - [64, 36], num2str(fi));
    %enhframe = insertText(enhframe, [32, 18], ['case '  num2str(movlabel(fi))]);

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

%% make movie with lines
svpath = checkpath([savepath '/tail_with_line/']);
% figure;
% hdim = imshow(zeros(size(frame1)));
for ii = 1:size(alltailpos,3)
    fi = tind(ii);
    fi
    frame = 255 - imread(fpaths{fi});
    frame = double(frame);
    enhframe = 1 - mat2gray(frame, vlim);
    %enhframe = adapthisteq(enhframe);
    
    refx = squeeze(alltailpos(:,1,ii));
    refy = squeeze(alltailpos(:,2,ii));

    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe1 = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe1 = insertShape(enhframe1, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    enhframe1 = insertText(enhframe1, [sizex sizey] - [64, 36], num2str(fi));
    %enhframe1 = insertText(enhframe1, [32, 18], ['case '  num2str(movlabel(fi))]);
    
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


%% compute magnitude  % check once more
% basevec = [selx(controlpoint) sely(controlpoint)] - [selx(1) sely(1)];
%basevec = [1 0] - [0 0];
b = regress(sely(1:end-1)', [ones(npoints-1,1) selx(1:end-1)']);
temp = [0 1] .* b(2) + b(1);
basevec = [1 temp(2)] - [0 temp(1)];
movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(1); sely(1)], 1, length(tind));

% distance of control points
dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));

% theta, anti-clockwise as positive
theta = atan2d(movevec(2,:), movevec(1,:)) - atan2d(basevec(2), basevec(1));
theta = -theta;
figure; plot(theta);

theta00 = atan2d(sely(end) - sely(1), selx(end) - selx(1)) - atan2d(basevec(2), basevec(1));
theta00 = -theta00;

% % % ddd = dd; % ddd = 1 - theta;
% % % ddd(theta<0) = -ddd(theta<0);
% % % figure; plot(ddd);

alltheta = theta00 * ones(nframes, 1);
alltheta(movind) = theta;
figure; plot(alltheta)

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


%%
% % % thres = 64;
% % % %backint = [];
% % % conn = bwlabel(backint > thres);
% % % connum = hist(conn, unique(conn));
% % % [maxv, maxi] = max(connum(2:end));
% % % infered_region = conn == 10;
% % % figure; plot(infered_region);




