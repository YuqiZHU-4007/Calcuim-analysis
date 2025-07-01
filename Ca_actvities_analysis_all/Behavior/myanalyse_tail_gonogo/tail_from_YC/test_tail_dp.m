
% copy out moved images

%TODO: save params

%% load images
[imgname, dirpath] = uigetfile('*.seq');
imgpath = [dirpath imgname];

% [~, name, ext] = fileparts(imgname);
savepath = checkpath([dirpath '/' imgname '_find_tail_dp/']);

%% options for movement detections
% dirpath = fileparts(fpaths{1});
% [~, dirname] = fileparts(dirpath);
% savepath = checkpath([dirpath '/../' imgname '_find_tail/']); %% /.. stands for parental folder
% savepath = checkpath('Z:/data/20150625/rec5338806_sig/tail_new/');

npoints = 11;   % npoints-1 segments
frameid1 = 1;   % frameid used for draw the stationary tail
controlpoint = 6;  % 

% conserved options --------------------------------------------

winsize = 241;  % temporal window for movement detection, better be odd
boxwidth = 50;  % spatial zone for movement detection
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
% frame1 = imread(fpaths{2});
[frame1, nframes] = readtailseq(imgpath, frameid1);
[sizey, sizex] = size(frame1);
npix = sizey * sizex;

% figure; imshow(normim(frame1))

%% interactive  %%getting tail axis poisition by drawing
selxy = ask_for_tail(frame1, npoints);
selx = selxy(:,1);  sely = selxy(:,2);

save([savepath '/settings_and_constants.mat'], 'npoints', ...
                                               'winsize', ...
                                               'boxwidth', ...
                                               'isplot', ...
                                               'controlpoint', ...
                                               'params', ...
                                               'vlim', ...
                                               'frame1', ...
                                               'selx', 'sely');

                                           
                                           
                                           
                                           
%% estimate backim from one single image
winsize = 15;
[ind] = roi_line_segment(frame1, [selx(1) sely(1)], [selx(end) sely(end)], winsize);
[y, x] = ind2sub(size(frame1), ind);

mask = zeros(size(frame1));
mask(ind) = 1;
figure; imshow(mask);
figure; imshow(mask .* im2double(frame1));

% theta = 30;
% rotmat1 = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
% rotmat2 = [cosd(-theta) sind(-theta); -sind(-theta) cosd(-theta)];
backim = double(frame1);
for ii = 1:length(ind)
    pos0 = [x(ii) y(ii)];
%     pos0 = [x(ii) - selx(1) y(ii) - sely(1)];
%     rp1 = round(pos0 * rotmat1 + [selx(1) sely(1)]);
%     rp2 = round(pos0 * rotmat2 + [selx(1) sely(1)]);
%     

    direc = [- (sely(2) - sely(1)), selx(2) - selx(1)];
    direc = direc / norm(direc) * (winsize * 2 + 1);
    
    rp1 = round(pos0 + direc);
    rp2 = round(pos0 - direc);
    
    outbox1 = any(rp1) < 0 || rp1(1) > size(frame1,2) || rp1(2) > size(frame1,1);
    outbox2 = any(rp1) < 0 || rp1(1) > size(frame1,2) || rp1(2) > size(frame1,1);
    if outbox1 & outbox1
        error('');
    elseif outbox1
        rp = rp2;
    elseif outbox2
        rp = rp1;
    else
        rp = [rp1; rp2];
    end
    rind = sub2ind(size(frame1), rp(:,2), rp(:,1));
    backim(ind(ii)) = mean(frame1(rind));
end
    
hold on; plot(rp(:,1), rp(:,2), 'rx');

backim = normim(backim, [0, 1/100]);
figure; imshow((backim))



%%
vec = [selx(controlpoint) - selx(1) sely(controlpoint) - sely(1)];
ttheta = atan2d(vec(2), vec(1)) - atan2d(1, 1);  % [-180, 180]

resize = 1/8;  % caution

xyres = imresize_xy(engim, resize, [selx sely]);
xyrot = imrotate_xy(engimr, ttheta, xyres);

cumham = cumsum(abs(diff(xyrot(:,1))) + abs(diff(xyrot(:,2))));
hamminglen = round(cumham(end));

tail = [selx, sely];

svpath = checkpath([savepath '/find_tail/']);

for fi = movind(:)'
    tic;
    
    frame = readtailseq(imgpath, fi);
    frame = normim(frame, [0, 1/100]);
%     figure; imshow(frame);

    engim = abs(double(frame) - backim);
    % engim = imfilter(engim, ones(11));
    % selx = selxy(:,1);  sely = selxy(:,2);
%     figure; imshow(normim(engim));
%     hold on; plot(selx, sely, 'r.')

    engimr = imresize(engim, resize, 'box');
    
    xy0 = imrotate_xy(engimr, ttheta, xyres(1,:));
    
    engimr = imrotate(engimr, ttheta, 'crop');

%     engimr = imfilter(engimr, fspecial('average', 3));
%     figure; imshow(normim(engimr));

    
    tail = find_tail_dp(engimr, xy0, hamminglen);
    
%     figure; imshow(normim(engimr));
%     hold on; plot(tail(:,1), tail(:,2), 'r')

% % %     tailxr = tail(:,1); tailyr = tail(:,2);
% % %     tailx = (tailxr - size(engimr,2)/2) * resize + size(frame1,2)/2;
% % %     taily = (tailyr - size(engimr,1)/2) * resize + size(frame1,1)/2;
    
%     figure; imshow(normim(engim));
%     hold on; plot(tailx, taily, 'r.')

    enhframe = engimr;
%     enhframe = im2uint8(frame);
%     enhframe = imresize(enhframe, resize, 'box');
%     enhframe = imrotate(enhframe, ttheta, 'crop');
%     ind = sub2ind(size(enhframe), tail(:,2), tail(:,1));
%     enhframe(ind) = 255;
    posline = tail'; posline = posline(:)';
    enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
    imwrite(enhframe, [svpath '/' num2str(fi, '%.6d') '.tif']);
    
% % %     pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
% % %     enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
% % %     enhframe = insertShape(enhframe, 'circle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
% % %     % put info in the movie
% % %     enhframe = insertText(enhframe, [sizex sizey] - [64, 36], num2str(fi));
% % %     enhframe = insertText(enhframe, [32, 18], ['case '  num2str(movlabel(fi))]);

    ccp = ceil(hamminglen/2);
    vec = tail(ccp,:) - tail(1,:);
    ttheta = (atan2d(vec(2), vec(1)) - atan2d(1, 1)) + ttheta;
    
    toc;

end



