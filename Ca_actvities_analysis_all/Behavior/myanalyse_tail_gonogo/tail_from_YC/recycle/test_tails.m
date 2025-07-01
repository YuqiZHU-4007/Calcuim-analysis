
%% load images
ext = '.bmp';

imdir = uigetdir;

listdir = dir([imdir '/*' ext]);
fnames = {listdir.name};

nframes = length(fnames);

fpaths = fnames;
for fi = 1:nframes
    fpaths{fi} = [imdir '/' num2str(fi-1) '.bmp'];
end

% fpaths = fnames;
% for ii = 1:length(fnames)
%     fpaths{ii} = [imdir '/' fnames{ii}];
% end

frame1 = imread(fpaths{1});
[sizey, sizex] = size(frame1);
npix = sizey * sizex;

%% compute a pixdistrib
pixdistrib = zeros(sizey*sizex, 256);

tic;
for fi = 1:nframes
    if mod(fi, 10000) == 0, disp(fi); toc; end
    frame = imread(fpaths{fi});
%     maxframe = max(cat(3, frame, maxframe), [], 3);
%     minframe = min(cat(3, frame, minframe), [], 3); 
%     avgframe = avgframe + double(frame);
%     stdframe = stdframe + double(frame) .^ 2;
%     frame = medfilt2(frame, [3 3]);
    ind = sub2ind(size(pixdistrib), 1:npix, double(frame(:))'+1);
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
toc;

save('tailim_statistics.mat', 'pixdistrib');

%% 
pixdistrib_back = pixdistrib;
% pixdistrib = pixdistrib_back;

figure; bar(sum(pixdistrib))
vlim = [0, 100];

% normalize
toohigh = sum(pixdistrib_back(:,vlim(2)+1:end), 2);
pixdistrib(:,vlim(2)+1:end) = 0;
pixdistrib(:,vlim(2)) = pixdistrib(:,vlim(2)) + toohigh;

x = 0:255;

e = sum(repmat(x, [npix, 1]) .* pixdistrib/nframes, 2);

figure; imshow(uint8(round(reshape(e, size(frame1)))));

xx = repmat(x, [npix, 1]);
ee = repmat(e, [1, 256]);
pp = pixdistrib/nframes;

sd = sqrt(sum((xx - ee).^2 .*  pp, 2));

figure; imagesc(reshape(sd, size(frame1)));

figure; hist(sd);

noisesd = prctile(sd, 5);

% caution: settings
quantx = 0.0001;  % 
xshift = icdf('norm', quantx, 0, noisesd);

pixcdf = cumsum(pixdistrib, 2) / nframes;
backim = zeros(npix, 1);
tic;
for ii = 1:npix
    incrind = find(diff(pixcdf(ii,:)) > 0);
    %     ls = find(pixcdf(ii,:)==0, 1, 'last');
    %     rs = find(pixcdf(ii,:)==1, 1, 'first');
    
    if length(incrind) < 2,
        backim(ii) = pixcdf(ii, incrind);
    else
        backim(ii) = interp1(pixcdf(ii,incrind), incrind-1, quantx);
    end
end
toc;

backim = reshape(backim - xshift, size(frame1));

figure; imshow(uint8(round(backim)));



%% find tails
% params
npoints = 17;
% seglen  = 0;
theta0  = 45;

% interactive
enhframe1 = mat2gray(frame1, vlim);
enhframe1 = adapthisteq(enhframe1);
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

imshow(enhframe1);
hold on;
plot(smx, smy, 'b-*');

%
cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
refx = interp1([0 cumlen'], smx, linspace(0, cumlen(end), npoints));
refy = interp1([0 cumlen'], smy, linspace(0, cumlen(end), npoints));
seglen = cumlen(end)/npoints;
pos0 = [refx(1) refy(1)];

figure; imshow(enhframe1);
hold on;
plot(refx, refy, 'r-*');


%% detect motion;
framep = frame1;
dframe = zeros(nframes-1, 1);
dframetail = zeros(nframes-1, 1);
sframe = zeros(nframes-1, 1);
sframetail = zeros(nframes-1, 1);
tailind = roi_line_segment(frame1, [refx(1) refy(1)], [refx(end) refy(end)], 50);
mask = zeros(size(frame1));
mask(tailind) = 1;
figure; imshow(enhframe1 .* mask);
tailind2 = roi_line_segment(frame1, [refx(1) refy(1)], [refx(end) refy(end)], 12);

for fi = 2:nframes
    if mod(fi, 10000) == 1000, disp(fi); toc; end
    frame = imread(fpaths{fi});
    frame = medfilt2(frame, [3 3]);
    dframe(fi-1) = mean(abs(frame(:) - framep(:)));
    dframetail(fi-1) = mean(abs(frame(tailind) - framep(tailind)));
    sframe(fi-1) = mean(abs(frame(:) - frame1(:)));
    sframetail(fi-1) = mean(abs(frame(tailind) - frame1(tailind)));
    framep = frame;
end

%% find tail
% show
figure;
hdim = imshow(zeros(size(frame1)));
hold on;
hdline = plot(refx, refy, 'r-*');

% frame0 = imread(fpaths{123950});
% frame0 = medfilt2(frame0, [3 3]);

theta0 = 45;

% main loop 
alltailpos = zeros(npoints, 2, nframes);
% alltheta   = zeros(npoints, nframes);
for fi = 123951:124100  %%12021:12170  %123951:124100  %1:nframes
    if mod(fi, 10000) == 0, disp(fi); toc; end
    frame = imread(fpaths{fi});
    frame = medfilt2(frame, [3 3]);
%     frame(frame > vlim(2)) = vlim(2);
    
    engim = abs(double(frame) - backim);
    engim = normim(engim, [1/100, 1/100]);  % caution: inline params
    
%     figure; imshow(engim)
    
    [refx(1), refy(1)]
    [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, seglen);
    
    x = posarray(:,1); y = posarray(:,2);
    
%     windowWidth = 5;
%     polynomialOrder = 2;
%     smx = sgolayfilt(x, polynomialOrder, windowWidth);
%     smy = sgolayfilt(y, polynomialOrder, windowWidth);
%     figure; imshow(engim);
%     hold on; plot(x, y, 'r.')
%     hold on; plot(smx, smy)

    smx = x; smy = y;
    
    cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
    refx = interp1([0 cumlen'], smx, linspace(0, seglen*npoints, npoints), 'cubic');
    refy = interp1([0 cumlen'], smy, linspace(0, seglen*npoints, npoints), 'cubic');
    
    %figure; 
%     hdim = imshow(engim);
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    set(hdim, 'cdata', enhframe); drawnow;
    delete(hdline);
    hold on; hdline = plot(refx, refy, 'r-*');
    
    posarray = [refx(:) refy(:)];
    
    alltailpos(:,:,fi) = posarray;
    
    theta0 = thetaarray(2);
    
    %pause(0.1)
end

%% show
ind = find(dframetail > 0.5);
ind = find(sframetail > 1 | dframetail > 0.42);
movets = zeros(1, nframes);
movets(ind) = 1;
figure; plot(movets)

ind = 11870:12450

% suspect region
ind = 45320:45850;    % moved
ind = 126200:126900;  % swinged
ind = 191900:192500;  % moved
ind = 194800:195800;  % swinged

figure; 
hdim = imshow(frame1);
for fi = ind(:)'
    frame = imread(fpaths{fi});
    set(hdim, 'cdata', frame);
    drawnow;
    xlabel(fi);
    ylabel(movlabel(fi));
end

%% movement detection
Fs = 300;
Fs = 1000/3.333;

winsize = 301;
hwin    = 150;

tstack = zeros(sizey, sizex, winsize);
for fi = 1:winsize
    frame = imread(fpaths{fi});
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
    if mod(fi, 10000) == 4000, disp(fi); toc; end
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
toc;


moving = (stdwintail > 0.5e-7);
movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));
movind = find(moving);

save('movement_detection.mat', 'stdwintail', 'stdwin', 'stdwintail2', 'tailind');

%% pixdistrib moved
tic;
pixdistrib = zeros(sizey*sizex, 256);
for fi = movind(:)'
    frame = imread(fpaths{fi});
    ind = sub2ind(size(pixdistrib), 1:npix, double(frame(:))'+1);
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
toc;


%% 
pixdistrib_back = pixdistrib;
% pixdistrib = pixdistrib_back;

figure; bar(sum(pixdistrib))
vlim = [0, 99];

% normalize
toohigh = sum(pixdistrib_back(:,vlim(2)+1:end), 2);
pixdistrib(:,vlim(2)+1:end) = 0;
pixdistrib(:,vlim(2)) = pixdistrib(:,vlim(2)) + toohigh;

x = 0:255;

e = sum(repmat(x, [npix, 1]) .* pixdistrib/length(movind), 2);

figure; imshow(uint8(round(reshape(e, size(frame1)))));

xx = repmat(x, [npix, 1]);
ee = repmat(e, [1, 256]);
pp = pixdistrib/length(movind);

sd = sqrt(sum((xx - ee).^2 .*  pp, 2));

figure; hist(sd);

noisesd = prctile(sd, 5);

% caution: settings
quantx = 0.05;  % 
xshift = icdf('norm', quantx, 0, noisesd);

pixcdf = cumsum(pixdistrib, 2) / length(movind);
backim = zeros(npix, 1);
tic;
for ii = 1:npix
    incrind = find(diff(pixcdf(ii,:)) > 0);
    %     ls = find(pixcdf(ii,:)==0, 1, 'last');
    %     rs = find(pixcdf(ii,:)==1, 1, 'first');
    
    if length(incrind) < 2,
        backim(ii) = pixcdf(ii, incrind);
    else
        backim(ii) = interp1(pixcdf(ii,incrind), incrind-1, quantx, 'linear', 0);
    end
end
toc;

backim = reshape(backim - xshift, size(frame1));

figure; imshow(uint8(round(backim)));





%% find tail
for movi = 1:length(movlabel)
    finds = find(movlabel == movi);
    
    gs = fspecial('disk', 31);
    
    tstack = zeros(sizey, sizex, length(finds));
    for ii = 1:length(finds)
        fi = finds(ii);
        frame = imread(fpaths{min([nframes, fi+hwin+1])});
        frame = medfilt2(frame, [3 3]);
        frame = double(frame);
        
%         tstack(:,:,ii) = abs(frame - imfilter(frame, gs));
        tstack(:,:,ii) = frame;
    end
    
    meanim = mean(tstack, 3);
    bs = mean(meanim(tailind2));
    
    tstack = abs(tstack - bs);

    meanim = mean(tstack, 3);
%     meanim = min(tstack, [], 3);
%     meanim = prctile(tstack, 1, 3);
    meanim = imfilter(meanim, gs);
    figure; imshow(normim(meanim))
    
    mvf = zeros(length(finds), 1);
    engstack = zeros(sizey, sizex, length(finds));
    for ii = 1:length(finds)
        fi = finds(ii);
        frame = imread(fpaths{min([nframes, fi+hwin+1])});
%         frame = medfilt2(frame, [3 3]);
        frame = double(frame);
        engstack(:,:,ii) = abs(frame - backim);
    end
    
end

engstack = tstack;

figure; 
hdim = imshow(normim(frame1));
for fi = 1:size(engstack,3)
    frame = engstack(:,:,fi);
    frame = normim(frame, [1/1000, 1/1000]);
    set(hdim, 'cdata', frame);
    drawnow;
    xlabel(fi);
end


%% find tail
% show
figure;
hdim = imshow(zeros(size(frame1)));
hold on;
hdline = plot(refx, refy, 'r-*');

% frame0 = imread(fpaths{123950});
% frame0 = medfilt2(frame0, [3 3]);

theta0 = 45;

% main loop 
alltailpos = zeros(npoints, 2, nframes);
% alltheta   = zeros(npoints, nframes);
for fi = movind(:)'
    frame = imread(fpaths{fi});
    frame = double(frame);
    engim = abs(frame - backim);
    engim = normim(engim, [1/100, 1/100]);  % caution: inline params
         
    figure; imshow(engim)
    
    [refx(1), refy(1)]
    [posarray, thetaarray] = find_tail(engim, pos0, theta0, npoints, seglen);
    
    x = posarray(:,1); y = posarray(:,2);
    
%     windowWidth = 5;
%     polynomialOrder = 2;
%     smx = sgolayfilt(x, polynomialOrder, windowWidth);
%     smy = sgolayfilt(y, polynomialOrder, windowWidth);
%     figure; imshow(engim);
%     hold on; plot(x, y, 'r.')
%     hold on; plot(smx, smy)

    smx = x; smy = y;
    
    cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
    refx = interp1([0 cumlen'], smx, linspace(0, seglen*npoints, npoints), 'cubic');
    refy = interp1([0 cumlen'], smy, linspace(0, seglen*npoints, npoints), 'cubic');
    
    %figure; 
%     hdim = imshow(engim);
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    set(hdim, 'cdata', enhframe); drawnow;
    delete(hdline);
    hold on; hdline = plot(refx, refy, 'r-*');
    
    posarray = [refx(:) refy(:)];
    
    alltailpos(:,:,fi) = posarray;
    
    theta0 = thetaarray(2);
    
    %pause(0.1)
end

% svpath = checkpath('E:/vc/data/fish/');
% for fi = movind(2916:end)'
%     fi
%     copyfile(fpaths{fi}, [svpath '/' num2str(fi) '.bmp']);
% end





