
% movement detection -> background image -> tail movement 


%% load images
ext = '.bmp';

imdir = uigetdir;

listdir = dir([imdir '/*' ext]);

% fnames = {listdir.name};
% fpaths = fnames;
% for ii = 1:length(fnames)
%     fpaths{ii} = [imdir '/' fnames{ii}];
% end

% movind = cellfun(@(x) str2num(x(1:end-4)), fnames) + 1;

nframes = length(listdir);

fpaths = cell(nframes, 1);
for fi = movind(:)'    % caution: fmt dependent
    fpaths{fi} = [imdir '/' num2str(fi-1) '.bmp'];
end

frame1 = imread(fpaths{movind(1)});
[sizey, sizex] = size(frame1);
npix = sizey * sizex;


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
for fi = movind(1:end)'
    fi
    frame = imread(fpaths{fi});
    ind = sub2ind(size(pixdistrib), 1:npix, double(frame(:))'+1);
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
toc;

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

save('backim.mat', 'backim');


%% find tail
% show
figure;
hdim = imshow(zeros(size(frame1)));
hold on;
hdline = plot(selx, sely, 'r-*');

% main loop 
thetap = theta0;
alltailpos = zeros(npoints, 2, nframes);
% alltheta   = zeros(npoints, nframes);
for fi = movind(150:end)  %53101:53520 %movind(2000:end)'
    frame = imread(fpaths{fi});
    frame = double(frame);
    engim = abs(frame - backim);
    engim = normim(engim, [1/100, 1/100]);  % caution: inline params
         
%     figure; imshow(engim)
    pos0
    [posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen);  % fixed init point
    
    x = posarray(:,1); y = posarray(:,2);
        
    windowWidth = 5;
    polynomialOrder = 3;
    smx = sgolayfilt(x, polynomialOrder, windowWidth);
    smy = sgolayfilt(y, polynomialOrder, windowWidth);
%     figure; imshow(engim);
%     hold on; plot(x, y, 'r.')
%     hold on; plot(smx, smy)

%     smx = x; smy = y;
    
    cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
    refx = interp1([0 cumlen'], smx, linspace(0, taillen, npoints), 'cubic');
    refy = interp1([0 cumlen'], smy, linspace(0, taillen, npoints), 'cubic');
    
    %figure; 
%     hdim = imshow(engim);
    delete(hdline);
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    set(hdim, 'cdata', enhframe); drawnow;
    hold on; hdline = plot(refx, refy, 'r-*');
    xlabel(fi);
    pause(0.01)
    
    posarray = [refx(:) refy(:)];
    
    alltailpos(:,:,fi) = posarray;
    
    thetap = thetaarray(2);
    
    %pause(0.1)
end

% svpath = checkpath('E:/vc/data/fish/');
% for fi = movind(2916:end)'
%     fi
%     copyfile(fpaths{fi}, [svpath '/' num2str(fi) '.bmp']);
% end



%% must have guess across time
moving = false(1, nframes);
moving(movind) = true;
moving = imerode(moving, strel('disk', 108));
movlabel = bwlabel(imclose(moving(:)', ones(1, 301)));  % caution

% show
figure;
hdim = imshow(zeros(size(frame1)));
% hold on;
% hdline = plot(selx, sely, 'r-*');
thetap = theta0;

tind = find(movlabel > 0);
% tte = zeros(length(tind), 2);
% ttm = zeros(length(tind), 2);
alltailpos = zeros(npoints, 2, length(tind));
for ii = 1:length(tind)  %53101:53520 %movind(2000:end)'
    fi = tind(ii);
    fi
    frame = imread(fpaths{fi});
    frame = double(frame);
    engim = abs(frame - backim);
    engim = normim(engim, [1/100, 1/100]);  % caution: inline params
         
%     figure; imshow(engim)
%     pos0
    [posarray, thetaarray] = find_tail(engim, pos0, thetap, npoints, seglen);  % fixed init point
    
    x = posarray(:,1); y = posarray(:,2);

    windowWidth = 5;
    polynomialOrder = 3;
    smx = sgolayfilt(x, polynomialOrder, windowWidth);
    smy = sgolayfilt(y, polynomialOrder, windowWidth);
%     figure; imshow(engim);
%     hold on; plot(x, y, 'r.')
%     hold on; plot(smx, smy)

%     smx = x; smy = y;
    
    cumlen = cumsum(sqrt(diff(smx).^2 + diff(smy).^2));
    refx = interp1([0 cumlen'], smx, linspace(0, taillen, npoints), 'pchip');
    refy = interp1([0 cumlen'], smy, linspace(0, taillen, npoints), 'pchip');
    
    %figure; 
%     hdim = imshow(engim);
%     delete(hdline);
%     hold on; hdline = plot(refx, refy, 'r-*');

    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe = insertShape(enhframe, 'filledcircle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
    enhframe = insertText(enhframe, [sizex sizey] - [64, 36], num2str(fi));
    enhframe = insertText(enhframe, [32, 18], ['case '  num2str(movlabel(fi))]);
%     set(hdim, 'cdata', enhframe); drawnow;
%     xlabel(fi);
%     svpath = checkpath('Z:/tailmovie2/');
%     imwrite(enhframe, [svpath '/rec5_mov_' num2str(fi, '%.6d') '.tif']);
    
    posarray = [refx(:) refy(:)];
    thetap = thetaarray(2);
    
%     tte(ii,:) = posarray(end,:);
%     ttm(ii,:) = posarray(10,:);
    alltailpos(:,:,ii) = posarray;
end

alltailpos_sm = alltailpos;
for jj = 11:npoints
    windowWidth = 7;
    polynomialOrder = 3;
    smx = sgolayfilt(squeeze(alltailpos(jj,1,:)), polynomialOrder, windowWidth);
    smy = sgolayfilt(squeeze(alltailpos(jj,2,:)), polynomialOrder, windowWidth);
    alltailpos_sm(jj,1,:) = smx;
    alltailpos_sm(jj,2,:) = smy;
end

controlpoint = 10;

figure;
hdim = imshow(zeros(size([frame1])));
for ii = 1:size(alltailpos,3)
    fi = tind(ii);
    fi
    frame = imread(fpaths{fi});
    frame = double(frame);
    enhframe = mat2gray(frame, vlim);
    enhframe = adapthisteq(enhframe);
    
    refx = squeeze(alltailpos(:,1,ii));
    refy = squeeze(alltailpos(:,2,ii));

    pos = [refx(:) refy(:)]; posline = pos'; posline = posline(:)';
    enhframe1 = insertShape(enhframe, 'line', posline, 'color', 'r');
    enhframe1 = insertShape(enhframe1, 'filledcircle', [pos repmat(2, size(pos,1), 1)], 'color', 'r');
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
    enhframe = enhframe1
    enhframe = enhframe2;
    
    svpath = checkpath('Z:/testtail2/');
    imwrite(enhframe, [svpath '/' num2str(ii) '.tif']);
    
    set(hdim, 'cdata', enhframe); drawnow;
end

basevec = [selx(controlpoint) sely(controlpoint)] - [selx(1) sely(1)];
movevec = squeeze(alltailpos(controlpoint,:,:)) - repmat([selx(1); sely(1)], 1, length(tind));
theta = basevec * movevec ./ norm(basevec) ./ normmove;
normmove = sqrt(sum(movevec.^2));
./ sqrt(movevec' * movevec);  %dot(basevec, movevec);

difvec = squeeze(alltailpos(controlpoint,:,:))' - repmat([selx(controlpoint) sely(controlpoint)], length(tind), 1);

dd = dist([selx(controlpoint) sely(controlpoint)], squeeze(alltailpos(controlpoint,:,:)));
ddd = dd;
ddd(difvec(:,2) >= 0 & difvec(:,1) >= 0) = -ddd(difvec(:,2) >= 0 & difvec(:,1) >= 0);
figure; plot(ddd);
ddd = ddd';


%% moving magnitude
% frmind = 
% degree = 
movmag = zeros(1, nframes);
movmag(frmind) = degree;
m = mean(movmag(imat), 2);


%%% 
sig = dfdf(88810,:);
movets = movets(1:nvol)';
figure; plot(sig);
hold on;
plot(movets/50, 'r');










