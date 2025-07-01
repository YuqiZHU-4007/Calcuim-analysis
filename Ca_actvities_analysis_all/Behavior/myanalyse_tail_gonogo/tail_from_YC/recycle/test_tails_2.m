

imdir = uigetdir;

listdir = dir([imdir '/*.tiff']);
fnames = {listdir.name};

nframes = length(fnames);

fpaths = fnames;
for ii = 1:length(fnames)
    fpaths{ii} = [imdir '/' fnames{ii}];
end

nkeep = floor(0.1 * nframes);
ncache = 1000;

frame1 = imread(fpaths{1});
[sizey, sizex] = size(frame1);
npix = sizey * sizex;

keepim = zeros(sizey, sizex, nkeep + ncache, 'uint8');

tic;
for fi = 1:nkeep
    frame = imread(fpaths{fi}); 
    keepim(:,:,fi) = frame;
end
toc;

tic;
cachei = 1;
for fi = nkeep+1:nframes
    fi
    frame = imread(fpaths{fi});
    if cachei == ncache
        keepim = sort(keepim, 3);
        cachei = 1;
    else
        keepim(:,:,nkeep+cachei) = frame;
        cachei = cachei + 1;
    end
end
toc;



% compute a min and a max
maxframe = 0 * ones(sizey, sizex, 'uint8');
minframe = 255 * ones(sizey, sizex, 'uint8');
avgframe = zeros(sizey, sizex);
stdframe = zeros(sizey, sizex);

% matlabpool(6);

pixdistrib = zeros(sizey*sizex, 256);

tic;
for fi = 1:nframes
    if mod(fi, 10000) == 0, disp(fi); toc; end
    frame = imread(fpaths{fi});
%     maxframe = max(cat(3, frame, maxframe), [], 3);
%     minframe = min(cat(3, frame, minframe), [], 3); 
%     avgframe = avgframe + double(frame);
%     stdframe = stdframe + double(frame) .^ 2;
    frame = medfilt2(frame, [3 3]);
    ind = sub2ind(size(pixdistrib), 1:npix, frame(:)'+1);
    pixdistrib(ind) = pixdistrib(ind) + 1;
end
toc;

avgframe = avgframe / nframes;
stdframe = stdframe / nframes - (avgframe.^2);




figure; imshow(frame1);
h = impoly();
pos = wait(h);

% interpolate for t



imwrite(maxframe, 'maxframe.tif');
imwrite(minframe, 'minframe.tif');
imwrite(avgframe, 'avgframe.tif');
imwrite(stdframe, 'stdframe.tif');

save('tailim_statistics.mat', 'avgframe', 'stdframe');


%% inspect
pixdistrib_back = pixdistrib;


figure; bar(sum(pixdistrib))

% normalize
pixdistrib = pixdistrib_back;
toohigh = sum(pixdistrib_back(:,91:end), 2);
pixdistrib(:,91:end) = 0;
pixdistrib(:,90) = pixdistrib(:,90) + toohigh;


i2 = sub2ind(size(cdata), 283, 312);
i3 = sub2ind(size(frame1), 78, 356);
figure; bar(pixdistrib(i3,:));

ii = i3

x = 0:255;
% e = sum(x .* pixdistrib(ii, :)/nframes);
% sd = sum((x - e) .^ 2 .*  pixdistrib(ii, :)/nframes);

e = sum(repmat(x, [npix, 1]) .* pixdistrib/nframes, 2);

xx = repmat(x, [npix, 1]);
ee = repmat(e, [1, 256]);
pp = pixdistrib/nframes;

sd = sqrt(sum((xx - ee).^2 .*  pp, 2));


sdsd = repmat(sd, [1, 256]);
sk =  sum(((xx - ee)./sdsd).^3 .* pp, 2);

figure; hist(sk(:))

skmap = reshape(sk, [sizey sizex]);
skmap(isnan(skmap)) = 0;
figure; hist(skmap(:))
skmap = normim(skmap, [1/200, 1/200]);

figure; imagesc(skmap)



[~, ski] = sort(skmap(:));
figure; bar(pixdistrib(ski(end),:))











