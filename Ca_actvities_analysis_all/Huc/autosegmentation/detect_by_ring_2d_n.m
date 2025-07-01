function [points, scoremap, radmap] = detect_by_ring(slice, radii, thres, mask)
% [points, scoremap, radmap] = detect_by_ring(slice, radii, thres)
%
% points:   logical image
% scoremap: score image
% radmap  : radius map
%

if nargin < 4
    mask = true(size(slice));
end

if nargin < 3
    thres = [];  %4/256;
end

ntry = length(radii);
fimArray = zeros(ntry, numel(slice));
rgimArray = zeros(ntry, numel(slice));
%varimArray = zeros(ntry, numel(slice));
for ri = 1:length(radii)
    rad = radii(ri);
    ft = double(strel('disk', rad).getnhood);
    fim = imfilter(slice, ft);
% % %     bim = rangefilt(slice, ft);
    fimArray(ri,:) = fim(:);
% % %     rgimArray(ri,:) = bim(:);
    %varim = stdfilt(slice, ft);
    %varimArray(ri,:) = varim(:);
end

des = zeros(ntry^2/2-ntry/2, numel(slice));
ki = 1;
for ii = 1:ntry-1
    for jj = ii+1:ntry
        nin = sum(sum(strel('disk', radii(ii)).getnhood));
        nex = sum(sum(strel('disk', radii(jj)).getnhood));
        %difring = (fimArray(jj,:) - fimArray(ii,:)) ./ (nex-nin) - fimArray(ii,:) ./ nin;
        difring =  (fimArray(ii,:) ./ nin).^2 - ((fimArray(jj,:) - fimArray(ii,:)) ./ (nex-nin)).^2;
% % %         varring = rgimArray(jj,:);
        %varring = (varimArray(jj,:) .* (nex - 1) - varimArray(ii,:) .* (nin - 1)) / (nex - nin - 1);
        des(ki,:) = difring;  % ./ varring;
        %des(ki,:) = difring ./ (fimArray(jj)./nin);  % dfdf style
        %des(ki,:) = (fimArray(jj,:) - fimArray(ii,:)) - fimArray(ii,:);
        ki = ki + 1;
    end
end

if size(des, 1) == 1
    fim = des;
    radin = radii(1) * ones(size(slice));
    radex = radii(2) * ones(size(slice));
else
    [fim, who] = max(des);
    [ii, jj] = UTind2sub(who, size(des,1));
    radin = radii(ii);
    radex = radii(jj);
    radin = reshape(radin, size(slice));
    radex = reshape(radex, size(slice));
end

fim = reshape(fim, size(slice));
%radmap = reshape(radmap, size(slice));
radmap = cat(3, radin, radex);

ft = strel('disk', radii(1)+1).getnhood;
ffim = ordfilt2(fim, sum(ft(:)), ft);
points = (ffim == fim);
%points = imextendedmax(fim, 0, 8);
points(~mask) = false;

if isempty(thres)
    disp aaa;
    temp = ffim(points);
    temp = temp(temp>0);
    thres = median(temp);
end
    
points(fim < thres) = false;
scoremap = fim;
%radmap = radmap;

end

function [i, j] = UTind2sub(k, N)
% [i, j] = UTind2sub(k, N)
% N is the length of the 1-d representation of the matrix
%

    n = (sqrt(8*N+1) + 1) / 2;

    b = - 2*n - 1;
    ac = 2*n + 2*k;

    i = ceil((-b - sqrt(b^2 - 4*ac)) / 2 - 1);

    j = k - (2*n-i).*(i-1)/2 + i;

end



