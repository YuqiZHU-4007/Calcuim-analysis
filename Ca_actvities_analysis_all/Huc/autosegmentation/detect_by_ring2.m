function [points, fim, radmap] = detect_by_ring2(slice, radii, thres, spacing)

if nargin < 4
    spacing = radii(1)*2-1;
end
if nargin < 3
    thres = 5;
end

ntry = length(radii);
fimArray = zeros(ntry, numel(slice));
for ri = 1:length(radii)
    rad = radii(ri);
    ft = double(strel('disk', rad).getnhood);
    fim = imfilter(slice, ft);
    fimArray(ri,:) = fim(:);
end

des = zeros(ntry^2/2-ntry/2, numel(slice));
ki = 1;
for ii = 1:ntry-1
    for jj = ii+1:ntry
        ni = sum(sum(strel('disk', radii(ii)).getnhood));
        no = sum(sum(strel('disk', radii(jj)).getnhood));
        des(ki,:) = (fimArray(jj,:) - fimArray(ii,:)) ./ (no-ni) - fimArray(ii,:) ./ ni;
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
%radmap = radii(who);

fim = reshape(fim, size(slice));
%radmap = reshape(radmap, size(slice));
radmap = cat(3, radin, radex);

ft = strel('disk', spacing).getnhood;
ffim = ordfilt2(fim, sum(ft(:)), ft);
points = (ffim == fim);
%points = imextendedmax(fim, 0, 8);

points(fim <= thres) = false;


