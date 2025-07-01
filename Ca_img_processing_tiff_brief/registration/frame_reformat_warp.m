function [targim, xshift, yshift] = frame_reformat_warp(targim, params, warpingSettings)
% [correctedImage] = frame_reformat_warp(targim, params, warpingSettings)
%

[height, width] = size(targim);

%% warp
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
%maximumShift = warpingSettings(3);

[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));
dxArray = params(:,1);
dyArray = params(:,2);

%% reformat
xshift = zeros(size(targim));
yshift = zeros(size(targim));

movedtileind = find(any(params~=0, 2));
for q = movedtileind(:)'  %%determine the maximal shift in the adjacent area of each pixel, as Ahrens'
%for q = 1:numel(xx)
    jj = xx(q); ii = yy(q);
    iRange = max(1, ii - squareSize):min(height, ii + squareSize);
    jRange = max(1, jj - squareSize):min(width, jj + squareSize);
    
    xMaxshiftSquare = xshift(iRange, jRange);
    yMaxshiftSquare = yshift(iRange, jRange);
    
    xMaxshiftSquare(abs(xMaxshiftSquare) < abs(dxArray(q))) = dxArray(q);
    xshift(iRange, jRange) = xMaxshiftSquare;
    
    yMaxshiftSquare(abs(yMaxshiftSquare) < abs(dyArray(q))) = dyArray(q);
    yshift(iRange, jRange) = yMaxshiftSquare;
end

movedpixind = find(xshift ~=0 | yshift ~=0);
%[y0, x0] = find(xshift ~=0 | yshift ~=0);
[y0, x0] = ind2sub(size(targim), movedpixind);
x1 = x0 - xshift(movedpixind);
y1 = y0 - yshift(movedpixind);

validind = find(~(x1<1 | x1>width | y1<1 | y1>height));
x0 = x0(validind); y0 = y0(validind);
x1 = x1(validind); y1 = y1(validind);

i0 = sub2ind(size(targim), y0, x0);
i1 = sub2ind(size(targim), y1, x1);

%correctedImage = targim;
targim(i0) = targim(i1);

%{
[x0, y0] = meshgrid(1:width, 1:height);
x1 = x0 - xshift;
y1 = y0 - yshift;
i0 = sub2ind(size(targim), y0, x0);
i1 = sub2ind(size(targim), y1, x1);
correctedImage = targim;
correctedImage(i0) = targim(i1);

shiftvec = (1:numel(targim))';
shiftvec(i0) = i1;
%}

%{
noShiftVector = (1:numel(targim))';
xshiftVector = xshift(:);
yShiftVector = yshift(:) * size(targim, 1);
xyShift = noShiftVector - xshiftVector - yShiftVector;
find(shiftvec~=xyShift)

imageVector = targim(:);
xyShift = xyShift .* (xyShift >= 1) .* (xyShift <= length(noShiftVector)) + noShiftVector .* (xyShift < 1) + noShiftVector .* (xyShift > length(noShiftVector));
correctedImageVector = imageVector(xyShift);
correctedImage = reshape(correctedImageVector, size(targim, 1), size(targim, 2));
%trans = [dxArray' dyArray'];
%}

