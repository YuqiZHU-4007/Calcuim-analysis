function [targvol, xshift, yshift] = frame_reformat_warp_gpu(targvol, params, warpingSettings)
% [correctedImage] = frame_reformat_warp(targvol, params, warpingSettings)
% gpuArray in, gpuArray out

[height, width, depth] = size(targvol);

%% warp
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
%maximumShift = warpingSettings(3);

for dd = 1:depth
    [xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));
%     dxArray = squeeze(params(dd,:,1));
%     dyArray = squeeze(params(dd,:,2));
    %% reformat
    xshift = gpuArray(zeros(size(targvol,1),size(targvol,2)));
    yshift = gpuArray(zeros(size(targvol,1),size(targvol,2)));
    
    movedtileind = find(any(squeeze(params(dd,:,:)), 2));
    for q = movedtileind(:)'
        %for q = 1:numel(xx)
        jj = xx(q); ii = yy(q);
        iRange = max(1, ii - squareSize):min(height, ii + squareSize);
        jRange = max(1, jj - squareSize):min(width, jj + squareSize);
        
        xMaxshiftSquare = xshift(iRange, jRange);
        yMaxshiftSquare = yshift(iRange, jRange);
        
        xMaxshiftSquare(abs(xMaxshiftSquare) < abs(params(dd,q,1))) = params(dd,q,1);
        xshift(iRange, jRange) = xMaxshiftSquare;
        
        yMaxshiftSquare(abs(yMaxshiftSquare) < abs(params(dd,q,2))) = params(dd,q,2);
        yshift(iRange, jRange) = yMaxshiftSquare;
    end
    
    movedpixind = find(xshift ~=0 | yshift ~=0);
    %[y0, x0] = find(xshift ~=0 | yshift ~=0);
    [y0, x0] = ind2sub(size(squeeze(targvol(:,:,dd))), movedpixind);
    x1 = x0 - xshift(movedpixind);
    y1 = y0 - yshift(movedpixind);
    
    validind = find(~(x1<1 | x1>width | y1<1 | y1>height));
    x0 = x0(validind); y0 = y0(validind);
    x1 = x1(validind); y1 = y1(validind);
    
    i0 = sub2ind(size(targvol), y0, x0, dd*ones(length(x0),1));
    i1 = sub2ind(size(targvol), y1, x1, dd*ones(length(x1),1));
    
    %correctedImage = targvol;
    targvol(i0) = targvol(i1);
end
%{
[x0, y0] = meshgrid(1:width, 1:height);
x1 = x0 - xshift;
y1 = y0 - yshift;
i0 = sub2ind(size(targvol), y0, x0);
i1 = sub2ind(size(targvol), y1, x1);
correctedImage = targvol;
correctedImage(i0) = targvol(i1);

shiftvec = (1:numel(targvol))';
shiftvec(i0) = i1;
%}

%{
noShiftVector = (1:numel(targvol))';
xshiftVector = xshift(:);
yShiftVector = yshift(:) * size(targvol, 1);
xyShift = noShiftVector - xshiftVector - yShiftVector;
find(shiftvec~=xyShift)

imageVector = targvol(:);
xyShift = xyShift .* (xyShift >= 1) .* (xyShift <= length(noShiftVector)) + noShiftVector .* (xyShift < 1) + noShiftVector .* (xyShift > length(noShiftVector));
correctedImageVector = imageVector(xyShift);
correctedImage = reshape(correctedImageVector, size(targvol, 1), size(targvol, 2));
%trans = [dxArray' dyArray'];
%}

