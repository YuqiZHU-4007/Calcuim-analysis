function [trans] = frame_register_warp_prepared_reference_try_to_reverse(pre_warp, targim, warpingSettings)
% [trans] = frame_register_warp(refrim, targim, warpingSettings)
%

%% registration
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
maximumShift = warpingSettings(3);

[height, width] = size(targim);

%jVector = (squareSize + 1):pointDensity:(size(image, 1) - squareSize);  % point coord along y
%kVector = (squareSize + 1):pointDensity:(size(image, 2) - squareSize);  % point coord along x
[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));

% target vars
dxArray = zeros(size(xx));  % length(jVector) * length(kVector) = number of control points
dyArray = zeros(size(xx));

%q = 0;
%for j = jVector
%    for k = kVector
%        q = q + 1;
for q = 1:numel(xx)
        jj = xx(q); ii = yy(q);
        I2 = targim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        %I1 = refrim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        %I1 = imageSquareReference;
        %I2 = imageSquare;
        % todo: caution: why not center the image?
        tile1fft = pre_warp(:,:,q);
		C = ifftshift(ifft2(tile1fft .* conj(fft2(I2))));
        %tic; C = normxcorr2(I1, I2); toc;  % this is way slow
        [dy, dx] = find(C == max(C(:)));
        dx = dx - squareSize - 2; %%for xlat, size are even and - 1024 -1, here size are  2*squareSize+1, thus -round((2*squareSize+1)/2-.5) -1, the same
        dy = dy - squareSize - 2;
        
        if length(dx) ~= 1
            dx = 0;
            dy = 0;
        end;
        dxArray(q) = dx;
        dyArray(q) = dy;
end
%    end;
%end;

% todo: now you can add adaptive grid size

%dxArray = dxArray .* (abs(dxArray) < maximumShift);
%dyArray = dyArray .* (abs(dyArray) < maximumShift);
dxArray(abs(dxArray) > maximumShift) = 0;
dyArray(abs(dyArray) > maximumShift) = 0;
trans = [dxArray(:) dyArray(:)];

