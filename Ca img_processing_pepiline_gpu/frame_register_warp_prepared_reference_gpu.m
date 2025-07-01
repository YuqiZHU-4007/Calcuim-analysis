function [trans] = frame_register_warp_prepared_reference_gpu(pre_warp1, pre_warp2, targvol, warpingSettings)
% [trans] = frame_register_warp_prepared_reference(pre_warp, targim, warpingSettings)
%

%% registration
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
maximumShift = warpingSettings(3);

[height, width, depth] = size(targvol);

%jVector = (squareSize + 1):pointDensity:(size(image, 1) - squareSize);  % point coord along y
%kVector = (squareSize + 1):pointDensity:(size(image, 2) - squareSize);  % point coord along x
[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));

% target vars

%q = 0;
%for j = jVector
%    for k = kVector
%        q = q + 1;
trans = gpuArray(zeros(depth, numel(xx), 2));

for dd = 1:depth    
    dxArray = gpuArray(zeros(numel(xx),1));  % length(jVector) * length(kVector) = number of control points
    dyArray = gpuArray(zeros(numel(xx),1));
    
for q = 1:ceil(numel(xx)/2)
        jj = xx(q); ii = yy(q); 
        tile2 = targvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), dd);
%         I1 = refrim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        tile1fft_gpu = pre_warp1(:,:,dd, q);
        %I1 = imageSquareReference;
        %I2 = imageSquare;
        % todo: caution: why not center the image?
%         C = ifftshift(ifft2(fft2(I1) .* conj(fft2(I2))));
%         tile2_gpu = gpuArray(tile2);
        tile2_fft = fft2(tile2);
        C_temp = ifft2(tile1fft_gpu .* conj(tile2_fft));
        C = ifftshift(C_temp);
%         C = ifftshift(ifft2(tile1fft .* conj(fft2(tile2))));
        %tic; C = normxcorr2(I1, I2); toc;  % this is way slow
        
        [dy, dx] = find(C == max(C(:)));
        dx = dx - squareSize - 2;
        dy = dy - squareSize - 2;
        
        if length(dx) ~= 1
            dx = 0;
            dy = 0;
        end;
        dxArray(q) = dx;
        dyArray(q) = dy;
end

for q = ceil(numel(xx)/2)+1:numel(xx)
        jj = xx(q); ii = yy(q); 
        tile2 = targvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), dd);
%         I1 = refrim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        tile1fft_gpu = pre_warp2(:,:,dd, ceil(numel(xx)/2)+q);
        %I1 = imageSquareReference;
        %I2 = imageSquare;
        % todo: caution: why not center the image?
%         C = ifftshift(ifft2(fft2(I1) .* conj(fft2(I2))));
%         tile2_gpu = gpuArray(tile2);
        tile2_fft = fft2(tile2);
        C_temp = ifft2(tile1fft_gpu .* conj(tile2_fft));
        C = ifftshift(C_temp);
%         C = ifftshift(ifft2(tile1fft .* conj(fft2(tile2))));
        %tic; C = normxcorr2(I1, I2); toc;  % this is way slow
        
        [dy, dx] = find(C == max(C(:)));
        dx = dx - squareSize - 2;
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
trans (dd, :, :)= [dxArray(:) dyArray(:)];
end %%depth
