function pre_regi = prepare_reference_image_for_registration(fxim, warpingSettings)
% pre_regi = prepare_reference_image_for_registration(fxim, warpingSettings)

%% xlat
I1 = double(fxim);
I1 = I1 - mean(I1(:));
pre_xlat = conj(fft2(I1));

%% warp
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
% maximumShift = warpingSettings(3);

[height, width] = size(fxim);

[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));

%pre_warp = cell(size(xx));
pre_warp = zeros(squareSize*2+1,  squareSize*2+1, numel(xx));
% pre_warp = fxim;
for q = 1:numel(xx)
        jj = xx(q); ii = yy(q);
        %tile2 = targim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        tile1 = fxim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
        temp = fft2(tile1);
        pre_warp(:,:,q) = temp;
end

% pre_regi = cat(3, pre_xlat, pre_warp);
pre_regi = {pre_xlat, pre_warp};

