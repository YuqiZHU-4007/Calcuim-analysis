function pre_regi_gpu = prepare_reference_image_for_registration_gpu(fxvol, warpingSettings)
% pre_regi = prepare_reference_image_for_registration(fxim,
% warpingSettings): outputing gpuArray

pre_regi_gpu = cell(1,2);
[height, width, depth] = size(fxvol);

%% xlat
% I1 = double(fxim);
I1 = fxvol - mean(fxvol(:)); %%already gpuArray
% I1_gpu = gpuArray(I1);
pre_regi_gpu{1} = conj(fft2(I1));  %%pre_xlat
clear I1

%% warp
pointDensity = warpingSettings(1);
squareSize   = warpingSettings(2);
% maximumShift = warpingSettings(3);

[xx, yy] = meshgrid((squareSize + 1):pointDensity:(width - squareSize), (squareSize + 1):pointDensity:(height - squareSize));

%pre_warp = cell(size(xx));
pre_warp = single(zeros(squareSize*2+1,  squareSize*2+1, depth, numel(xx)));
pre_regi_gpu{2} = gpuArray(pre_warp);
clear pre_warp_temp
% pre_warp = fxim;

for q = 1:numel(xx)  
        jj = xx(q); ii = yy(q);
        %tile2 = targim((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize));
%         tile1 = fxvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), :);
        pre_regi_gpu{2} (:,:,:,q) = fft2(fxvol((ii - squareSize):(ii + squareSize), (jj - squareSize):(jj + squareSize), :));        
end

% pre_regi = cat(3, pre_xlat, pre_warp);
% pre_regi_gpu = {pre_xlat, pre_warp};

end