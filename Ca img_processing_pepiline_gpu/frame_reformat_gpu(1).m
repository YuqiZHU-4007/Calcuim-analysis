function [correctedImage] = frame_reformat_gpu(targim, params, warpingSettings)
% [correctedImage] = frame_reformat(targim, params, warpingSettings)
%

imaff = frame_reformat_xlat_gpu(targim, params(:,1,:));
correctedImage = frame_reformat_warp_gpu(imaff, params(:,2:end,:), warpingSettings);


