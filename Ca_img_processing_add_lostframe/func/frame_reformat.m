function [correctedImage] = frame_reformat(targim, params, warpingSettings)
% [correctedImage] = frame_reformat(targim, params, warpingSettings)
%

imaff = frame_reformat_xlat(targim, params(1,:));
correctedImage = frame_reformat_warp(imaff, params(2:end,:), warpingSettings);


