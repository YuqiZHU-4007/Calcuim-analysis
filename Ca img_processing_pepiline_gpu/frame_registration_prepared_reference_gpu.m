function [mvreg, trans] = frame_registration_prepared_reference_gpu(pre_regi, targvol, warpingSettings)
% [mvreg, trans] = frame_registration_prepared_reference(prepared_refrim, targim, warpingSettings)
% 
% warpingSettings:
%   pointDensity = warpingSettings(1);
%   squareSize   = warpingSettings(2);
%   maximumShift = warpingSettings(3);
%
% unit: pixels 
%
% my suggestion: squareSize = pointDensity + maximumShift
%   e.g. [80 128 12]

% %%gpuArray in, gpuArray out

xlat  = frame_register_xlat_prepared_reference_gpu(pre_regi{1}, targvol);
volaff = frame_reformat_xlat_gpu(targvol, xlat);
trans = frame_register_warp1_prepared_reference_gpu(pre_regi{2}, volaff, warpingSettings);
mvreg = frame_reformat_warp_gpu(volaff, trans, warpingSettings);
trans = cat(2, nan(size(trans,1),1,size(trans,3)), trans);
trans(:,1,:) = xlat;
end




