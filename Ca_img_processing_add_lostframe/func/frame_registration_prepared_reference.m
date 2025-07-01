function [mvreg, trans] = frame_registration_prepared_reference(pre_regi, targim, warpingSettings)
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
%

xlat  = frame_register_xlat_prepared_reference_backup_afterward(pre_regi{1}, targim);
imaff = frame_reformat_xlat(targim, xlat);
trans = frame_register_warp_prepared_reference_try_to_reverse(pre_regi{2}, imaff, warpingSettings);
mvreg = frame_reformat_warp(imaff, trans, warpingSettings);
trans = [xlat; trans];
end




