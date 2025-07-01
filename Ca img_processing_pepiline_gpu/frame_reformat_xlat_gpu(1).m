function [targvol] = frame_reformat_xlat_gpu(targvol, params)
% [targim] = frame_reformat(targim, params)
%

%% xlat
% dx = params(:,1); dy = params(:,2);
for dd = 1:size(targvol,3)
if params(dd,2) > 0 && params(dd,1) > 0
    targvol(1:(end - params(dd,2)), 1:(end - params(dd,1)), dd) = targvol((1 + params(dd,2)):end, (1 + params(dd,1)):end, dd);
elseif params(dd,2) > 0 && params(dd,1) < 0
    targvol(1:(end - params(dd,2)), (1 - params(dd,1)):end, dd) = targvol((1 + params(dd,2)):end, 1:(end + params(dd,1)), dd);
elseif params(dd,2) < 0 && params(dd,1) > 0
    targvol((1 - params(dd,2)):end, 1:(end - params(dd,1)), dd) = targvol(1:(end + params(dd,2)), (1 + params(dd,1)):end, dd);
elseif params(dd,2) < 0 && params(dd,1) < 0
    targvol((1 - params(dd,2)):end, (1 - params(dd,1)):end, dd) = targvol(1:(end + params(dd,2)), 1:(end + params(dd,1)), dd);
end
end

