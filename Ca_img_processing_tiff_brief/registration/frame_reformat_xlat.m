function [targim] = frame_reformat_xlat(targim, params)
% [targim] = frame_reformat(targim, params)
%

%% xlat
dx = params(1); dy = params(2);
if dy >= 0 && dx >= 0
    targim(1:(end - dy), 1:(end - dx)) = targim((1 + dy):end, (1 + dx):end);
elseif dy >= 0 && dx <= 0
    targim(1:(end - dy), (1 - dx):end) = targim((1 + dy):end, 1:(end + dx));
elseif dy <= 0 && dx >= 0
    targim((1 - dy):end, 1:(end - dx)) = targim(1:(end + dy), (1 + dx):end);
elseif dy <= 0 && dx <= 0
    targim((1 - dy):end, (1 - dx):end) = targim(1:(end + dy), 1:(end + dx));
end

