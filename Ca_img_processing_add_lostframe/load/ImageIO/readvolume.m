function vol = readvolume(imio, ti, subsamplez)
% vol = readvolume(imio, ti, subsamplez = 1:nz)
%

if nargin < 3
    subsamplez = 1:imio.nz;
end

frame = readframe(imio, ti, subsamplez(1));
[height, width] = size(readframe(imio, 1, 1));
vol = zeros(height, width, numel(subsamplez), 'uint16');
vol(:,:,1) = frame;

for ii = 2:numel(subsamplez)
    zi = subsamplez(ii);
    frame = readframe(imio, ti, zi);
    vol(:,:,ii) = frame;
end
