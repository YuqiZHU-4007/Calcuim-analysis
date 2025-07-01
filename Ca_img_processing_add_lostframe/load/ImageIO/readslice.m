function tstack = readslice(imio, zi, subsamplet)
% tstack = readslice(imio, zi, subsamplet)
%

frame = readframe(imio, subsamplet(1), zi);
[height, width] = size(readframe(imio, 1, 1));
tstack = zeros(height, width, numel(subsamplet), 'uint16');
tstack(:,:,1) = frame;

for ii = 2:numel(subsamplet)
    ti = subsamplet(ii);
    frame = readframe(imio, ti, zi);
    tstack(:,:,ii) = frame;
end


