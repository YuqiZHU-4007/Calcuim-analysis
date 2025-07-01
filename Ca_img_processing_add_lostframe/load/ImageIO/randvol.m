function vol = randvol(imio, nsamples, subsamplez)
% vol = randvol(imio, nsamples = nt, subsamplez = 1:nz)
% 

if nargin < 3
    subsamplez = 1:imio.nz;
end

if nargin < 2
    nsamples = min([128 imio.nt]);
end

if nsamples > imio.nt
    error('');
end

[height, width] = size(readframe(imio, 1, 1));
depth = imio.nz;
nvol  = imio.nt;

vol = zeros(height, width, depth);
randtimes = sort(randsample(1:nvol, nsamples));
for zi = subsamplez
    zi
    tstack = zeros(height, width, nsamples, 'uint16');
    ii = 1;
    for ti = randtimes
        frame = readframe(imio, ti, zi);
        tstack(:,:,ii) = frame;
        ii = ii + 1;
    end
    tstack = reshape(tstack, [height*width nsamples]);
    tstack = sort(tstack, 2);
    ms = mean(tstack(:, ceil(nsamples/8):end-floor(nsamples/8)), 2);
    vol(:,:,zi) = reshape(ms, [height, width]);
end
vol = uint16(round(vol));

end
