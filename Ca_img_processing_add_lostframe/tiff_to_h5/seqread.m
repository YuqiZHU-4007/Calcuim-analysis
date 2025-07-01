function tstack = seqread(imgpath)
% tstack = seqread(imgpath)

[frame, nframes] = frameread(imgpath, 1);
[ysize, xsize] = size(frame);

tstack = zeros(ysize, xsize, nframes, class(frame));
for fi = 1:nframes
    frame = frameread(imgpath, fi);
    tstack(:,:,fi) = frame;
end

