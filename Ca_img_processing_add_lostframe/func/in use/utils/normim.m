function [imscale] = normim(im, tails)
% [imscale] = normim(im, tails = [percentmin, percentmax] = [0, 0])
% get an image and scales it between 0 and 1
% a alias of normim

if nargin < 2
    tails = [0, 0];
end

im=double(im);
ss = sort(im(:));
imin = floor(tails(1) * numel(im)) + 1;
imax = ceil((1-tails(2)) * numel(im));

imscale = mat2gray(im, [ss(imin) ss(imax)]);

