function [imscale] = rescalegd(im, tails)
% [imscale] = rescalegd(im, tails = [percentmin, percentmax] = [0, 0])
%get an image and scales it between 0 and 1

if nargin < 2
    tails = [0, 0];
end

im = double(im);
ss = sort(im(:));
imin = floor(tails(1) * numel(im)) + 1;
imax = ceil((1-tails(2)) * numel(im));

imscale = mat2gray(im, [ss(imin) ss(imax)]);%Convert matrix to grayscale image
%The returned matrix I contains values in the range 0.0 (black) to 1.0 (full intensity or white). 
%amin and amax are the values in A that correspond to 0.0 and 1.0 in I. Values less than amin become 0.0, and values greater than amax become 1.0.

