function [indmat, ind] = veccolon(vec, stp, len, totallen)
% [indmat, ind] = veccolon(vec, stp, len, totallen = Inf)
%
% example
% vec = 1:3, stp = 2, len = 3
% indmat = [1 3 5 
%           2 4 6  
%           3 5 7]
% ind = sort(unique(indmat(:))) within totallen
%

if nargin < 4
    totallen = Inf;
end

indmat = repmat(vec(:), [1 len]) + repmat(0:stp:stp*len-1, [length(vec), 1]);
ind = sort(unique(indmat(:)));
ind = ind(ind <= totallen);


