function [h, incre] = sepplot(sig, prct)
% h = sepplot(sig, prct = 95)

if nargin < 2
    prct = 99;
end

sig = sig';

[nsig, nt] = size(sig);

interval = prctile(sig(:), prct);
incre = repmat((0:interval:interval*nsig-interval)', [1, nt]);
sig = sig + repmat((0:interval:interval*nsig-interval)', [1, nt]);
h = plot(sig');
