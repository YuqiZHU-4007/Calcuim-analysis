function [h, incre] = myplot_all(t, y)
sig =y;
[nsig, nt] = size(sig);
interval = prctile(sig(:), 99.99);
incre = repmat((0:interval:interval*nsig-interval)', [1, nt]);
sig = sig + repmat((0:interval:interval*nsig-interval)', [1, nt]);
h = plot(t, sig,'b');